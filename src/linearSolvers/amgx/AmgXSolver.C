/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | openfoam-amgx
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of openfoam-amgx.

    openfoam-amgx is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    openfoam-amgx is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with openfoam-amgx.  If not, see <https://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "AmgXSolver.H"
#include "fvMesh.H"
#include "regIOobject.H"
#include "Time.H"
#include "Pstream.H"
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include "executionPolicy.H"
#include <fstream>
#include <map>
#include <sstream>

#define AMGX_CHECK(call)                                                     \
    do                                                                       \
    {                                                                        \
        const AMGX_RC code = (call);                                         \
        if (code != AMGX_RC_OK)                                              \
            FatalErrorInFunction << "AmgX call failed: " #call               \
                << " (code " << int(code) << ')' << exit(FatalError);        \
    } while (false)

namespace Foam
{
static_assert(sizeof(scalar) == sizeof(double), "AmgX requires double fields");
static_assert(sizeof(label) == sizeof(int), "AmgX requires int32 labels");

defineTypeNameAndDebug(AmgXSolver, 0);
lduMatrix::solver::addsymMatrixConstructorToTable<AmgXSolver>
    addAmgXSymmetricConstructor;
lduMatrix::solver::addasymMatrixConstructorToTable<AmgXSolver>
    addAmgXAsymmetricConstructor;

// The mesh registry destroys all State objects before its runtime is released.
// No CUDA or AmgX handles are retained in process-exit static destructors.
class AmgXRuntime
{
    AMGX_config_handle resourceConfig_ = nullptr;
    AMGX_resources_handle resources_ = nullptr;

public:
    AmgXRuntime() { AMGX_CHECK(AMGX_initialize()); }
    ~AmgXRuntime()
    {
        // Resources destroy process-wide pools and library handles. Every
        // field's solver, vectors and matrix must already have been released.
        bool failed = false;
        if (resources_)
            failed |= AMGX_resources_destroy(resources_) != AMGX_RC_OK;
        if (resourceConfig_)
            failed |= AMGX_config_destroy(resourceConfig_) != AMGX_RC_OK;
        if (failed)
            WarningInFunction << "AmgX resource destruction failed" << endl;
        if (AMGX_finalize() != AMGX_RC_OK)
            WarningInFunction << "AmgX finalisation failed" << endl;
    }

    AMGX_resources_handle resources(const std::string& config)
    {
        if (!resources_)
        {
            // Keep an independent configuration alive for the resource's
            // complete lifetime. Field solver configurations remain separate.
            AMGX_CHECK(AMGX_config_create(&resourceConfig_, config.c_str()));
            AMGX_CHECK(AMGX_resources_create_simple(&resources_, resourceConfig_));
        }
        return resources_;
    }

    static std::shared_ptr<AmgXRuntime> acquire()
    {
        static std::weak_ptr<AmgXRuntime> current;
        auto runtime = current.lock();
        if (!runtime)
        {
            runtime = std::make_shared<AmgXRuntime>();
            current = runtime;
        }
        return runtime;
    }
};

// Keep one library runtime across sequential mesh lifetimes in a Time
// context. Time clears owned registry objects before unloading its libraries.
class AmgXRuntimeOwner : public regIOobject
{
public:
    TypeName("membraneAmgXRuntime");
    std::shared_ptr<AmgXRuntime> runtime;

    explicit AmgXRuntimeOwner(const Time& time)
    : regIOobject
      (
          IOobject
          (
              "membraneAmgXRuntime", time.timeName(), time,
              IOobject::NO_READ, IOobject::NO_WRITE
          )
      ),
      runtime(AmgXRuntime::acquire())
    {}

    bool writeData(Ostream&) const override { return false; }
};
defineTypeNameAndDebug(AmgXRuntimeOwner, 0);

class AmgXMeshCache : public regIOobject
{
public:
    TypeName("membraneAmgXCache");
    std::shared_ptr<AmgXRuntime> runtime;
    std::map<std::string, std::unique_ptr<AmgXSolver::State>> states;

    explicit AmgXMeshCache(const objectRegistry& db)
    : regIOobject
      (
          IOobject
          (
              "membraneAmgXCache", db.time().timeName(), db,
              IOobject::NO_READ, IOobject::NO_WRITE
          )
      ),
      runtime(AmgXRuntime::acquire())
    {}

    bool writeData(Ostream&) const override { return false; }
};
defineTypeNameAndDebug(AmgXMeshCache, 0);

static void finishDownload()
{
    // Downloads to unified/system memory may be queued on the default stream.
    const cudaError_t code = cudaStreamSynchronize(0);
    if (code != cudaSuccess)
        FatalErrorInFunction << cudaGetErrorString(code) << exit(FatalError);
}

AmgXSolver::State::~State()
{
    bool failed = false;
    if (solver) failed |= AMGX_solver_destroy(solver) != AMGX_RC_OK;
    if (solution) failed |= AMGX_vector_destroy(solution) != AMGX_RC_OK;
    if (rhs) failed |= AMGX_vector_destroy(rhs) != AMGX_RC_OK;
    if (matrix) failed |= AMGX_matrix_destroy(matrix) != AMGX_RC_OK;
    if (config) failed |= AMGX_config_destroy(config) != AMGX_RC_OK;
    if (failed) WarningInFunction << "AmgX handle destruction failed" << endl;
}

AmgXSolver::AmgXSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName, matrix, interfaceBouCoeffs, interfaceIntCoeffs,
        interfaces, solverControls
    ),
    mode_("dDDI"),
    config_(),
    verifyCsr_(true),
    reuseHierarchy_(false),
    reuseUnchanged_(false)
{
    const dictionary& options = controlDict_.subDict("AmgX");
    mode_ = options.lookupOrDefault<word>("mode", "dDDI");
    verifyCsr_ = options.lookupOrDefault<bool>("verifyCsr", true);
    reuseHierarchy_ = options.lookupOrDefault<bool>("reuseHierarchy", false);
    reuseUnchanged_ = options.lookupOrDefault<bool>("reuseUnchangedMatrix", false);
    fileName path(options.get<fileName>("configFile"));
    path.expand();
    std::ifstream input(path.c_str());
    if (!input)
        FatalErrorInFunction << "Cannot read AmgX configuration " << path
            << exit(FatalError);
    std::ostringstream contents;
    contents << input.rdbuf();
    config_ = contents.str();
    if (config_.empty())
        FatalErrorInFunction << "Empty AmgX configuration " << path
            << exit(FatalError);
    (void)mode();
    if (minIter_ != 0)
        FatalErrorInFunction << "AmgX requires minIter 0" << exit(FatalError);
}

AMGX_Mode AmgXSolver::mode() const
{
    if (mode_ == "dDDI") return AMGX_mode_dDDI;
    if (mode_ == "dDFI") return AMGX_mode_dDFI;
    if (mode_ == "dFFI") return AMGX_mode_dFFI;
    FatalErrorInFunction << "Supported AmgX modes: dDDI, dDFI, dFFI"
        << exit(FatalError);
    return AMGX_mode_dDDI;
}

AmgXSolver::State& AmgXSolver::state() const
{
    const objectRegistry& db = matrix_.mesh().thisDb();
    const Time& time = db.time();
    const word runtimeName("membraneAmgXRuntime");
    if (!time.foundObject<AmgXRuntimeOwner>(runtimeName, false))
    {
        if (time.found(runtimeName))
            FatalErrorInFunction << "Object registry name collision: " << runtimeName
                << exit(FatalError);
        regIOobject::store(new AmgXRuntimeOwner(time));
    }
    const word name("membraneAmgXCache");
    if (!db.foundObject<AmgXMeshCache>(name, false))
    {
        if (db.found(name))
            FatalErrorInFunction << "Object registry name collision: " << name
                << exit(FatalError);
        regIOobject::store(new AmgXMeshCache(db));
    }
    auto& cache = db.lookupObjectRef<AmgXMeshCache>(name, false);
    auto& entry = cache.states[std::string(fieldName_)];
    const std::string signature = std::string(mode_) + '\n' + config_
        + (reuseHierarchy_ ? "\nreuseNumeric" : "\nfullSetup")
        + (reuseUnchanged_ ? "\nexactMatrixReuse" : "\nrefreshEachSolve");
    if (entry && entry->signature != signature) entry.reset();
    if (!entry)
    {
        entry.reset(new State);
        entry->signature = signature;
        initialise(*entry);
    }
    if
    (
        entry->cells != matrix_.diag().size()
     || entry->faces != matrix_.upper().size()
    )
        FatalErrorInFunction << "AmgX requires a static mesh" << exit(FatalError);
    return *entry;
}

void AmgXSolver::initialise(State& st) const
{
    st.runtime = AmgXRuntime::acquire();
    st.cells = matrix_.diag().size();
    st.faces = matrix_.upper().size();
    st.csr.reset(new ldu2csr(matrix_));
    AMGX_CHECK(AMGX_config_create(&st.config, config_.c_str()));
    const AMGX_resources_handle resources = st.runtime->resources(config_);
    AMGX_CHECK(AMGX_matrix_create(&st.matrix, resources, mode()));
    AMGX_CHECK(AMGX_vector_create(&st.solution, resources, mode()));
    AMGX_CHECK(AMGX_vector_create(&st.rhs, resources, mode()));
    AMGX_CHECK(AMGX_solver_create(&st.solver, resources, mode(), st.config));
}

void AmgXSolver::checkOperator(State& st, const direction cmpt) const
{
    if (!verifyCsr_ || st.operatorChecks >= 3) return;
    ++st.operatorChecks;
    scalarField probe(st.cells), expected(st.cells);
    forAll(probe, i) probe[i] = scalar((i % 97)*17 % 97)/97;
    matrix_.Amul(expected, probe, interfaceBouCoeffs_, interfaces_, cmpt);
    scalar error = 0, scale = 0;
    forAll(probe, i)
    {
        scalar value = 0;
        for (label j = st.csr->rowPtr()[i]; j < st.csr->rowPtr()[i+1]; ++j)
            value += st.csr->values()[j]*probe[st.csr->colInd()[j]];
        error += mag(value-expected[i]);
        scale += mag(expected[i]);
    }
    Info<< "AMGX_CSR_AUDIT relative_L1=" << error/max(scale, VSMALL) << endl;
    if (!(error <= 1e-12*max(scale, VSMALL)))
        FatalErrorInFunction << "CSR operator differs from OpenFOAM: relative L1 "
            << error/max(scale, VSMALL) << exit(FatalError);
}

void AmgXSolver::upload(State& st) const
{
    const bool refreshExisting = st.uploaded && reuseHierarchy_;
    const label nnz = st.csr->nNonZeros();
    const void* values = st.csr->values();
    if (!std::all_of
        (amgxExecution, st.csr->values(), st.csr->values()+nnz,
         [](scalar value) { return std::isfinite(value); }))
        FatalErrorInFunction << "Nonfinite matrix coefficient" << exit(FatalError);
    // The current double matrix must be exactly equal, with the same cached
    // graph, precision and configuration. No tolerance/hash approximation.
    // RHS upload, zero guess and true OpenFOAM residual checks still run.
    if (reuseUnchanged_ && st.uploaded && st.previousMatrix.size() == nnz
     && std::equal(amgxExecution, st.csr->values(), st.csr->values()+nnz,
                   st.previousMatrix.cbegin()))
    {
        ++st.unchangedSkips;
        if (st.unchangedSkips <= 3)
            Info<< "AMGX_MATRIX_UNCHANGED field=" << fieldName_
                << " skipped_setups=" << st.unchangedSkips << endl;
        return;
    }
    if (mode_ == "dDFI" || mode_ == "dFFI")
    {
        st.singleMatrix.setSize(nnz);
        std::transform
        (
            amgxExecution, st.csr->values(), st.csr->values()+nnz,
            st.singleMatrix.begin(), [](scalar value) { return float(value); }
        );
        if (!std::all_of
            (amgxExecution,
             st.singleMatrix.cbegin(), st.singleMatrix.cend(),
             [](float value) { return std::isfinite(value); }))
            FatalErrorInFunction << "Matrix exceeds float range" << exit(FatalError);
        values = st.singleMatrix.cbegin();
    }
    if (!st.uploaded)
    {
        AMGX_CHECK
        (
            AMGX_matrix_upload_all
            (
                st.matrix, st.cells, nnz, 1, 1,
                st.csr->rowPtr(), st.csr->colInd(), values, nullptr
            )
        );
        st.uploaded = true;
    }
    else
    {
        AMGX_CHECK
        (
            AMGX_matrix_replace_coefficients
            (st.matrix, st.cells, nnz, values, nullptr)
        );
    }
    // Always refresh the fine matrix above. When requested, resetup refreshes
    // numerical coarse operators while configuration controls structure reuse.
    // In particular, setup with structure_reuse_levels=-1 would return early;
    // only the numeric resetup path is used for that experimental setting.
    if (refreshExisting)
    {
        AMGX_CHECK(AMGX_solver_resetup(st.solver, st.matrix));
        ++st.resetupCalls;
    }
    else
    {
        AMGX_CHECK(AMGX_solver_setup(st.solver, st.matrix));
        ++st.setupCalls;
    }
    if (reuseUnchanged_)
    {
        st.previousMatrix.setSize(nnz);
        std::copy(amgxExecution, st.csr->values(), st.csr->values()+nnz,
                  st.previousMatrix.begin());
    }
    if (st.setupCalls+st.resetupCalls <= 3)
        Info<< "AMGX_HIERARCHY_AUDIT field=" << fieldName_
            << " setup=" << st.setupCalls << " resetup=" << st.resetupCalls << endl;
}

label AmgXSolver::correction
(
    State& st, scalarField& psi, const scalarField& rhs
) const
{
    // The zero-guess solve does not size a freshly created solution vector.
    // Allocate it explicitly and reset every correction to a known zero.
    AMGX_CHECK(AMGX_vector_set_zero(st.solution, st.cells, 1));
    if (mode_ == "dFFI")
    {
        const scalar scale = gMax(mag(rhs));
        if (!(scale > 0) || !std::isfinite(scale))
            FatalErrorInFunction << "Invalid residual correction scale"
                << exit(FatalError);
        st.singleRhs.setSize(st.cells);
        st.singleSolution.setSize(st.cells);
        std::transform
        (
            amgxExecution, rhs.cbegin(), rhs.cend(),
            st.singleRhs.begin(),
            [scale](scalar value) { return float(value/scale); }
        );
        AMGX_CHECK(AMGX_vector_upload(st.rhs, st.cells, 1, st.singleRhs.cbegin()));
        AMGX_CHECK(AMGX_solver_solve_with_0_initial_guess
            (st.solver, st.rhs, st.solution));
        AMGX_CHECK(AMGX_vector_download(st.solution, st.singleSolution.begin()));
        finishDownload();
        std::transform
        (
            amgxExecution,
            st.singleSolution.cbegin(), st.singleSolution.cend(), psi.begin(),
            [scale](float value) { return scale*scalar(value); }
        );
    }
    else
    {
        AMGX_CHECK(AMGX_vector_upload(st.rhs, st.cells, 1, rhs.cbegin()));
        AMGX_CHECK(AMGX_solver_solve_with_0_initial_guess
            (st.solver, st.rhs, st.solution));
        AMGX_CHECK(AMGX_vector_download(st.solution, psi.begin()));
        finishDownload();
    }
    AMGX_SOLVE_STATUS status;
    AMGX_CHECK(AMGX_solver_get_status(st.solver, &status));
    // A finite, improving inexact solve may be corrected by the double outer
    // iteration. Numerical breakdown is always an error.
    if (status == AMGX_SOLVE_FAILED || status == AMGX_SOLVE_DIVERGED)
        FatalErrorInFunction << "AmgX correction failed: status " << int(status)
            << exit(FatalError);
    int iterations = 0;
    AMGX_CHECK(AMGX_solver_get_iterations_number(st.solver, &iterations));
    return iterations;
}

solverPerformance AmgXSolver::solve
(
    scalarField& psi, const scalarField& source, const direction cmpt
) const
{
    if (Pstream::parRun())
        FatalErrorInFunction << "AmgX supports one serial process" << exit(FatalError);
    forAll(interfaces_, i)
        if (interfaces_.set(i))
            FatalErrorInFunction << "AmgX does not support coupled interfaces"
                << exit(FatalError);
    const auto& mesh = refCast<const fvMesh>(matrix_.mesh().thisDb());
    if (mesh.changing())
        FatalErrorInFunction << "AmgX requires a static mesh" << exit(FatalError);
    if (psi.empty() || source.size() != psi.size())
        FatalErrorInFunction << "Invalid AmgX field sizes" << exit(FatalError);

    solverPerformance performance(typeName, fieldName_);
    scalarField Apsi(psi.size()), work(psi.size());
    matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);
    const scalar norm = this->normFactor(psi, source, Apsi, work);
    scalarField residual(source-Apsi);
    performance.initialResidual() =
        gSumMag(residual, matrix_.mesh().comm())/norm;
    performance.finalResidual() = performance.initialResidual();
    const scalar target = max(tolerance_, relTol_*performance.initialResidual());
    if (!(norm > 0) || !std::isfinite(norm)
        || !std::isfinite(performance.initialResidual()) || !std::isfinite(target))
        FatalErrorInFunction << "Nonfinite OpenFOAM residual" << exit(FatalError);
    if (performance.initialResidual() <= target)
    {
        performance.checkConvergence(tolerance_, relTol_);
        return performance;
    }

    State& st = state();
    st.csr->updateValues(matrix_);
    checkOperator(st, cmpt);
    const bool flipped = matrix_.diag()[0] < 0;
    if (flipped)
        std::transform
        (
            amgxExecution,
            st.csr->values(), st.csr->values()+st.csr->nNonZeros(),
            st.csr->values(), [](scalar value) { return -value; }
        );
    upload(st);

    scalar previous = performance.finalResidual();
    label stagnant = 0;
    for (label outer = 0; outer < 12 && performance.nIterations() < maxIter_; ++outer)
    {
        scalarField delta(psi.size(), Zero), rhs(residual);
        if (flipped)
            std::transform
            (
                amgxExecution, rhs.cbegin(), rhs.cend(), rhs.begin(),
                [](scalar value) { return -value; }
            );
        performance.nIterations() += correction(st, delta, rhs);
        psi += delta;
        matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);
        residual = source-Apsi;
        performance.finalResidual() =
            gSumMag(residual, matrix_.mesh().comm())/norm;
        if (!std::isfinite(performance.finalResidual()))
            FatalErrorInFunction << "Nonfinite true double residual" << exit(FatalError);
        if (performance.finalResidual() <= target) break;
        stagnant = performance.finalResidual() > 0.99*previous ? stagnant+1 : 0;
        if (stagnant >= 2) break;
        previous = performance.finalResidual();
    }
    if (performance.finalResidual() > target || performance.nIterations() >= maxIter_)
        FatalErrorInFunction << "AmgX did not meet the original OpenFOAM target: "
            << performance.finalResidual() << " > " << target
            << "; iterations " << performance.nIterations() << exit(FatalError);
    performance.checkConvergence(tolerance_, relTol_);
    return performance;
}
}
