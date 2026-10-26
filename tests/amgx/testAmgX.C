// SPDX-License-Identifier: GPL-3.0-or-later
#include "fvCFD.H"
#include <cmath>
#include <memory>

int main(int argc, char* argv[])
{
    argList::noParallel();
    argList::addBoolOption("reuseUnchangedMatrix", "Skip coefficient upload/setup only for an exactly unchanged current matrix");
    argList::addBoolOption("reuseHierarchy", "Refresh numeric multigrid operators with existing structure");
    argList::addOption("config", "file", "AmgX FGMRES configuration JSON");
    argList::addOption("momentumConfig", "file", "Separate solver JSON for simultaneous momentum states");
    argList::addBoolOption("simultaneousFields", "Retain pressure and three momentum states with mixed modes across overlapping meshes");
    argList::addBoolOption("compact", "Use five updates and one overlapping-mesh cycle for instrumented checks");
    argList::addBoolOption("mixed", "Also check dDFI if supported by AmgX");
    #include "setRootCase.H"
    #include "createTime.H"
    const fileName config(args.get<fileName>("config"));
    const bool mixed = args.found("mixed");

    const bool simultaneous = args.found("simultaneousFields");
    const bool compact = args.found("compact");
    const label cycles = simultaneous && compact ? 1 : 2;
    const label updates = compact ? 5 : 10;
    const fileName momentumConfig = args.found("momentumConfig")
        ? args.get<fileName>("momentumConfig") : config;
    auto verify = [&](fvMesh& mesh, const label cycle)
    {
        lduMatrix matrix(mesh);
        FieldField<Field, scalar> boundaryCoeffs(mesh.boundary().size());
        lduInterfaceFieldPtrsList interfaces(mesh.boundary().size());
        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled())
                FatalErrorInFunction << "Verification requires uncoupled patches"
                    << exit(FatalError);
            boundaryCoeffs.set
            (patchi, new scalarField(mesh.boundary()[patchi].size(), Zero));
        }
        const labelUList& upperAddr = mesh.lduAddr().upperAddr();
        const labelUList& lowerAddr = mesh.lduAddr().lowerAddr();
        const wordList modes(mixed ? wordList({"dDDI", "dFFI", "dDFI", "dDDI"})
                                   : wordList({"dDDI", "dFFI", "dDDI"}));
        const wordList fields
        (
            simultaneous
          ? wordList({"verificationPressure", "verificationUx", "verificationUy", "verificationUz"})
          : wordList({"verification"})
        );
        forAll(modes, pass)
        {
            forAll(fields, fieldi)
            {
                const word selectedMode = simultaneous && fieldi == 0
                    ? word("dDDI") : modes[pass];
                // Nonsymmetric, strictly diagonally dominant graph operator.
                // Change the coefficients and sign as well as precision, retaining
                // the same field name to exercise cached-state invalidation.
                forAll(matrix.diag(),i) matrix.diag()[i] = 1.0 + 0.1*pass;
                forAll(matrix.upper(),i) matrix.upper()[i] = -(0.2 + 0.03*pass);
                forAll(matrix.lower(),i) matrix.lower()[i] = -(0.4 + 0.02*pass);
                forAll(upperAddr, f)
                {
                    matrix.diag()[upperAddr[f]] -= matrix.lower()[f];
                    matrix.diag()[lowerAddr[f]] -= matrix.upper()[f];
                }
                if (pass % 2)
                {
                    forAll(matrix.diag(),i) matrix.diag()[i] *= -1;
                    forAll(matrix.upper(),i) matrix.upper()[i] *= -1;
                    forAll(matrix.lower(),i) matrix.lower()[i] *= -1;
                }
                scalarField exact(mesh.nCells()), source(mesh.nCells());
                forAll(exact, i) exact[i] = 0.7 + 0.2*std::sin(scalar(i+1));
                forAll(source,i) source[i] = matrix.diag()[i]*exact[i];
                    forAll(upperAddr,f) {
                        source[lowerAddr[f]] += matrix.upper()[f]*exact[upperAddr[f]];
                        source[upperAddr[f]] += matrix.lower()[f]*exact[lowerAddr[f]];
                    }
                dictionary controls;
                controls.add("solver", word("AmgX"));
                controls.add("tolerance", scalar(1e-12));
                controls.add("relTol", scalar(0));
                controls.add("maxIter", label(2000));
                dictionary options;
                options.add("configFile", simultaneous && fieldi > 0 ? momentumConfig : config);
                options.add("mode", selectedMode);
                options.add("verifyCsr", true);
                options.add("reuseUnchangedMatrix", args.found("reuseUnchangedMatrix"));
                options.add
                (
                    "reuseHierarchy", args.found("reuseHierarchy")
                 && (!simultaneous || fieldi == 0 || !args.found("momentumConfig"))
                );
                controls.add("AmgX", options);
                autoPtr<lduMatrix::solver> solver = lduMatrix::solver::New
                (
                    fields[fieldi], matrix, boundaryCoeffs, boundaryCoeffs,
                    interfaces, controls
                );
                // More than three solves per state exercise value updates after
                // the automatic initial CSR operator probes have finished.
                for (label update = 0; update < updates; ++update)
                {
                    if (update % 2 == 0) { forAll(matrix.diag(),i) matrix.diag()[i] += 0.01; }
                    else { forAll(exact,i) exact[i] += 0.03*std::cos(scalar(i+2)); }
                    forAll(source,i) source[i] = matrix.diag()[i]*exact[i];
                    forAll(upperAddr,f) {
                        source[lowerAddr[f]] += matrix.upper()[f]*exact[upperAddr[f]];
                        source[upperAddr[f]] += matrix.lower()[f]*exact[lowerAddr[f]];
                    }
                    scalarField solution(mesh.nCells(), Zero);
                    const solverPerformance perf = solver->solve(solution, source);
                    scalar error=0; forAll(solution,i) error=std::max(error,std::abs(solution[i]-exact[i]));
                    if (perf.finalResidual() > 1e-12 || !(error < 1e-8))
                        FatalErrorInFunction << "AmgX manufactured-solution error "
                            << error << exit(FatalError);
                    Info<< "AMGX_VERIFICATION cycle=" << cycle << " mode=" << selectedMode
                        << " field=" << fields[fieldi]
                        << " update=" << update << " error=" << error << endl;
                }
            }
        }
    };
    // Keep populated caches on two meshes alive together. Destroy the second,
    // then solve again on the first to expose premature global-pool cleanup.
    for (label cycle = 0; cycle < cycles; ++cycle)
    {
        fvMesh mesh
        (
            IOobject
            (
                polyMesh::defaultRegion, runTime.timeName(), runTime,
                IOobject::MUST_READ
            )
        );
        verify(mesh, cycle);
        if (simultaneous)
        {
            {
                fvMesh second
                (
                    IOobject
                    (
                        polyMesh::defaultRegion, runTime.timeName(), runTime,
                        IOobject::MUST_READ, IOobject::NO_WRITE, false
                    )
                );
                verify(second, cycle + 2);
            }
            verify(mesh, cycle + 4);
        }
    }
    Info<< "AMGX_VERIFICATION passed" << endl;
    return 0;
}
