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

#ifndef membraneAmgXSolver_H
#define membraneAmgXSolver_H

#include "lduMatrix.H"
#include "ldu2csr.H"
#include <amgx_c.h>
#include <memory>
#include <string>

namespace Foam
{
class AmgXRuntime;

// Serial, static-mesh AmgX solver with OpenFOAM's double-precision stopping rule.
class AmgXSolver : public lduMatrix::solver
{
public:
    struct State
    {
        // Declared first: the runtime outlives all handles and buffers.
        std::shared_ptr<AmgXRuntime> runtime;
        AMGX_config_handle config = nullptr;
        AMGX_matrix_handle matrix = nullptr;
        AMGX_vector_handle solution = nullptr;
        AMGX_vector_handle rhs = nullptr;
        AMGX_solver_handle solver = nullptr;
        autoPtr<ldu2csr> csr;
        List<float> singleMatrix, singleSolution, singleRhs;
        std::string signature;
        label cells = 0;
        label faces = 0;
        label operatorChecks = 0;
        label setupCalls = 0;
        label resetupCalls = 0;
        label unchangedSkips = 0;
        List<scalar> previousMatrix;
        bool uploaded = false;

        State() = default;
        ~State();
        State(const State&) = delete;
        State& operator=(const State&) = delete;
    };

private:
    word mode_;
    std::string config_;
    bool verifyCsr_;
    bool reuseHierarchy_;
    bool reuseUnchanged_;

    AMGX_Mode mode() const;
    State& state() const;
    void initialise(State&) const;
    void upload(State&) const;
    void checkOperator(State&, const direction) const;
    label correction(State&, scalarField&, const scalarField&) const;

public:
    TypeName("AmgX");

    AmgXSolver
    (
        const word& fieldName,
        const lduMatrix& matrix,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const FieldField<Field, scalar>& interfaceIntCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const dictionary& solverControls
    );

    ~AmgXSolver() override = default;
    AmgXSolver(const AmgXSolver&) = delete;
    AmgXSolver& operator=(const AmgXSolver&) = delete;

    solverPerformance solve
    (
        scalarField& psi,
        const scalarField& source,
        const direction cmpt = 0
    ) const override;
};
}
#endif
