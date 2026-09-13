/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::explicitFOmembraneVelocityFvPatchVectorField

Description
    Forward-osmosis water flux for an asymmetric membrane in AL-FS orientation.
    Both sides must be coincident, equal-area faces in the same conformal patch.
    Each pair must remain local to one MPI rank. forwardDirection is required
    and points from feed into draw, independently of patch-face ordering.

    Example in 0/U:
    membrane
    {
        type explicitFOmembraneVelocity;
        forwardDirection (0 0 1);
        eq advanced;
        slip noSlip;
        value uniform (0 0 0);
    }

    A [m/(Pa s)], B [m/s] and K [s/m] are read from transportProperties,
    not from this patch dictionary. m_A is a mass fraction [kg/kg].
    eq advanced includes reverse salt permeability; eq simple is the
    high-rejection approximation. The supported osmotic mode assumes negligible
    hydraulic pressure difference and draw surface concentration >= feed.
    Reversed driving or invalid pairing fails explicitly.

    Draw velocity is scaled to conserve water mass across the density jump.
    Optional slip applies only on the draw side: alpha is a length [m], and
    aRelax is its dimensionless relaxation factor. Use noSlip for the paper cases.
    The legacy m_A none mode instead uses A times the pressure difference [Pa].

SourceFiles
    explicitFOmembraneVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef explicitFOmembraneVelocityFvPatchVectorField_H
#define explicitFOmembraneVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class explicitFOmembraneVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class explicitFOmembraneVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Transport Properties dictionary
        IOdictionary transProps_;

        //- Name of the pressure field
        word pName_;

        //- Name of the concentration field
        word m_AName_;

        //- "Forward" direction of membrane
        vector forwardDirection_;

        //- Name of the flux equation
        word fluxEqName_;

        //- A coefficient
        scalar A_;

        //- Solute permeation coefficient
        scalar B_;

        //- K coefficient
        scalar K_;
        
        //- Slip or noSlip
        word slipName_;

        //- Non-dimensioned Porosity Coefficient
        scalar alpha_;

        //- Slip relaxation
        scalar aRelax_;

        //- Osmotic pressure coefficient
        dimensionedScalar pi_mACoeff_;

        //- Solution density
        dimensionedScalar rho0_;

        //- Density m_A coefficient
        dimensionedScalar rho_mACoeff_;

        //- face indexing list
        List<label> fm_;

        //- feed-side index list
        List<label> fs_;

        //- Reuse roots only for exactly unchanged concentrations. Coefficients
        //  are fixed for this patch-field instance; mapping invalidates the cache.
        scalarField cachedFeed_;
        scalarField cachedDraw_;
        scalarField cachedFlux_;

public:

    //- Runtime type information
    TypeName("explicitFOmembraneVelocity");


    // Constructors

        //- Construct from patch and internal field
        explicitFOmembraneVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        explicitFOmembraneVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given explicitFOmembraneVelocityFvPatchVectorField
        //  onto a new patch
        explicitFOmembraneVelocityFvPatchVectorField
        (
            const explicitFOmembraneVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        explicitFOmembraneVelocityFvPatchVectorField
        (
            const explicitFOmembraneVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new explicitFOmembraneVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        explicitFOmembraneVelocityFvPatchVectorField
        (
            const explicitFOmembraneVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new explicitFOmembraneVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the transport properties dictionary
            const IOdictionary& transportProperties() const
            {
                return transProps_;
            }

            //- Return the name of the pressure field
            const word& pName() const
            {
                return pName_;
            }

            //- Return the name of the concentration field
            const word& m_AName() const
            {
                return m_AName_;
            }

            //- Return the "forward" direction of the membrane
            const vector& forwardDirection() const
            {
                return forwardDirection_;
            }

            //- Return the name of the flux equation
            const word& fluxEqName() const
            {
                return fluxEqName_;
            }

            //- Return the value of the A coefficient
            scalar A() const
            {
                return A_;
            }

            //- Return the value of the solute permeation coefficient
            scalar B() const
            {
                return B_;
            }

            //- Return the value of the K coefficient
            scalar K() const
            {
                return K_;
            }
           
            //- Return name of slip condition.
            const word& slipName() const
            {
                return slipName_;
            }

            //- Return the value of the alpha
            scalar alpha() const
            {
                return alpha_;
            }

            //- Return the value of the slip relaxation factor
            scalar aRelax() const
            {
                return aRelax_;
            }

            //- Return the osmotic pressure coefficient
            dimensionedScalar pi_mACoeff() const
            {
                return pi_mACoeff_;
            }

            //- Return the fluid density
            dimensionedScalar rho0() const
            {
                return rho0_;
            }

            //- Return the density m_A coefficient
            dimensionedScalar rho_mACoeff() const
            {
                return rho_mACoeff_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;

private:
        //- initialise the object
        void initialise();

        //- compute the cell-index mapping
        void calcFaceMapping();

        //- Monotone, bracketed membrane flux solve; accumulates iterations.
        scalar solveFlux(const scalar& feed, const scalar& draw, scalar& iterations);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

