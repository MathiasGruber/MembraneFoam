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
    Foam::explicitROmembraneVelocityFvPatchVectorField

Description
    Foam::explicitROmembraneVelocityFvPatchVectorField

    this boundary condition is designed to simulate a membrane 
    across which the velocity is set by the pressure difference calculated 
    from the previous timestep. The concentration can also be set. The full 
    boundary condition set up with typical values is given below:
    {
        type 	explicitROmembraneVelocity;
        value	uniform (0 0 0);// initialisation value
        K	2e-10;		// K coefficient
        m_A	none;		// field of concentration in moles (can be set
                                // to "none" or left out if the osmotic 
                                // pressure is not to taken into account
    }
    Density changes across the membrane are taken into account.

    The values of hydrodynamic pressure do not take gravity explicitly into 
    account, so it is assumed that the membrane has a thickness of zero.
 
SourceFiles
    explicitROmembraneVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef explicitROmembraneVelocityFvPatchVectorField_H
#define explicitROmembraneVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class explicitROmembraneVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class explicitROmembraneVelocityFvPatchVectorField
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

        //- K coefficient
        scalar K_;

        //- Osmotic pressure coefficient
        dimensionedScalar pi_mACoeff_;

        //- Solution density
        dimensionedScalar rho0_;

        //- Density m_A coefficient
        dimensionedScalar rho_mACoeff_;

        //- face indexing list
        List<label> fm_;

public:

    //- Runtime type information
    TypeName("explicitROmembraneVelocity");


    // Constructors

        //- Construct from patch and internal field
        explicitROmembraneVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        explicitROmembraneVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given explicitROmembraneVelocityFvPatchVectorField
        //  onto a new patch
        explicitROmembraneVelocityFvPatchVectorField
        (
            const explicitROmembraneVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        explicitROmembraneVelocityFvPatchVectorField
        (
            const explicitROmembraneVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new explicitROmembraneVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        explicitROmembraneVelocityFvPatchVectorField
        (
            const explicitROmembraneVelocityFvPatchVectorField&,
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
                new explicitROmembraneVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

            virtual bool fixesValue() const
            {
                return true;
            }

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

            //- Return the value of the K coefficient
            scalar K() const
            {
                return K_;
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
    //- initialise the object and compute the cell-index mapping
    void calcFaceMapping();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
