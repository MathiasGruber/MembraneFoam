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
    Foam::explicitFOmembraneSoluteFvPatchScalarField

Description
    Foam::explicitFOmembraneSoluteFvPatchScalarField

    This boundary condition is designed to simulate an assymetric membrane 
    across which a scalar is transported. This is designed to work with the
    explicitFOmembraneSoluteVelocityFvPatchVectorField. The full boundary 
    condition set up with typical values is given below:
    {
        type 	    explicitFOmembraneSolute;
        value	    uniform 0;
        U	U; 		                        // velocity field (default value of "U")
        B	1e-7; 		                    // solute permeation coefficient 
    }

    The following parameters are set in the "constant/transportProperties"
    dictionary and are shared with the solver.

    pi_mACoeff	    pi_mACoeff	    [ 0 0 0 0 0 0 0 ]	    805.1e5;
    mu0		        mu0 		    [ 1 -1 -1 0 0 0 0 ]	    0.89e-3;
    mu_mACoeff	    mu_mACoeff 	    [ 1 -1 -1 0 0 0 0 ] 	1.63; 
    D_AB_Min	    D_AB_Min 	    [ 0 2 -1 0 0 0 0 ]	    1.45e-9;
    D_AB_Coeff	    D_AB_Coeff 	    [ 0 2 -1 0 0 0 0 ]	    1.61e-9;
    D_AB_mACoeff	D_AB_mACoeff 	[ 0 0 0 0 0 0 0 ]	    14; 
    rho0		    rho0		    [ 1 -3 0 0 0 0 0 ]	    997.1;
    rho_mACoeff	    rho_mACoeff	    [ 0 0 0 0 0 0 0 ]	    0.696;
 
SourceFiles
    explicitFOmembraneSoluteFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef explicitFOmembraneSoluteFvPatchScalarField_H
#define explicitFOmembraneSoluteFvPatchScalarField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class explicitFOmembraneSoluteFvScalarPatch Declaration
\*---------------------------------------------------------------------------*/

class explicitFOmembraneSoluteFvPatchScalarField
:
      public fvPatchField<scalar>
{
    // Private data

        //- Transport Properties dictionary
        IOdictionary transProps_;

        //- Name of the velocity field
        word UName_;

        //- Pure water permeability
        scalar A_;

        //- Rejection coefficient
        scalar B_;
    
        // Diffusion resistivity
        scalar K_;

        //- Diffusivity coefficient
        dimensionedScalar D_AB_Min_; 

        //- Diffusivity coefficient
        dimensionedScalar D_AB_Coeff_; 

        //- Diffusivity coefficient
        dimensionedScalar D_AB_mACoeff_; 

        //- Density of fluid
        dimensionedScalar rho0_;

        //- Density m_A coefficient
        dimensionedScalar rho_mACoeff_;

        //- Osmotic pressure coefficient
        dimensionedScalar pi_mACoeff_;

        //- Face indexing list
        List<label> fm_;

        //- valueInternalCoefficients
        Field<scalar> VIC_;

        //- valueBoundaryCoefficients
        Field<scalar> VBC_;

        //- gradientInternalCoefficients
        Field<scalar> GIC_;

        //- gradientBoundaryCoefficients
        Field<scalar> GBC_;

        //- current Time index
        label curTimeIndex_;

public:

    //- Runtime type information
    TypeName("explicitFOmembraneSolute");


    // Constructors

        //- Construct from patch and internal field
        explicitFOmembraneSoluteFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        explicitFOmembraneSoluteFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given explicitFOmembraneSoluteFvPatchScalarField
        //  onto a new patch
        explicitFOmembraneSoluteFvPatchScalarField
        (
            const explicitFOmembraneSoluteFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        explicitFOmembraneSoluteFvPatchScalarField
        (
            const explicitFOmembraneSoluteFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new explicitFOmembraneSoluteFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        explicitFOmembraneSoluteFvPatchScalarField
        (
            const explicitFOmembraneSoluteFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new explicitFOmembraneSoluteFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the Transport Properties dictionary
            const IOdictionary& transportProperties() const
            {
                return transProps_;
            }

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return the rejection coefficient
            scalar B() const
            {
                return B_;
            }

            //- Return the diffusion resistivity
            scalar K() const
            {
                return K_;
            }

            //- Return the minimum diffusivity
            dimensionedScalar D_AB_Min() const
            {
                return D_AB_Min_;
            }

            //- Return the main diffusivity coefficient
            dimensionedScalar D_AB_Coeff() const
            {
                return D_AB_Coeff_;
            }

            //- Return the m_A diffusivity coefficient
            dimensionedScalar D_AB_mACoeff() const
            {
                return D_AB_mACoeff_;
            }

            //- Return the solution density
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
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar> > gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;

    // Member operators
/*
        virtual void operator=(const UList<scalar>&) {}

        virtual void operator=(const fvPatchField<scalar>&) {}
        virtual void operator+=(const fvPatchField<scalar>&) {}
        virtual void operator-=(const fvPatchField<scalar>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<scalar>&) {}
        virtual void operator-=(const Field<scalar>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const scalar&) {}
        virtual void operator+=(const scalar&) {}
        virtual void operator-=(const scalar&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
*/
private:
    //- initialise the object and compute the cell-index mapping
    void calcFaceMapping();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

