/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "explicitFOmembraneSoluteFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

explicitFOmembraneSoluteFvPatchScalarField::explicitFOmembraneSoluteFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchScalarField(p, iF),
    transProps_
    (
        IOobject
        (
            "transportProperties",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    UName_("U"),
    A_(readScalar(transProps_.lookup("A"))),
    B_(readScalar(transProps_.lookup("B"))),
    K_(readScalar(transProps_.lookup("K"))),
    D_AB_Min_(transProps_.lookup("D_AB_Min")),
    D_AB_Coeff_(transProps_.lookup("D_AB_Coeff")),
    D_AB_mACoeff_(transProps_.lookup("D_AB_mACoeff")),
    rho0_(transProps_.lookup("rho0")),
    rho_mACoeff_(transProps_.lookup("rho_mACoeff")),
    pi_mACoeff_(transProps_.lookup("pi_mACoeff")),
    fm_(p.size()),
    VIC_(p.size()),
    VBC_(p.size()),
    GIC_(p.size()),
    GBC_(p.size()),
    curTimeIndex_(-1)
{
    calcFaceMapping();
}

explicitFOmembraneSoluteFvPatchScalarField::explicitFOmembraneSoluteFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<scalar>(p, iF, dict),
    transProps_
    (
        IOobject
        (
            "transportProperties",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    A_(readScalar(transProps_.lookup("A"))),
    B_(readScalar(transProps_.lookup("B"))),
    K_(readScalar(transProps_.lookup("K"))),
    D_AB_Min_(transProps_.lookup("D_AB_Min")),
    D_AB_Coeff_(transProps_.lookup("D_AB_Coeff")),
    D_AB_mACoeff_(transProps_.lookup("D_AB_mACoeff")),
    rho0_(transProps_.lookup("rho0")),
    rho_mACoeff_(transProps_.lookup("rho_mACoeff")),
    pi_mACoeff_(transProps_.lookup("pi_mACoeff")),
    fm_(p.size()),
    VIC_(p.size()),
    VBC_(p.size()),
    GIC_(p.size()),
    GBC_(p.size()),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        /* Already Transfered, no need to do it again
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        */
    }
    else
    {
        // initialise the field to 0.0 if no information is given
        fvPatchField<scalar>::operator=(pTraits<scalar>::zero);
    }

    calcFaceMapping();
}

explicitFOmembraneSoluteFvPatchScalarField::explicitFOmembraneSoluteFvPatchScalarField
(
    const explicitFOmembraneSoluteFvPatchScalarField& empsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<scalar>(empsf, p, iF, mapper),
    transProps_
    (
        IOobject
        (
            "transportProperties",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    UName_(empsf.UName_),
    A_(empsf.A_),    
    B_(empsf.B_),
    K_(empsf.K_),
    D_AB_Min_(empsf.D_AB_Min_),
    D_AB_Coeff_(empsf.D_AB_Coeff_),
    D_AB_mACoeff_(empsf.D_AB_mACoeff_),
    rho0_(empsf.rho0_),
    rho_mACoeff_(empsf.rho_mACoeff_),
    pi_mACoeff_(empsf.pi_mACoeff_),
    fm_(p.size()),
    VIC_(p.size()),
    VBC_(p.size()),
    GIC_(p.size()),
    GBC_(p.size()),
    curTimeIndex_(-1)
{
}



explicitFOmembraneSoluteFvPatchScalarField::explicitFOmembraneSoluteFvPatchScalarField
(
    const explicitFOmembraneSoluteFvPatchScalarField& empsf
)
:
    fvPatchField<scalar>(empsf),
    transProps_(empsf.transProps_),
    UName_(empsf.UName_),
    A_(empsf.A_),    
    B_(empsf.B_),
    K_(empsf.K_),
    D_AB_Min_(empsf.D_AB_Min_),
    D_AB_Coeff_(empsf.D_AB_Coeff_),
    D_AB_mACoeff_(empsf.D_AB_mACoeff_),
    rho0_(empsf.rho0_),
    rho_mACoeff_(empsf.rho_mACoeff_),
    pi_mACoeff_(empsf.pi_mACoeff_),
    fm_(empsf.fm_),
    VIC_(empsf.size()),
    VBC_(empsf.size()),
    GIC_(empsf.size()),
    GBC_(empsf.size()),
    curTimeIndex_(-1)
{}

explicitFOmembraneSoluteFvPatchScalarField::explicitFOmembraneSoluteFvPatchScalarField
(
    const explicitFOmembraneSoluteFvPatchScalarField& empsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(empsf, iF),
    transProps_(empsf.transProps_),
    UName_(empsf.UName_),
    A_(empsf.A_),    
    B_(empsf.B_),
    K_(empsf.K_),
    D_AB_Min_(empsf.D_AB_Min_),
    D_AB_Coeff_(empsf.D_AB_Coeff_),
    D_AB_mACoeff_(empsf.D_AB_mACoeff_),
    rho0_(empsf.rho0_),
    rho_mACoeff_(empsf.rho_mACoeff_),
    pi_mACoeff_(empsf.pi_mACoeff_),
    fm_(empsf.fm_),
    VIC_(empsf.size()),
    VBC_(empsf.size()),
    GIC_(empsf.size()),
    GBC_(empsf.size()),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitFOmembraneSoluteFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
}


void explicitFOmembraneSoluteFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& psf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(psf, addr);
}


void explicitFOmembraneSoluteFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Execute the change to the openFraction only once per time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        // geometric help variables
        const vectorField& Sf = patch().Sf();
        
        // Delta coeffs
        // OpenFoam 2.2 and below: const scalarField& deltas = 1.0/patch().deltaCoeffs();
        tmp<scalarField> temp0 = 1.0/patch().deltaCoeffs();
        const scalarField& deltas = temp0();
        
        // Normal vectors
        tmp<vectorField> tvfnf = patch().nf();
        const vectorField& vfnf = tvfnf();

        // Velocity field along surface normal
        const fvPatchField<vector>& upvf = patch().lookupPatchField<volVectorField, vector>(UName_);
        
        // Velocity field
        // OpenFoam 2.2 and below: scalarField magU = max( mag( cmptMultiply(upvf,vfnf) ), VSMALL );
        tmp<scalarField> temp1 = max( mag( cmptMultiply(upvf,vfnf) ), VSMALL );
        const scalarField magU = temp1();

        /* Used in test of BCs. Only strictly required for debugging */
        const fvPatchScalarField& m_A = patch().lookupPatchField<volScalarField, scalar>("m_A");
        tmp<scalarField> temp2 = m_A.patchInternalField();
        // const scalarField& m_AInternal = temp2();

        // Total salt flux and area
        // scalar totalWeightFlux      = 0;
        // scalar totalArea            = 0;

        // Variables used
        scalar A                    = 0;
        scalar Js                   = 0;
        scalar B                    = 0;
        scalar rho                  = 0;

        // Debug variables
        // scalar newGrad              = 0;
        // scalar DrawMassInbalance    = 0;
        // scalar FeedmassInbalance    = 0;

        forAll(patch(), facei)
        {
            /* 
                Determine side based on surface normals & flow direction
                Afterwards implement following BC: D_AB*(dM_A/dn)+B*m_A = Js
            */
            rho = rho0_.value() * 
                  ( 1.0 + rho_mACoeff_.value()*operator[](facei) );
            A   = rho*max(D_AB_Coeff_*
                  ( 1.0 - D_AB_mACoeff_ * operator[](facei)) , D_AB_Min_ ).value();

            if ( (upvf[facei]&Sf[facei]) <= 0.0 )
            {
                // Variable B and salt flux:
                B  = rho*magU[facei];
                Js = ( B_ * magU[facei] ) / ( pi_mACoeff_.value() / 1000.0 * A_ );

                // Set coefficients            
                VIC_[facei] = 1.0 / (1.0 + B*deltas[facei]/A);
                VBC_[facei] = -Js / (A/deltas[facei] + B);
                GIC_[facei] = -1.0 / (A/B + deltas[facei]);
                GBC_[facei] = Js/( A + B*deltas[facei] );  

                /* Debugging:
                   Following tests whether the equation 
                   A*(dM_A/dn)+B*m_A = Js is satisfied
                */
                // newGrad = GIC_[facei]*m_AInternal[facei]+GBC_[facei];
                // DrawMassInbalance = A*newGrad+operator[](facei)*B - Js;
            }
            else
            {
                // Variable B and salt flux:
                B  = -rho*magU[facei];
                Js = ( B_ * magU[facei] ) / ( pi_mACoeff_.value() / 1000.0 * A_ );
           
                // Set coefficients            
                VIC_[facei] = 1.0 / (1.0 + B*deltas[facei]/A);
                VBC_[facei] = Js / (A/deltas[facei] + B);
                GIC_[facei] = -1.0 / (A/B + deltas[facei]);
                GBC_[facei] = Js/( A + B*deltas[facei] );  

                /* Debugging:
                   Following tests whether the equation 
                   A*(dM_A/dn)+B*m_A = Js is satisfied
                */
                // newGrad = GIC_[facei]*m_AInternal[facei]+GBC_[facei];
                // FeedmassInbalance = A*newGrad + operator[](facei)*B - Js; 

                // Add to salt flux through membrane
                // totalWeightFlux += patch().magSf()[facei] * Js;    
            }
            
        }

        // totalArea = sum(patch().magSf()) / 2.0;

        // Info << "Salt Flux: "
        //      << (totalWeightFlux*1e3*3600/totalArea) << " g/(h*m2) " 
        //      << ", Draw/Feed Flux Balance Residual: " << DrawMassInbalance << " / " << FeedmassInbalance
        //      << endl;

        
        // Set time index
        curTimeIndex_ = this->db().time().timeIndex();
    }
}


void explicitFOmembraneSoluteFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    operator==(this->patchInternalField() * VIC_ + VBC_);
}


tmp<Field<scalar> > explicitFOmembraneSoluteFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new Field<scalar>(VIC_));
}


tmp<Field<scalar> > explicitFOmembraneSoluteFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new Field<scalar>(VBC_));
}


tmp<Field<scalar> > explicitFOmembraneSoluteFvPatchScalarField::gradientInternalCoeffs() const
{
    return tmp<scalarField>(new Field<scalar>(GIC_));
}


tmp<Field<scalar> > explicitFOmembraneSoluteFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return tmp<scalarField>(new Field<scalar>(GBC_));
}

void explicitFOmembraneSoluteFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "m_A", "m_A", UName_);

    this->writeEntry("value", os);
}


void explicitFOmembraneSoluteFvPatchScalarField::calcFaceMapping()
{
    // set up the face-index mapping based on cell centres
    const vectorField& cfvf = patch().Cf();
    forAll(cfvf, facei)
    {
        for(label i=0; i<cfvf.size(); i++)
        {
            if (facei!=i)
            {
                if (mag(cfvf[facei]-cfvf[i])<1e-9)
                {
                    fm_[facei]=i;
                    if(debug)
                    {
                        Info << "patch face " << facei << " -> " << i << endl;
                    }
                    break;
                }
            }
        }
    }
}

makePatchTypeField
(
    fvPatchScalarField,
    explicitFOmembraneSoluteFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

