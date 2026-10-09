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

#include "explicitROmembraneSoluteFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "membraneFaceMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

explicitROmembraneSoluteFvPatchScalarField::explicitROmembraneSoluteFvPatchScalarField
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
    R_(0),
    D_AB_Min_(transProps_.lookup("D_AB_Min")),
    D_AB_Coeff_(transProps_.lookup("D_AB_Coeff")),
    D_AB_mACoeff_(transProps_.lookup("D_AB_mACoeff")),
    rho0_(transProps_.lookup("rho0")),
    rho_mACoeff_(transProps_.lookup("rho_mACoeff")),
    fm_(p.size()),
    VIC_(p.size(), 1),
    VBC_(p.size(), 0),
    GIC_(p.size(), 0),
    GBC_(p.size(), 0)
{
    calcFaceMapping();
}

explicitROmembraneSoluteFvPatchScalarField::explicitROmembraneSoluteFvPatchScalarField
(
    const explicitROmembraneSoluteFvPatchScalarField& empsf,
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
    R_(empsf.R_),
    D_AB_Min_(empsf.D_AB_Min_),
    D_AB_Coeff_(empsf.D_AB_Coeff_),
    D_AB_mACoeff_(empsf.D_AB_mACoeff_),
    rho0_(empsf.rho0_),
    rho_mACoeff_(empsf.rho_mACoeff_),
    fm_(p.size()),
    VIC_(p.size(), 1),
    VBC_(p.size(), 0),
    GIC_(p.size(), 0),
    GBC_(p.size(), 0)
{}

explicitROmembraneSoluteFvPatchScalarField::explicitROmembraneSoluteFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<scalar>(p, iF),
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
    R_(readScalar(dict.lookup("R"))),
    D_AB_Min_(transProps_.lookup("D_AB_Min")),
    D_AB_Coeff_(transProps_.lookup("D_AB_Coeff")),
    D_AB_mACoeff_(transProps_.lookup("D_AB_mACoeff")),
    rho0_(transProps_.lookup("rho0")),
    rho_mACoeff_(transProps_.lookup("rho_mACoeff")),
    fm_(p.size()),
    VIC_(p.size(), 1),
    VBC_(p.size(), 0),
    GIC_(p.size(), 0),
    GBC_(p.size(), 0)
{
    if (!(R_ >= 0 && R_ <= 1))
        FatalIOErrorInFunction(dict) << "R must be between zero and one" << exit(FatalIOError);
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        // initialise the field to 0.0 if no information is given
        fvPatchField<scalar>::operator=(pTraits<scalar>::zero);
    }

    calcFaceMapping();
}

explicitROmembraneSoluteFvPatchScalarField::explicitROmembraneSoluteFvPatchScalarField
(
    const explicitROmembraneSoluteFvPatchScalarField& empsf
)
:
    fvPatchField<scalar>(empsf),
    transProps_(empsf.transProps_),
    UName_(empsf.UName_),
    R_(empsf.R_),
    D_AB_Min_(empsf.D_AB_Min_),
    D_AB_Coeff_(empsf.D_AB_Coeff_),
    D_AB_mACoeff_(empsf.D_AB_mACoeff_),
    rho0_(empsf.rho0_),
    rho_mACoeff_(empsf.rho_mACoeff_),
    fm_(empsf.fm_),
    VIC_(empsf.size(), 1),
    VBC_(empsf.size(), 0),
    GIC_(empsf.size(), 0),
    GBC_(empsf.size(), 0)
{}

explicitROmembraneSoluteFvPatchScalarField::explicitROmembraneSoluteFvPatchScalarField
(
    const explicitROmembraneSoluteFvPatchScalarField& empsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(empsf, iF),
    transProps_(empsf.transProps_),
    UName_(empsf.UName_),
    R_(empsf.R_),
    D_AB_Min_(empsf.D_AB_Min_),
    D_AB_Coeff_(empsf.D_AB_Coeff_),
    D_AB_mACoeff_(empsf.D_AB_mACoeff_),
    rho0_(empsf.rho0_),
    rho_mACoeff_(empsf.rho_mACoeff_),
    fm_(empsf.fm_),
    VIC_(empsf.size(), 1),
    VBC_(empsf.size(), 0),
    GIC_(empsf.size(), 0),
    GBC_(empsf.size(), 0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitROmembraneSoluteFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<scalar>::autoMap(m);
    VIC_.setSize(size(), 1);
    VBC_.setSize(size(), 0);
    GIC_.setSize(size(), 0);
    GBC_.setSize(size(), 0);
    calcFaceMapping();
}


void explicitROmembraneSoluteFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& psf,
    const labelList& addr
)
{
    fvPatchField<scalar>::rmap(psf, addr);
    VIC_.setSize(size(), 1);
    VBC_.setSize(size(), 0);
    GIC_.setSize(size(), 0);
    GBC_.setSize(size(), 0);
    calcFaceMapping();
}


void explicitROmembraneSoluteFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    // Geometric help variables
    const vectorField& Sf = patch().Sf();
    
    // OpenFoam 2.2 and below: const scalarField& deltas = 1.0/patch().deltaCoeffs();
    tmp<scalarField> temp0 = 1.0/patch().deltaCoeffs();
    const scalarField& deltas = temp0();
    
    // Get the velocity field
    const fvPatchField<vector>& upvf = patch().lookupPatchField<volVectorField, vector>(UName_);
    
    // OpenFoam 2.2 and below: scalarField magU = max( mag(upvf & Sf) / mag(Sf) , VSMALL);
    tmp<scalarField> temp1 = max( mag(upvf & Sf) / mag(Sf) , VSMALL);
    const scalarField magU = temp1();
    
    forAll(patch(), facei)
    {
        // test the velocity
        if ((upvf[facei]&Sf[facei])>=0.0)
        {
            // the diffusivity is calculated from the current values for m_A, 
            // in order to linearise the boundary condition
            scalar D_AB = max(D_AB_Coeff_*(1.0-D_AB_mACoeff_*operator[](facei)),D_AB_Min_).value();

            if (magU[facei]*R_*deltas[facei] >= D_AB)
                FatalErrorInFunction << "Unresolved RO boundary layer on face " << facei
                    << ": refine the wall-normal mesh (U_n*R*delta/D must be < 1)"
                    << exit(FatalError);

            // upstream side (flow is directed out of the adjacent cell)
            VIC_[facei] = 1.0 / (1.0 - magU[facei]*R_*deltas[facei]/D_AB);
            VBC_[facei] = 0.0;
            GIC_[facei] = magU[facei]*R_ / (D_AB - magU[facei]*R_*deltas[facei]);
            GBC_[facei] = 0.0;
        }
        else
        {
            scalar D_AB = max(D_AB_Coeff_*(1.0-D_AB_mACoeff_*operator[](facei)),D_AB_Min_).value();

            // downstream side (flow is directed into the adjacent cell)
            VIC_[facei] = 1.0 / (1.0 + magU[facei]*deltas[facei]/D_AB);
            VBC_[facei] = (1.0-R_)*operator[](fm_[facei])*magU[facei]*deltas[facei] / (D_AB + magU[facei]*deltas[facei]);
            GIC_[facei] = -magU[facei] / (D_AB + magU[facei]*deltas[facei]);
            GBC_[facei] = (1.0-R_)*operator[](fm_[facei])*magU[facei] / (D_AB + magU[facei]*deltas[facei]);
        }
    }
}


void explicitROmembraneSoluteFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    operator==(this->patchInternalField() * VIC_ + VBC_);
}


tmp<Field<scalar> > explicitROmembraneSoluteFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new Field<scalar>(VIC_));
}


tmp<Field<scalar> > explicitROmembraneSoluteFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new Field<scalar>(VBC_));
}


tmp<Field<scalar> > explicitROmembraneSoluteFvPatchScalarField::gradientInternalCoeffs() const
{
    return tmp<scalarField>(new Field<scalar>(GIC_));
}


tmp<Field<scalar> > explicitROmembraneSoluteFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return tmp<scalarField>(new Field<scalar>(GBC_));
}

void explicitROmembraneSoluteFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("U", UName_);
    os.writeKeyword("R") << R_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}


void explicitROmembraneSoluteFvPatchScalarField::calcFaceMapping()
{
    membraneFaceMapping(patch(), fm_);
}

makePatchTypeField
(
    fvPatchScalarField,
    explicitROmembraneSoluteFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
