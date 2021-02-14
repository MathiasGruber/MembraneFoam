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

\*---------------------------------------------------------------------------*/

#include "explicitROmembraneVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::explicitROmembraneVelocityFvPatchVectorField::explicitROmembraneVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    transProps_(this->db().lookupObject<IOdictionary>("transportProperties")),
    pName_("p"),
    m_AName_("m_A"),
    K_(0.0),
    pi_mACoeff_("pi_mACoeff", dimless, transProps_),
    rho0_("rho0", dimDensity, transProps_),
    rho_mACoeff_("rho_mACoeff", dimless, transProps_),
    fm_(p.size())
{
    calcFaceMapping();
}


Foam::explicitROmembraneVelocityFvPatchVectorField::explicitROmembraneVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
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
    pName_(dict.lookupOrDefault<word>("p", "p")),
    m_AName_(dict.lookupOrDefault<word>("m_A", "m_A")),
    K_(readScalar(dict.lookup("K"))),
    pi_mACoeff_("pi_mACoeff", dimless, transProps_),
    rho0_("rho0", dimDensity, transProps_),
    rho_mACoeff_("rho_mACoeff", dimless, transProps_),
    fm_(p.size())
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        // initialise the field to (0, 0, 0) if no information is given
        fvPatchField<vector>::operator=(vector(0, 0, 0));
    }

    calcFaceMapping();
}


Foam::explicitROmembraneVelocityFvPatchVectorField::explicitROmembraneVelocityFvPatchVectorField
(
    const explicitROmembraneVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    transProps_(ptf.transProps_),
    pName_(ptf.pName_),
    m_AName_(ptf.m_AName_),
    K_(ptf.K_),
    pi_mACoeff_(ptf.pi_mACoeff_),
    rho0_(ptf.rho0_),
    rho_mACoeff_(ptf.rho_mACoeff_),
    fm_(p.size())
{
    calcFaceMapping();
}


Foam::explicitROmembraneVelocityFvPatchVectorField::explicitROmembraneVelocityFvPatchVectorField
(
    const explicitROmembraneVelocityFvPatchVectorField& emvpvf
)
:
    fixedValueFvPatchVectorField(emvpvf),
    transProps_(emvpvf.transProps_),
    pName_(emvpvf.pName_),
    m_AName_(emvpvf.m_AName_),
    K_(emvpvf.K_),
    pi_mACoeff_(emvpvf.pi_mACoeff_),
    rho0_(emvpvf.rho0_),
    rho_mACoeff_(emvpvf.rho_mACoeff_),
    fm_(emvpvf.fm_)
{}


Foam::explicitROmembraneVelocityFvPatchVectorField::explicitROmembraneVelocityFvPatchVectorField
(
    const explicitROmembraneVelocityFvPatchVectorField& emvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(emvpvf, iF),
    transProps_(emvpvf.transProps_),
    pName_(emvpvf.pName_),
    m_AName_(emvpvf.m_AName_),
    K_(emvpvf.K_),
    pi_mACoeff_(emvpvf.pi_mACoeff_),
    rho0_(emvpvf.rho0_),
    rho_mACoeff_(emvpvf.rho_mACoeff_),
    fm_(emvpvf.fm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::explicitROmembraneVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::explicitROmembraneVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(pvf, addr);
}


void Foam::explicitROmembraneVelocityFvPatchVectorField::updateCoeffs()
{
    if(updated())
    {
        return;
    }

    {
        // get the pressure field
        const fvPatchField<scalar>& ppsf = patch().lookupPatchField<volScalarField, scalar>(pName_);
        // get the surface normals
        tmp<vectorField> tvfnf = patch().nf();
        const vectorField& vfnf = tvfnf();

        if(m_AName_=="none")
        {
            forAll(ppsf, facei)
            {
                // set the velocity
                operator[](facei) = vfnf[facei] * (K_*(ppsf[facei]-ppsf[fm_[facei]]));
            }
        }
        else
        {
            const fvPatchScalarField& m_A = patch().lookupPatchField<volScalarField, scalar>(m_AName_);
            forAll(ppsf, facei)
            {
                            
                // set the velocity assuming the velocity is calculated for the upstream side
                vector v =
                    vfnf[facei] * K_ * 
                    (
                        (ppsf[facei]-ppsf[fm_[facei]])
                      - pi_mACoeff_.value()*(m_A[facei]-m_A[fm_[facei]])
                    );
                // test to see this is on the downstream side
                if ((v&vfnf[facei])<0.0)
                {
                    scalar rho_w = rho0_.value() * (1.0 + rho_mACoeff_.value()*m_A[fm_[facei]]);
                    scalar rho_p = rho0_.value() * (1.0 + rho_mACoeff_.value()*m_A[facei]);
                    // correct the downstream velocity
                    v *= rho_w / rho_p;
                }
                operator[](facei) = v;
            }
        }
    }
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::explicitROmembraneVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    os.writeEntryIfDifferent<word>("m_A", "m_A", m_AName_);
    os.writeKeyword("K") << K_ << token::END_STATEMENT << nl;
//    os.writeKeyword("pi_mACoeff") << pi_mACoeff_.value() << token::END_STATEMENT << nl;
//    os.writeKeyword("rho0") << rho0_ << token::END_STATEMENT << nl;
//    os.writeKeyword("rho_mACoeff") << rho_mACoeff_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


void Foam::explicitROmembraneVelocityFvPatchVectorField::calcFaceMapping()
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
                    if (debug)
                    {
                        Info << "patch face " << facei << " -> " << i << endl;
                    }
                    break;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        explicitROmembraneVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
