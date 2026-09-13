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

#include "explicitFOmembraneVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "membraneFaceMapping.H"
#include "fluxModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::explicitFOmembraneVelocityFvPatchVectorField::explicitFOmembraneVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    transProps_(this->db().lookupObject<IOdictionary>("transportProperties")),
    pName_("p"),
    m_AName_("m_A"),
    forwardDirection_(pTraits<vector>::zero),
    fluxEqName_("simple"),
    A_(readScalar(transProps_.lookup("A"))),
    B_(readScalar(transProps_.lookup("B"))),
    K_(readScalar(transProps_.lookup("K"))),
    slipName_("noSlip"),
    alpha_(1.0),
    aRelax_(1.0),
    pi_mACoeff_(transProps_.lookup("pi_mACoeff")),
    rho0_(transProps_.lookup("rho0")),
    rho_mACoeff_(transProps_.lookup("rho_mACoeff")),
    fm_(p.size()),
    fs_(p.size()/2)
{
    initialise();
}


Foam::explicitFOmembraneVelocityFvPatchVectorField::explicitFOmembraneVelocityFvPatchVectorField
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
    forwardDirection_(dict.lookupOrDefault<vector>("forwardDirection", pTraits<vector>::zero)),
    fluxEqName_(dict.lookupOrDefault<word>("eq", "simple")),
    A_(readScalar(transProps_.lookup("A"))),
    B_(readScalar(transProps_.lookup("B"))),
    K_(readScalar(transProps_.lookup("K"))),
    slipName_(dict.lookupOrDefault<word>("slip", "noSlip")),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 1.0)),
    aRelax_(dict.lookupOrDefault<scalar>("aRelax", 1.0)),
    pi_mACoeff_(transProps_.lookup("pi_mACoeff")),
    rho0_(transProps_.lookup("rho0")),
    rho_mACoeff_(transProps_.lookup("rho_mACoeff")),
    fm_(p.size()),
    fs_(p.size()/2)
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
        fvPatchField<vector>::operator=(pTraits<vector>::zero);
    }

    initialise();
}


Foam::explicitFOmembraneVelocityFvPatchVectorField::explicitFOmembraneVelocityFvPatchVectorField
(
    const explicitFOmembraneVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    transProps_(ptf.transProps_),
    pName_(ptf.pName_),
    m_AName_(ptf.m_AName_),
    forwardDirection_(ptf.forwardDirection_),
    fluxEqName_(ptf.fluxEqName_),
    A_(ptf.A_),
    B_(ptf.B_),
    K_(ptf.K_),
    slipName_(ptf.slipName_),
    alpha_(ptf.alpha_),
    aRelax_(ptf.aRelax_),
    pi_mACoeff_(ptf.pi_mACoeff_),
    rho0_(ptf.rho0_),
    rho_mACoeff_(ptf.rho_mACoeff_),
    fm_(p.size()),
    fs_(p.size()/2)
{
    initialise();
}


Foam::explicitFOmembraneVelocityFvPatchVectorField::explicitFOmembraneVelocityFvPatchVectorField
(
    const explicitFOmembraneVelocityFvPatchVectorField& efomvpvf
)
:
    fixedValueFvPatchVectorField(efomvpvf),
    transProps_(efomvpvf.transProps_),
    pName_(efomvpvf.pName_),
    m_AName_(efomvpvf.m_AName_),
    forwardDirection_(efomvpvf.forwardDirection_),
    fluxEqName_(efomvpvf.fluxEqName_),
    A_(efomvpvf.A_),
    B_(efomvpvf.B_),
    K_(efomvpvf.K_),
    slipName_(efomvpvf.slipName_),
    alpha_(efomvpvf.alpha_),
    aRelax_(efomvpvf.aRelax_),
    pi_mACoeff_(efomvpvf.pi_mACoeff_),
    rho0_(efomvpvf.rho0_),
    rho_mACoeff_(efomvpvf.rho_mACoeff_),
    fm_(efomvpvf.fm_),
    fs_(efomvpvf.fs_)
{}


Foam::explicitFOmembraneVelocityFvPatchVectorField::explicitFOmembraneVelocityFvPatchVectorField
(
    const explicitFOmembraneVelocityFvPatchVectorField& efomvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(efomvpvf, iF),
    transProps_(efomvpvf.transProps_),
    pName_(efomvpvf.pName_),
    m_AName_(efomvpvf.m_AName_),
    forwardDirection_(efomvpvf.forwardDirection_),
    fluxEqName_(efomvpvf.fluxEqName_),
    A_(efomvpvf.A_),
    B_(efomvpvf.B_),
    K_(efomvpvf.K_),
    slipName_(efomvpvf.slipName_),
    alpha_(efomvpvf.alpha_),
    aRelax_(efomvpvf.aRelax_),
    pi_mACoeff_(efomvpvf.pi_mACoeff_),
    rho0_(efomvpvf.rho0_),
    rho_mACoeff_(efomvpvf.rho_mACoeff_),
    fm_(efomvpvf.fm_),
    fs_(efomvpvf.fs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::explicitFOmembraneVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    initialise();
}


void Foam::explicitFOmembraneVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(pvf, addr);
    initialise();
}


void Foam::explicitFOmembraneVelocityFvPatchVectorField::updateCoeffs()
{
    if(updated())
    {
        return;
    }

    {
        // get the surface normals
        tmp<vectorField> tvfnf = patch().nf();
        const vectorField& vfnf = tvfnf();

        if(m_AName_=="none")
        {
            // get the pressure field
            const fvPatchField<scalar>& ppsf = patch().lookupPatchField<volScalarField, scalar>(pName_);

            // p is in Pa; A is a volumetric permeability in m/(Pa s).
            forAll(ppsf, facei)
            {
                // set the velocity
                operator[](facei) = vfnf[facei] * (A_*(ppsf[facei]-ppsf[fm_[facei]]));
            }
            
        }
        else
        {
            scalar feedMem       = 0;                    // Variable for feed membrane m_A
            scalar drawMem       = 0;                    // Variable for draw membrane m_A
            scalar i             = 0;                    // Total iterations counter
            scalar totalMassFlux = 0;                    // Total flux through membrane
            
            vector slipUinternal(0,0,0);                 // Slip BC component
            vector slipUboundary(0,0,0);                 // Slip BC component
            vector maxSlip(0,0,0);                       // Maximum Slip
            vector slipPrev(0,0,0);                      // Prev slip
            scalar dUdy = 0;                             // Strain rate

            // Get the mass fraction field
            const fvPatchScalarField& m_A = patch().lookupPatchField<volScalarField, scalar>(m_AName_);
            const scalarField& magSf = patch().magSf();
            
            // Get the current internal velocity field
            tmp<vectorField> temp            = this->patchInternalField();
            const vectorField& internalU     = temp();
            
            // Get cell-centre distances
            // OpenFoam 2.2 and below: const scalarField deltas = 1.0/patch().deltaCoeffs();
            tmp<scalarField> temp0 = 1.0/patch().deltaCoeffs();
            const scalarField& deltas = temp0();

            const bool resetFluxCache = cachedFeed_.size() != fs_.size();
            if (resetFluxCache)
            {
                cachedFeed_.setSize(fs_.size(), 0);
                cachedDraw_.setSize(fs_.size(), 0);
                cachedFlux_.setSize(fs_.size(), 0);
            }

            forAll(fs_, facei)
            {
                label fsi = fs_[facei];
                label dsi = fm_[fsi];
                feedMem = m_A[fsi];
                drawMem = m_A[dsi];

                if
                (
                    resetFluxCache
                 || cachedFeed_[facei] != feedMem
                 || cachedDraw_[facei] != drawMem
                )
                {
                    cachedFlux_[facei] = solveFlux(feedMem, drawMem, i);
                    cachedFeed_[facei] = feedMem;
                    cachedDraw_[facei] = drawMem;
                }
                const scalar flux = cachedFlux_[facei];
                
                // Calculate the velocity for the assymetric membrane
                vector v = vfnf[fsi] * flux;
                
                // Set the feed-side velocity
                operator[](fsi) = v;

                // Correct the velocity due to density change
                v *= (1.0 + rho_mACoeff_.value() * feedMem) / (1.0 + rho_mACoeff_.value() * drawMem);

                // Total flux and area
                totalMassFlux += flux * rho0_.value()*(1.0 + rho_mACoeff_.value() * feedMem) * magSf[dsi];

				// Slip boundary condition
                slipUboundary = vector::zero;
				if( slipName() == "slip" ){
                    slipUinternal = internalU[dsi] - (internalU[dsi] & vfnf[dsi]) * vfnf[dsi];
                    slipPrev = operator[](dsi) - (operator[](dsi) & vfnf[dsi])*vfnf[dsi];
					if( mag(slipUinternal) > SMALL ){
                        dUdy = (mag(slipUinternal) - mag(slipPrev) ) / deltas[dsi];
                        slipUboundary =  alpha()*dUdy * (slipUinternal/mag(slipUinternal)); 
                        slipUboundary = slipPrev + ( slipUboundary - slipPrev )*aRelax_;
						if( mag(slipUboundary) > mag(maxSlip) ){
							maxSlip = slipUboundary;
						}
					}
				}
                         
                // Set the draw-side velocity
                operator[](dsi) = v+slipUboundary;
            }
         
            Info << patch().name() << ": " << "Bracketed flux solve - Total iterations = " << i
                 << "\n    Water flux, " << fluxEqName_ << ": " << gSum(scalarField(1, totalMassFlux))/max(gSum(magSf)/2, VSMALL) * 3600 << " kg/(h*m2)"
                 << "\n    Draw/Feed m_A estimate: " << drawMem << " / " << feedMem 
                 << "\n    Max Slip Velocity: " << maxSlip << " with slip Coeff: " << alpha() << " and under-relax factor: " << aRelax_
                 << "\n    A: " << A_ << " / B: " << B_ << " / K: " << K_
                 << endl;
        }
    }
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::explicitFOmembraneVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("p", pName_);
    os.writeEntry("m_A", m_AName_);
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("B") << B_ << token::END_STATEMENT << nl;
    os.writeKeyword("K") << K_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    os.writeKeyword("eq") << fluxEqName_ << token::END_STATEMENT << nl;
    os.writeKeyword("aRelax") << aRelax_ << token::END_STATEMENT << nl;
    os.writeKeyword("forwardDirection") << forwardDirection_ << token::END_STATEMENT << nl;
    os.writeKeyword("slip") << slipName_ << token::END_STATEMENT << nl;
    os.writeKeyword("pi_mACoeff") << pi_mACoeff_.value() << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


void Foam::explicitFOmembraneVelocityFvPatchVectorField::initialise()
{
    calcFaceMapping();
    cachedFeed_.clear();
    cachedDraw_.clear();
    cachedFlux_.clear();

    if (fluxEqName_ != "simple" && fluxEqName_ != "advanced")
        FatalErrorInFunction << "eq must be simple or advanced" << exit(FatalError);
    if (slipName_ != "noSlip" && slipName_ != "slip")
        FatalErrorInFunction << "slip must be noSlip or slip" << exit(FatalError);
    if (mag(forwardDirection_) < VSMALL)
        FatalErrorInFunction << "Specify forwardDirection from feed into draw. "
            << "Inferring orientation from a field during construction is ambiguous."
            << exit(FatalError);
    fs_.setSize(patch().size()/2);
    const auto normals = patch().nf();
    label count = 0;
    forAll(fm_,i)
    {
        if (i > fm_[i]) continue;
        const scalar alignment = normals()[i] & forwardDirection_;
        if (mag(alignment) < SMALL)
            FatalErrorInFunction << "forwardDirection is tangent to membrane" << exit(FatalError);
        fs_[count++] = alignment > 0 ? i : fm_[i];
    }
    if (count != fs_.size())
        FatalErrorInFunction << "Incomplete membrane pairing" << exit(FatalError);
}

void Foam::explicitFOmembraneVelocityFvPatchVectorField::calcFaceMapping()
{
    membraneFaceMapping(patch(), fm_);
}


Foam::scalar Foam::explicitFOmembraneVelocityFvPatchVectorField::solveFlux
(
    const scalar& feedMem,
    const scalar& drawMem,
    scalar& i
)
{

    if (min(feedMem, drawMem) < -1e-12)
        FatalErrorInFunction
            << "Negative membrane mass fraction on patch " << patch().name()
            << ": feed=" << feedMem << ", draw=" << drawMem
            << ". Refine the mesh or reduce the time step." << exit(FatalError);
    int iterations = 0;
    try
    {
        const scalar result = membrane::flux(max(feedMem, scalar(0)), max(drawMem, scalar(0)), A_, B_, K_,
            pi_mACoeff_.value(), fluxEqName_ == "advanced", &iterations);
        i += iterations;
        return result;
    }
    catch (const std::exception& error)
    {
        FatalErrorInFunction << error.what()
            << " on patch " << patch().name()
            << ": feed=" << feedMem << ", draw=" << drawMem
            << exit(FatalError);
    }
    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        explicitFOmembraneVelocityFvPatchVectorField
    );
}

// ************************************************************************* //

