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
    
    // if this constructor is called by paraFoam, m_A is not part of the object registry and initialise is not executed
    List<word> objRegNames = db().names();
    forAll(objRegNames, iCounter)
    {
        if( objRegNames[iCounter] == m_AName_ )
        {
            initialise();
        }
    }
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
	// if this constructor is called by decomposePar, m_A is not part of the object registry and initialise is not executed
    List<word> objRegNames = db().names();
    forAll(objRegNames, iCounter)
    {
        if( objRegNames[iCounter] == m_AName_ )
        {
            initialise();
        }
    }
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
}


void Foam::explicitFOmembraneVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(pvf, addr);
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

            scalar rho = rho0_.value();
            forAll(ppsf, facei)
            {
                // set the velocity
                operator[](facei) = vfnf[facei] * (A_*(ppsf[facei]-ppsf[fm_[facei]])) / rho;
            }
            
        }
        else
        {
            scalar maxBound      = 1e-3;                 // Max flux, advancedPesDiff can't cope with higher fluxes!
            scalar minBound      = 1e-10;                // Min flux 1e-10
            scalar xacc          = 1e-10;                // Accuracy for ridder' method
            scalar feedMemSol       = 0;                 // Variable for feed membrane m_A
            scalar drawMemSol       = 0;                 // Variable for draw membrane m_A
            scalar feedMemPres       = 0;                // Variable for feed membrane pressure
            scalar drawMemPres       = 0;                // Variable for draw membrane pressure
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

            // Get the pressure field
	    	const fvPatchField<scalar>& p = patch().lookupPatchField<volScalarField, scalar>(pName_);
            
            // Get the current internal velocity field
            const fvPatchVectorField& Ufield = patch().lookupPatchField<volVectorField, vector>("U");
            tmp<vectorField> temp            = Ufield.patchInternalField();
            const vectorField& internalU     = temp();
            //tmp<vectorField> temp2           = Ufield.snGrad();
            //const vectorField& gradU         = temp2();
            
            // Get cell-centre distances
            // OpenFoam 2.2 and below: const scalarField deltas = 1.0/patch().deltaCoeffs();
            tmp<scalarField> temp0 = 1.0/patch().deltaCoeffs();
            const scalarField& deltas = temp0();

            forAll(fs_, facei)
            {
                label fsi = fs_[facei];
                label dsi = fm_[fsi];
                feedMemSol = m_A[fsi];
                drawMemSol = m_A[dsi];
                feedMemPres = p[fsi];
                drawMemPres = p[dsi];

                // User ridder's method to solve for the flux. Save flux in rtn.
                scalar flux = ridderSolve(feedMemSol, drawMemSol, feedMemPres, drawMemPres, minBound, maxBound, xacc, i); 
                
                // Calculate the velocity for the assymetric membrane
                vector v = vfnf[fsi] * flux;
                
                // Set the feed-side velocity
                operator[](fsi) = v;

                // Correct the velocity due to density change
                v *= (1.0 + rho_mACoeff_.value() * feedMemSol) / (1.0 + rho_mACoeff_.value() * drawMemSol);

                // Total flux and area
                totalMassFlux += flux * rho0_.value()*(1.0 + rho_mACoeff_.value() * feedMemSol) * magSf[dsi];

				// Slip boundary condition
				if( slipName() == "slip" ){
                    slipUinternal = internalU[dsi] - (internalU[dsi] & vfnf[dsi]) * vfnf[dsi];
                    slipPrev = operator[](dsi) - cmptMultiply (cmptMultiply ( vfnf [ dsi ] , vfnf [ dsi ] ) , operator[](dsi) );
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
         
            Info << patch().name() << ": " << "Ridders' Method - Total iterations = " << i
                 << "\n    Water flux, " << fluxEqName_ << ": " << totalMassFlux/(sum(magSf)/2) * 3600 << " kg/(h*m2)" 
                 << "\n    Draw/Feed m_A estimate: " << drawMemSol << " / " << feedMemSol 
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
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "m_A", "m_A", m_AName_);
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

    // fill out the fs_ list so that it contains the indices of all the "feed-side" faces
    if (mag(forwardDirection_) < VSMALL)
    {
        Info << patch().name() << ": forward direction specified by mass fraction\n" << endl;

        // the forward direction of the membrane has not been defined by the user
        // so the mass fraction will be used to determine the feed side
        const fvPatchScalarField& m_Apsf = patch().lookupPatchField<volScalarField, scalar>(m_AName_);
        tmp<scalarField> tm_Asf = m_Apsf.patchInternalField();
        const scalarField& m_Asf = tm_Asf();
        forAll(fs_, facei)
        {
            fs_[facei] = (m_Asf[facei] > m_Asf[fm_[facei]]) ? fm_[facei] : facei;
        }
    }
    else
    {
        Info << patch().name() << ": forward direction specified by user\n" << endl;

        tmp<vectorField> tvfnf = patch().nf();
        const vectorField& vfnf = tvfnf();
        forAll(fs_, facei)
        {
            fs_[facei] = ((vfnf[facei] & forwardDirection_) < 0.0) ? fm_[facei] : facei;
        }
    }
}

void Foam::explicitFOmembraneVelocityFvPatchVectorField::calcFaceMapping()
{
    // set up the face-index mapping based on cell centres
    const vectorField& cfvf = patch().Cf();
    forAll(cfvf, facei)
    {
        forAll(cfvf, i)
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


Foam::scalar Foam::explicitFOmembraneVelocityFvPatchVectorField::fluxEquation( const scalar& Jvalue, 
                                                                         const scalar& feedm_A, 
                                                                         const scalar& drawm_A,
                                                                         const scalar& feedP,
                                                                         const scalar& drawP )
{

    if( fluxEqName_ == "simple" || B() < SMALL )
    {
        /*- Implicit flux equation,
        Valid when B is low compared to other terms, i.e. high rejection.
        See "Modelling Water Flux in Forward Osmosis: Implications for Improved Membrane Design
        American institude of Chemical Engineers, vol 53, No. 7, p. 1736-1744
        */
        return Jvalue - A()*( pi_mACoeff().value() * ( drawm_A*exp(-Jvalue* K() ) - feedm_A ) );
    }
    else if( fluxEqName_ == "noFlux" )
    {
        return 0;
    }
    else if( fluxEqName_ == "advanced" )
    {
        /*- Implicit flux equation
        Valid at any B-value.
        See "Coupled effects of internal concentration polarization and fouling on flux 
             behaviors of forward osmosis membranes during humic acid filtration"
        Journal of Membrane Science 354 (2010) 123-133
        */
        scalar numerator    = A()*pi_mACoeff().value()*drawm_A + B();
        scalar denominator  = A()*pi_mACoeff().value()*feedm_A + Jvalue + B();
        
        // To avoid floating point exceptions
        if( denominator > SMALL ){
            return Jvalue - ( 1/K() ) * log( numerator / denominator );
        }
        else{
            return 0;
        }
    }
    else if( fluxEqName_ == "advancedPresDiff" )
    {
        /*- Implicit flux equation
        Like "advanced" flux equation but assuming a non-zero hydraulic pressure difference across the membrane.
        */

        scalar expJK = exp( Jvalue*K() );
        scalar dP = drawP - feedP;
        scalar numerator = pi_mACoeff().value()*feedm_A * Jvalue * ( drawm_A/feedm_A - expJK );
        scalar denominator = ( Jvalue + B() ) * expJK - B();
             
        // To avoid floating point exceptions
        if( denominator > SMALL ){
            return Jvalue -  A()  * ( numerator/denominator - dP );
        }
        else{
            return 0;
        }
    }
    else
    {
        FatalErrorIn
        (
            "In the file: explicitFOmembraneVelocity.C"
        ) << "No flux model was selected. Select either of the following models in 0/U, eq = {simple, advanced, advancedPresDiff}. " << abort(FatalError);
    }
    return 0;
}


Foam::scalar Foam::explicitFOmembraneVelocityFvPatchVectorField::ridderSolve( const scalar& feedMemSol,
                                                                        const scalar& drawMemSol,
                                                                        const scalar& feedMemPres,
                                                                        const scalar& drawMemPres,
                                                                        const scalar& minBound,
                                                                        const scalar& maxBound,
                                                                        const scalar& xacc,
                                                                        scalar& i )
{

    // Function of boundaries
    scalar fl = fluxEquation( minBound , feedMemSol , drawMemSol , feedMemPres , drawMemPres );
    scalar fh = fluxEquation( maxBound , feedMemSol , drawMemSol , feedMemPres , drawMemPres );

    if( (fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0) )
    {
        i++;

        // Save bounds in new variables
        scalar xl = minBound;
        scalar xh = maxBound;

        // An unlikely value, to simplify logic below
        scalar ans = -1.11e-30;

        // Variables used
        scalar xm, fm, s, xnew, fnew;

        // iteration counter
        for( int j=0 ; j<50 ; j++ )
        {
            xm = 0.5*(xl+xh);
            fm = fluxEquation( xm , feedMemSol , drawMemSol , feedMemPres , drawMemPres );

            // First of two function evaluations
            s = sqrt( fm*fm - fl*fh );
            if( s < SMALL )
            {
                return ans;
            }

            // Update the formula and check answer
            if (fl >= fh)
            {
                xnew = xm + (xm - xl)*fm/s;
            }
            else
            {
                xnew = xm - (xm - xl)*fm/s;
            }
            if ( mag( xnew - ans ) <= xacc)
            {
                return ans;
            }

            ans = xnew;

            fnew = fluxEquation( ans , feedMemSol , drawMemSol , feedMemPres , drawMemPres );

            if ( mag(fnew) < SMALL )
            {
                return ans;
            }

            // Bookkeeping to keep root bracketed on next iteration
            if (checkSign(fm,fnew))
            {
                xl = xm;
                xh = fm;
                xh = ans;
                fh = fnew;
            }
            else if (checkSign(fl,fnew) )
            {
                xh = ans;
                fh = fnew;
            } 
            else if (checkSign(fh,fnew) )
            {
                xl = ans;
                fl=fnew;
            } 
            else 
            {
                FatalErrorIn
                (
                    "In the file: explicitFOmembraneVelocity.C"
                )   << "Error in search logic" << abort(FatalError);
            }
            // Check the bounds
            if ( mag( xh - xl ) <= xacc)
            {
                return ans;
            }
        }
    }        
    else 
    {
        if( mag(fl) < SMALL )
        {
            return minBound;
        }
        if( mag(fh) < SMALL )
        {
            return maxBound;
        }
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

