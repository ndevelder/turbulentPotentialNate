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

#include "turbulentPotentialNate.H"
#include "addToRunTimeSelectionTable.H"
#include "backwardsCompatibilityWallFunctions.H"
#include "components.H"
#include "fvCFD.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulentPotentialNate, 0);
addToRunTimeSelectionTable(RASModel, turbulentPotentialNate, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


tmp<volScalarField> turbulentPotentialNate::Ts() const
{
	if(tslimiter_ == "true")
	{
        return max(k_/epsilon_, minTS());
	}
	
    return (k_/epsilon_);
}

tmp<volScalarField> turbulentPotentialNate::TsEh() const
{
	if(tslimiter_ == "true")
	{
        return max(1.0/epsHat_, minTS());
	}
	
    return (1.0/epsHat_);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentPotentialNate::turbulentPotentialNate
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),


    cEp1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEp1",
            coeffDict_,
            1.45
        )
    ),
    cEp2con_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEp2con",
            coeffDict_,
            1.83
        )
    ),
    cEp3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEp3",
            coeffDict_,
            0.15
        )
    ),
    cD1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cD1",
       	    coeffDict_,
            0.5
        )
    ),
    cD2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cD2",
       	    coeffDict_,
            0.88
        )
    ),
    cVv1_
    (
     	dimensioned<scalar>::lookupOrAddToDict
        (
            "cVv1",
            coeffDict_,
            2.0
        )
    ),
    cTv1_
    (
     	dimensioned<scalar>::lookupOrAddToDict
        (
            "cTv1",
            coeffDict_,
            0.0
        )
    ),
    cP1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP1",
            coeffDict_,
            2.0
        )
    ),
    cP2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP2",
            coeffDict_,
            0.6
        )
    ),
    cP3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP3",
            coeffDict_,
            0.12
        )
    ),
    cP4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cP4",
            coeffDict_,
            0.85714
        )
    ),
    cPphi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cPphi",
            coeffDict_,
            3.2
        )
    ),
    cMu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cMu",
            coeffDict_,
            0.21
        )
    ),

    cT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cT",
            coeffDict_,
            0.0033
        )
    ),

    cPr_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cPr",
            coeffDict_,
            1.0
        )
    ),

    cEhm_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEhm",
            coeffDict_,
            10.0
        )
    ),
	
    cEhR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cEhR",
            coeffDict_,
            1.0
        )
    ),

	gT1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gT1",
            coeffDict_,
            0.0
        )
    ),
	
	gT2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gT2",
            coeffDict_,
            0.0
        )
    ),
	gT3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gT3",
            coeffDict_,
            0.0
        )
    ),	
	cNF_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cNF",
            coeffDict_,
            10.0
        )
    ),
    
    	cPw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cPw",
            coeffDict_,
            18.0
        )
    ),
    
    sigmaKInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaKInit",
            coeffDict_,
            1.0
        )
    ),

    sigmaEpsInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEpsInit",
            coeffDict_,
            0.833
        )
    ),

    sigmaEpsVisc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEpsVisc",
            coeffDict_,
            1.0
        )
    ),

    sigmaPhiInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPhiInit",
            coeffDict_,
            0.33
        )
    ),

    sigmaPsiInit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaPsiInit",
            coeffDict_,
            1.0
        )
    ),

    psiNuFrac_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "psiNuFrac",
            coeffDict_,
            1.0
        )
    ),


   solveK_
   (
       coeffDict_.lookup("solveK")
   ),

   solveEps_
   (
       coeffDict_.lookup("solveEps")
   ),

   solvePsi_
   (
       coeffDict_.lookup("solvePsi")
   ),

   solvePhi_
   (
       coeffDict_.lookup("solvePhi")
   ),

   solveNut_
   (
       coeffDict_.lookup("solveNut")
   ),

   eqnSigmaK_
   (
       coeffDict_.lookup("eqnSigmaK")
   ),

   eqnSigmaEps_
   (
       coeffDict_.lookup("eqnSigmaEps")
   ),

   eqnSigmaPhi_
   (
       coeffDict_.lookup("eqnSigmaPhi")
   ),

   eqnSigmaPsi_
   (
       coeffDict_.lookup("eqnSigmaPsi")
   ),

   eqncEp2_
   (
       coeffDict_.lookup("eqncEp2")
   ),

   eqnEpsHat_
   (
       coeffDict_.lookup("eqnEpsHat")
   ),
   
   timeScaleEps_
   (
       coeffDict_.lookup("timeScaleEps")
   ),
   prodType_
   (
       coeffDict_.lookup("prodType")
   ),
   debugWrite_
   (
       coeffDict_.lookup("debugWrite")
   ),
   tslimiter_
   (
       coeffDict_.lookup("tslimiter")
   ),
   psiProd_
   (
       coeffDict_.lookup("psiProd")
   ),
    y_(mesh_),


    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    gradk_
    (
        IOobject
        (
            "gradk",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(k_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    nutNorm_
    (
        IOobject
        (
            "nutNorm",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (nut_/max(nut_))
    ),
    tpphi_
    (
        IOobject
        (
            "tpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    tpphisqrt_
    (
        IOobject
        (
            "tpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt(tpphi_)
    ),
    vorticity_
    (
        IOobject
        (
            "vorticity",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(U_)
    ),
	phis_
    (
        IOobject
        (
            "phis",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    vorticityTmp_
    (
        IOobject
        (
            "vorticityTmp",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::curl(U_)
    ),
    ivorticity_
    (
        IOobject
        (
            "ivorticity",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("iv", dimensionSet(0,0,1,0,0,0,0), vector(1.0,1.0,1.0))
    ),
    tppsi_
    (
        IOobject
        (
            "tppsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    uGrad_
    (
        IOobject
        (
            "uGrad",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(U_)
    ),
	epsHat_
    (
        IOobject
        (
            "epsHat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsHat", dimensionSet(0,0,-1,0,0,0,0), 1.0)
    ),
    eHrC_
    (
        IOobject
        (
            "eHrC",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        cEhR_*dimensionedScalar("eHrC", dimensionSet(0,1,0,0,0,0,0), 1.0)
    ),
    kol_
    (
        IOobject
        (
            "kol",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kol", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    ),
    kSafe_
    (
        IOobject
        (
            "kSafe",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (k_)
    ),
    kSqrt_
    (
        IOobject
        (
            "kSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    gradkSqrt_
    (
        IOobject
        (
            "gradkSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(kSqrt_)
    ),
    nutSafe_
    (
        IOobject
        (
            "nutSafe",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (nut_)
    ),
    epsilonSafe_
    (
        IOobject
        (
            "epsilonSafe",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (epsilon_)
    ),
    sigmaK_
    (
        IOobject
        (
            "sigmaK",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaKInit_
    ),
    sigmaEps_
    (
        IOobject
        (
            "sigmaEps",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaEpsInit_
    ),
    sigmaPhi_
    (
        IOobject
        (
            "sigmaPhi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaPhiInit_
    ),
    sigmaPsi_
    (
        IOobject
        (
            "sigmaPsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        sigmaPsiInit_
    ),
    cEp2_
    (
        IOobject
        (
            "cEp2",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (cEp2con_ - 0.16*exp(-0.1*sqr(k_)/(nu()*epsilon_)))
    ),
    tpProd_
    (
        IOobject
        (
            "tpProd",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (2*nut_*magSqr(symm(fvc::grad(U_))))
    ),
    cP1eqn_
    (
        IOobject
        (
            "cP1eqn",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        2.0*(0.5+0.5*((tpProd_)/epsilon_))
    ),
    dimRat_
    (
        IOobject
        (
            "dimRat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (psiReal() & psiReal())/(k_*phiReal())
    ),
    gradTpphi_
    (
        IOobject
        (
            "gradTpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(tpphi_)
    ),
    gradTppsi_
    (
        IOobject
        (
            "gradTppsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(tppsi_)
    ),
    tpProdSqr_
    (
        IOobject
        (
            "tpProdSqr",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(tppsi_ & vorticity_)
    ),
    tpProd3d_
    (
        IOobject
        (
            "tpProd3d",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (mag(psiReal() ^ vorticity_))
    ),
	phiPressureStrain
    (
        IOobject
        (
            "phiPressureStrain",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiPressureStrain", dimensionSet(0,2,-3,0,0,0,0), 1.0)
    ),
    phiPressureDiff
    (
        IOobject
        (
            "phiPressureDiff",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiPressureDiff", dimensionSet(0,2,-3,0,0,0,0), 1.0)
    ),
    phiDiss
    (
        IOobject
        (
            "phiDiss",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiDiss", dimensionSet(0,2,-3,0,0,0,0), 1.0)
    ),
    phiViscDiff
    (
        IOobject
        (
            "phiViscDiff",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiViscDiff", dimensionSet(0,2,-3,0,0,0,0), 1.0)
    ),
    phiTurbDiff
    (
        IOobject
        (
            "phiTurbDiff",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiTurbDiff", dimensionSet(0,2,-3,0,0,0,0), 1.0)
    )
{

    Info<< "Made it past constructors " << endl;

    // Calculate eddy viscosity
    if(solveNut_ == "true")
    {
		if(timeScaleEps_ == "epsilon" || timeScaleEps_ != "epsHat")
		{
            nut_ = cMu_*tpphi_*Ts();
        }
        
        if(timeScaleEps_ == "epsHat")
		{
            nut_ = cMu_*tpphi_*TsEh();
        }
        
    }

    kSafe_ = max(k_, dimensionedScalar("minK", k_.dimensions(), 1.0e-15));

    Info<< "solveK is: " <<solveK_ <<endl;
    Info<< "solveEps is: " <<solveEps_ <<endl;
    Info<< "solvePhi is: " <<solvePhi_ <<endl;
    Info<< "solvePsi is: " <<solvePsi_ <<endl;

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotentialNate::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (
			  (2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)
			)
        )
    );
}

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotentialNate::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

// Term that is directly added to the momentum equation
tmp<fvVectorMatrix> turbulentPotentialNate::divDevReff(volVectorField& U) const
{
    return
    (
       fvc::grad(tpphi_)
     + fvc::curl(tppsi_)
     + fvc::laplacian(nut_, U, "laplacian(nuEff,U)")
     - fvm::laplacian(nuEff(), U)
    );
}



bool turbulentPotentialNate::read()
{
    if (RASModel::read())
    {
        cEp1_.readIfPresent(coeffDict());
        cEp2con_.readIfPresent(coeffDict());
        cEp3_.readIfPresent(coeffDict());
        cP1_.readIfPresent(coeffDict());
        cP2_.readIfPresent(coeffDict());
        cP3_.readIfPresent(coeffDict());
        cPw_.readIfPresent(coeffDict());
        cPphi_.readIfPresent(coeffDict());
        cMu_.readIfPresent(coeffDict());
		cPphi_.readIfPresent(coeffDict());
		cEhm_.readIfPresent(coeffDict());
		cEhR_.readIfPresent(coeffDict());
		cPr_.readIfPresent(coeffDict());
		cD1_.readIfPresent(coeffDict());
		cD2_.readIfPresent(coeffDict());
        cVv1_.readIfPresent(coeffDict());
        cTv1_.readIfPresent(coeffDict());
		cT_.readIfPresent(coeffDict());
		sigmaKInit_.readIfPresent(coeffDict());
        sigmaEpsInit_.readIfPresent(coeffDict());
        sigmaEpsVisc_.readIfPresent(coeffDict());
        sigmaPhiInit_.readIfPresent(coeffDict());
		sigmaPsiInit_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void turbulentPotentialNate::correct()
{

    if (mesh_.changing())
    {
        bound(k_, dimensionedScalar("minK", k_.dimensions(), 1.0e-10));
        bound(epsilon_, dimensionedScalar("minEps", epsilon_.dimensions(), 1.0e-10));
		bound(tpphi_,dimensionedScalar("minTpphi", tpphi_.dimensions(), 1.0e-10));
    }
    
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }
    
    
    if (mesh_.changing())
    {
        y_.correct();
    }

    // Set the time scale using either epsilon or epsHat
    volScalarField T("TimeScale",Ts());
    volScalarField Teh("TimescaleEpsHat",TsEh());
	
    if(timeScaleEps_ == "epsilon" || timeScaleEps_ != "epsHat")
    {
        T = Ts();
        bound(T, dimensionedScalar("minT", T.dimensions(), 1.0e-10));
    }
        
    if(timeScaleEps_ == "epsHat")
    {
        T = TsEh();
        bound(T, dimensionedScalar("minT", T.dimensions(), 1.0e-10));
    }
        
	
    // Vorticity
    vorticity_ = fvc::curl(U_);
	uGrad_ = fvc::grad(U_);	
    //Info<< "Made it past vorticity" << endl;
	
	
	// Production initialize
	const volScalarField S2(2*magSqr(symm(fvc::grad(U_))));
    volScalarField G("RASModel::G", nut_*S2);
	
	
    // Production set
	if(prodType_ == "strain")
	{
		Info<< "Using nut*S^2 for production" << endl;
		const volScalarField S2(2*magSqr(symm(fvc::grad(U_))));
		volScalarField G("RASModel::G", nut_*S2);
		tpProd_ = tppsi_ & vorticity_;
	}
	else
	{
		Info<< "Using psi*vorticity for production" << endl;
		tpProd_ = tppsi_ & vorticity_;
		G = tpProd_;
	}
    //Info<< "Made it production" << endl;

	
    // Epsilon-hat (Epsilon - Epsilon Wall)
	if(eqnEpsHat_ == "mod")
	{
        epsHat_ = (epsilon_)/(k_ + cEhm_*nu()*mag(gradkSqrt_));
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), 1.0e-10));
	}
	else if(eqnEpsHat_ == "dif")
	{
        epsHat_ = (epsilon_ - 2.0*nu()*sqr(mag(gradkSqrt_)))/k_;
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), 1.0e-10));
	}
	else if(eqnEpsHat_ == "eps")
	{
        epsHat_ = epsilon_/k_;
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), 1.0e-10));
	}
	else if(eqnEpsHat_ == "rough")
	{
        epsHat_ = epsilon_/k_;
        bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), 1.0e-10));
	}
	else
	{
        Info<< "No EpsHat Model Chosen" <<endl;
	    epsHat_ = (epsilon_)/(k_ + cEhm_*nu()*mag(fvc::grad(kSqrt_)));
	    bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), 1.0e-10));
	}

	
	// Sigma equations 
    if(eqnSigmaK_ == "true")
    {
	    sigmaK_ = 0.67 + 0.33*(tpProd_/(epsHat_*k_));
	}

    if(eqnSigmaEps_ == "true")
    {
	    sigmaEps_ = 0.33 + 0.5*(tpProd_/(epsHat_*k_));
    }

    if(eqnSigmaPhi_ == "true")
    {
	    sigmaPhi_ = 0.21 + 0.12*(tpProd_/(epsHat_*k_));
	}

    if(eqnSigmaPsi_ == "true")
    {
	    sigmaEps_ = 0.33 + 0.4*(tpProd_/(epsHat_*k_));
    }

    epsilonSafe_ = max(epsilon_, dimensionedScalar("minEps", epsilon_.dimensions(), 1.0e-15));

    if(eqncEp2_ == "true")
    {
        cEp2_ = cEp2con_ - 0.16*exp(-0.25*sqr(k_)/(nu()*epsilonSafe_));
    }
    else
    {
        cEp2_ =  cEp2con_;
    }

    cP1eqn_ = 2.0*(0.5+0.5*((tpProd_)/epsilonSafe_));


    //Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
       cEp1_*G*epsHat_ 
     + fvm::Sp(-1.0*cEp2_*epsHat_,epsilon_)
    );

    if(solveEps_ == "true")
    {
    epsEqn().relax();
    solve(epsEqn);
    bound(epsilon_,dimensionedScalar("minEps", epsilon_.dimensions(), 1.0e-15));
    }


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (

        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      + fvc::laplacian(nu()/cMu_,tpphi_)
      - fvm::Sp(epsilon_/k_,k_)
    );

    if(solveK_ == "true")
    {
    kEqn().relax();
    solve(kEqn);
    bound(k_,dimensionedScalar("minK", k_.dimensions(), 1.0e-15));
    }


    // K-related Terms
    kSqrt_ = sqrt(k_);
    bound(kSqrt_,dimensionedScalar("minKsqrt", kSqrt_.dimensions(), 3.16e-8));
    kSqrt_.correctBoundaryConditions();

    gradk_ = fvc::grad(k_);
    gradkSqrt_ = fvc::grad(kSqrt_);
    
    tpphisqrt_ = sqrt(tpphi_);
    volVectorField gradPhiSqrt("gradPhiSqrt", fvc::grad(tpphisqrt_));
       
    volScalarField GdK("GdK", G/k_);
    
    volVectorField gradPhiOverK("gradPhiOverK",fvc::grad(PhiOverK()));

    // Phi equation
    tmp<fvScalarMatrix> tpphiEqn
    (
        fvm::ddt(tpphi_)
      + fvm::div(phi_, tpphi_)
      + fvm::SuSp(-fvc::div(phi_), tpphi_)
      - fvm::laplacian(DphiEff(), tpphi_)
      ==
      // Pressure Strain
      // cPphi_*nutFrac()*(1.0 - Alpha())*((2.0*k_/3.0) - tpphi_)*epsHat_
      //+ fvm::SuSp( -0.06*((2.0*k_/(3.0*tpphi_)) - (tpphi_/tpphi_))*epsHat_,tpphi_ )
      (0.5 + cPphi_*(2.0*Alpha() - 1.0)*nutFrac())*tpphi_*epsHat_
      + cP2_*GdK*tpphi_
      // Pressure diffusion
      + fvm::Sp(-1.0*(cP4_+ cP2_)*GdK,tpphi_) 
      + (cP4_+ cP2_)*Alpha()*((tppsi_ & tppsi_)/((nut_*k_)*(1.0+cPw_/reTau())))*tpphi_
      // Dissipation
      + fvm::Sp(-2.0*Alpha()*epsHat_,tpphi_)
      + gT2_*fvm::Sp(-2.0*nu()*(gradPhiSqrt & gradPhiSqrt)/tpphi_,tpphi_)
      // Diffusion and Transition
      + gT1_*fvm::Sp(-1.0*(nu())*(gradk_ & gradPhiOverK)/tpphi_,tpphi_) 
      + gT3_*fvm::Sp(-1.0*(nut_)*(gradk_ & gradPhiOverK)/tpphi_,tpphi_)	       
      + cT_*(1.5*tpphi_-k_)*GdK*sqrt((nut_/nu()))
    );

    if(solvePhi_ == "true")
    {
    tpphiEqn().relax();
    solve(tpphiEqn);
    bound(tpphi_,dimensionedScalar("minTpphi", tpphi_.dimensions(), 1.0e-15));
    }

	// Phi output debug terms
	phiPressureStrain = (0.5 + cPphi_*(2.0*Alpha() - 1.0)*nutFrac())*tpphi_*epsHat_ + cP2_*GdK*tpphi_;
	phiPressureDiff = -1.0*(cP4_+ cP2_)*GdK*tpphi_ + (cP4_+ cP2_)*Alpha()*((tppsi_ & tppsi_)/((nut_*k_)*(1.0+cPw_/reTau())))*tpphi_;
	phiDiss = -2.0*Alpha()*epsHat_*tpphi_ - gT2_*2.0*nu()*(gradPhiSqrt & gradPhiSqrt);
	phiViscDiff = fvc::laplacian(nu(), tpphi_);
	phiTurbDiff = fvc::laplacian(sigmaPhi_*nut_, tpphi_);


    gradTpphi_ = fvc::grad(tpphi_);

	volTensorField gradPsiOverK("gradPsiOverK",fvc::grad(PsiOverK()));
	
	
	// Choice of production term for psi equation
	volScalarField psProd("psProd",GdK);
	
	if(psiProd_ == "psi")
	{
		volScalarField psProd("psProd",(tppsi_ & vorticity_)/k_);
	}


    // Psi Equation
    tmp<fvVectorMatrix> tppsiEqn
    (
        fvm::ddt(tppsi_)
      + fvm::div(phi_, tppsi_)
      + fvm::Sp(-fvc::div(phi_), tppsi_)
      - fvm::laplacian(DpsiEff(), tppsi_)

      ==
        (1.0 - cP2_)*tpphi_*vorticity_
      - (2*Alpha() - cP2_)*psProd*tppsi_ 
      + cMu_*(2*Alpha() - 1.0)*tpphi_*vorticity_
      - gT1_*2.0*(nu())*(gradk_ & gradPsiOverK) 
      - gT3_*2.0*(nut_)*(gradk_ & gradPsiOverK)     
	  + fvm::Sp(-1.0*(cP1_*nutFrac()*epsHat_*(1.0-Alpha())),tppsi_)   
      + fvm::Sp(-1.0*Alpha()*(epsilon_/k_),tppsi_)

      + gT2_*fvm::Sp(-2.0*nu()*(gradkSqrt_ & gradPhiSqrt)/(sqrt(k_*tpphi_)),tppsi_)
      + cT_*sqrt((nut_/nu()))*vorticity_*k_
    );

    if(solvePsi_ == "true")
    {
    tppsiEqn().relax();
    solve(tppsiEqn);
    }

    gradTppsi_ = fvc::grad(tppsi_);
    
    
    volScalarField psiZ(tppsi_.component(2));
	
	
	    // Calculate eddy viscosity
    if(solveNut_ == "true")
    {
        if(timeScaleEps_ == "epsilon" || timeScaleEps_ != "epsHat")
		{
            nut_ = cMu_*tpphi_*Ts();
            nut_.correctBoundaryConditions();
        }
        
        if(timeScaleEps_ == "epsHat")
		{
            nut_ = cMu_*tpphi_*TsEh();
            nut_.correctBoundaryConditions();
        }
        
    }
		
	
	Info<< "Maximum nut: " << gMax(nut_) << " Maximum K: " << gMax(k_) << " Maximum Epsilon: " << gMax(epsilon_) <<endl;
    Info<< "Maximum Phi: " << gMax(tpphi_) << " Maximum Psi_z: " << gMax(psiZ) << " Maximum Production: " << gMax(G) <<endl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
