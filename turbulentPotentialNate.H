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

Class
    Foam::incompressible::RASModels::turbulentPotentialNate

Description
    Implementation of the Turbulent Potential turbulence model for incompressible
    flows.

    Turbulence model described in:
    @verbatim
        Wang, C. and J. Blair Perot
        "Prediction of turbulent transition in boundary layers using the turbulent potential model"
        Journal of Turbulence,
        2002
    @endverbatim

    The default model coefficients correspond to the following:
    @verbatim
        turbulentPotentialNateCoeffs
        {

        }
    @endverbatim

SourceFiles
    turbulentPotentialNate.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentPotentialNate_H
#define turbulentPotentialNate_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class turbulentPotentialNate Declaration
\*---------------------------------------------------------------------------*/

class turbulentPotentialNate
:
    public RASModel
{
    // Private data

        // Model coefficients
            dimensionedScalar cEp1_;
            dimensionedScalar cEp2con_;
            dimensionedScalar cEp3_;
            dimensionedScalar cD1_;
            dimensionedScalar cD2_;
            dimensionedScalar cKp_;
            dimensionedScalar cVv1_;
            dimensionedScalar cTv1_;
            dimensionedScalar cP1_;
            dimensionedScalar cP2_;
            dimensionedScalar cP3_;
            dimensionedScalar cP4_;
            dimensionedScalar cPphi_;
            dimensionedScalar cMu_;
            dimensionedScalar cT_;
            dimensionedScalar cPr_;
            dimensionedScalar cEhm_;
            dimensionedScalar cEhR_;
			dimensionedScalar gT1_;
			dimensionedScalar gT2_;
			dimensionedScalar gT3_;
			dimensionedScalar cNF_;
			dimensionedScalar cPw_;
            dimensionedScalar sigmaKInit_;
            dimensionedScalar sigmaEpsInit_;
            dimensionedScalar sigmaEpsVisc_;
            dimensionedScalar sigmaPhiInit_;
            dimensionedScalar sigmaPsiInit_;
            dimensionedScalar psiNuFrac_;
            word              solveK_;
            word              solveEps_;
            word              solvePsi_;
            word              solvePhi_;
            word              solveNut_;
            word              eqnSigmaK_;
            word              eqnSigmaEps_;
            word              eqnSigmaPhi_;
            word              eqnSigmaPsi_;
            word              eqncEp2_;
       	    word              eqnEpsHat_;
            word              timeScaleEps_;
			word              prodType_;
			word              debugWrite_;
			word              tslimiter_;
			word              psiProd_;

        //- Wall distance field
        //  Note: different to wall distance in parent RASModel
            wallDist y_;

        // Fields

            volScalarField k_;
            volVectorField gradk_;
            volScalarField epsilon_;
            volScalarField nut_;
            volScalarField nutNorm_;
            volScalarField tpphi_;
            volScalarField tpphisqrt_;
            volVectorField vorticity_;
			volVectorField phis_;
            volVectorField vorticityTmp_;
            volVectorField ivorticity_;
            volVectorField tppsi_;
            volTensorField uGrad_;
            volScalarField epsHat_;
            volScalarField eHrC_;
			volScalarField kol_;
			volScalarField kSafe_;
			volScalarField kSqrt_;
			volVectorField gradkSqrt_;
			volScalarField nutSafe_;
			volScalarField epsilonSafe_;
			volScalarField sigmaK_;
			volScalarField sigmaEps_;
			volScalarField sigmaPhi_;
			volScalarField sigmaPsi_;
			volScalarField cEp2_;
			volScalarField tpProd_;
			volScalarField cP1eqn_;
			volScalarField dimRat_;
			volVectorField gradTpphi_;
			volTensorField gradTppsi_;
			volScalarField tpProdSqr_;
			volScalarField tpProd3d_;
			volScalarField phiPressureStrain;
			volScalarField phiPressureDiff;
			volScalarField phiDiss;
			volScalarField phiViscDiff;
			volScalarField phiTurbDiff;


    // Private member functions

    //- Return time scale, Ts using epsilon
    tmp<volScalarField> Ts() const;
    
    //- Return time scale, TsEh using epsHat
    tmp<volScalarField> TsEh() const;


public:

    //- Runtime type information
    TypeName("turbulentPotentialNate");


    // Constructors

        //- Construct from components
        turbulentPotentialNate
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor

        virtual ~turbulentPotentialNate()
        {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }
        
        //- Return the kolmogorov time scale
        tmp<volScalarField> minTS() const
		{
            return tmp<volScalarField>
            (
				new volScalarField("minTS",6.0*sqrt(nu()/(epsilon_ + epsilonSmall_)))
			);
		}

        //- Return phiS
        tmp<volSymmTensorField> phiS() const
        {
            return tmp<volSymmTensorField>
            (
                new volSymmTensorField("phiS", tpphi_*symm(uGrad_))
            );
        }
		
	//- Return divergence of phiS
        tmp<volVectorField> divphiS() const
        {
            return tmp<volVectorField>
            (
                new volVectorField("divphiS", fvc::div(phiS()))
            );
        }
		
	//- Return s
        tmp<volVectorField> sreal() const
        {
            return tmp<volVectorField>
            (
                new volVectorField("sreal", phis_/tpphi_)
            );
        }


        //- Return the phi/k
        tmp<volScalarField> PhiOverK() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("PhiOverK", tpphi_/(k_ + k0_))
            );
        }

		//- Return the psi/k
        tmp<volVectorField> PsiOverK() const
        {
            return tmp<volVectorField>
            (
                new volVectorField("PsiOverK", tppsi_/(k_ + k0_))
            );
        }

        //- Return the gradient of K
        tmp<volScalarField> gradsqK() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("gradsqK", mag(fvc::grad(sqrt(k_))))
            );
        }

        //- Return the gradient of Phi
        tmp<volVectorField> gradtpphi() const
        {
            return tmp<volVectorField>
            (
                new volVectorField("gradtpphi", fvc::grad(tpphi_))
            );
        }

        //- Return gradient of sqrt of Phi
        tmp<volVectorField> gradtpphisqrt() const
        {
            return tmp<volVectorField>
            (
                new volVectorField("gradtpphisqrt", fvc::grad(tpphisqrt_))
            );
        }

        //- Return the local turbulent reynolds number
        tmp<volScalarField> reTau() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("reTau", k_*k_/(nu()*(epsilon_ + epsilonSmall_)))
            );
        }


        //- Return the epsilon hat
        tmp<volScalarField> EpsHat() const
        {
            return epsHat_;
        }

        //- Return the nut fraction
        tmp<volScalarField> nutFrac() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("nutFrac", nut()/(nut()+cNF_*nu()))
            );
        }


        //- Return Alpha
        tmp<volScalarField> Alpha() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("Alpha", 1.0/(1.0 + 1.5*PhiOverK()))
            );
        }

        //- Return cEp2
        tmp<volScalarField> cEp2() const
        {
                return cEp2_;
        }

        //- Return the effective diffusivity for k
        virtual tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_*sigmaK_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_*sigmaEps_ + nu()*sigmaEpsVisc_)
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DphiEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DphiEff", nut_*sigmaPhi_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DpsiEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DpsiEff", nut_*sigmaPsi_ + psiNuFrac_*nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DEff", nut_ + nu())
            );
        }

        //- Return the Reynolds number
        tmp<volScalarField> tpReynolds() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("tpReynolds", k_*k_/((epsilon_ + epsilonSmall_)*nu()))
            );
        }

        //- Return the vorticity
        virtual tmp<volVectorField> vorticity() const
        {
            return vorticity_;
        }

        //- Return the production
        virtual tmp<volScalarField> psiproduction() const
        {
            return tpProd_;
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the turbulent potential psi
        virtual tmp<volVectorField> psiReal() const
        {
            return tppsi_;
        }

        //- Return the turbulent potential psi
        virtual tmp<volScalarField> phiReal() const
        {
            return tpphi_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
