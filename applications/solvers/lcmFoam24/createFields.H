#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "createPhi.H"
 
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);


Info<< "Reading phaseProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

autoPtr<twoPhaseChangeModel> phaseChangePtr
(
    twoPhaseChangeModel::New(mixture)
);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.rho()
);
rho.oldTime();

pressureReference pressureReference(p, p, pimple.dict());
mesh.schemes().setFluxRequired(p.name());


volScalarField& alpha1(mixture.alpha1());
mesh.schemes().setFluxRequired(alpha1.name());


#include "createAlphaFluxes.H"


// Construct incompressible turbulence model
incompressibleInterPhaseTransportModel turbulence
(
    U,
    phi,
    alphaPhi1,
    mixture
);





// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"







//Fields for compressibleLcmFoam

dimensionedScalar rho1_val ("rho1_val",dimensionSet(1, -3, 0, 0, 0, 0 ,0), 1110.0);
dimensionedScalar nu1_val ("nu1_val",dimensionSet(0, 2, -1, 0, 0, 0 ,0), 297e-06);
dimensionedScalar c_sound_resin ("c_sound_resin",dimensionSet(0, 1, -1, 0, 0, 0 ,0), 1470.0);
dimensionedScalar p1_0 ("p1_0",dimensionSet(1, -1, -2, 0, 0, 0 ,0), 100000);
dimensionedScalar rho1_0 ("rho1_0",dimensionSet(1, -3, 0, 0, 0, 0 ,0), 1110.0);

dimensionedScalar rho2_val ("rho2_val",dimensionSet(1, -3, 0, 0, 0, 0 ,0), 1.225);
dimensionedScalar nu2_val ("nu2_val",dimensionSet(0, 2, -1, 0, 0, 0 ,0), 14.6e-06);
dimensionedScalar R_s ("R_s",dimensionSet(0, 2, -2, -1, 0, 0 ,0), 287.0);
dimensionedScalar T_val ("T_val",dimensionSet(0, 0, 0, 1, 0, 0 ,0), 293.0);
dimensionedScalar p2_0 ("p2_0",dimensionSet(1, -1, -2, 0, 0, 0 ,0), 100000);
dimensionedScalar rho2_0 ("rho2_0",dimensionSet(1, -3, 0, 0, 0, 0 ,0), 1.225);

dimensionedScalar mu1_val=rho1_val*nu1_val;
dimensionedScalar mu2_val=rho2_val*nu2_val;

dimensionedScalar dummy1("dummy1",dimensionSet(-1,1,2,0,0,0,0),scalar(1)); 
dimensionedScalar pminval("pminval",dimensionSet(1,-1,-2,0,0,0,0),scalar(1000)); 

volVectorField U_
(
    IOobject
    (
        "U_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U
);
surfaceScalarField Ufaces_
(
    IOobject
    (
        "Ufaces_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U_) & mesh.Sf()
);
volScalarField alpha1_
(
    IOobject
    (
        "alpha1_",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField alpha2_
(
    IOobject
    (
        "alpha2_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1-alpha1_
);
volScalarField psi1_
(
    IOobject
    (
        "psi1_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    1/pow(c_sound_resin,2)
);
volScalarField psi2_
(
    IOobject
    (
        "psi2_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    1/(R_s*T_val)
    //0.01*psi1_
);
volScalarField rho1_
(
    IOobject
    (
        "rho1_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //mesh, rho1_val
    rho1_0+psi1_*(p-p1_0)
);
volScalarField rho2_
(
    IOobject
    (
        "rho2_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //mesh, rho2_val
    max(pminval,p)/(R_s*T_val)
);
volScalarField rhoeff_
(
    IOobject
    (
        "rhoeff_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1_*rho1_+alpha2_*rho2_
);
volScalarField mueff_
(
    IOobject
    (
        "mueff_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alpha1_*mu1_val+alpha2_*mu2_val)  //volume averaged
);
/*
*/






