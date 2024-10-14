{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    adjustPhi(phiHbyA, U, p);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAUf, MRF);




 while (pimple.correctNonOrthogonal())
 {
     // Pressure corrector
     fvScalarMatrix pEqn
     (
         fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
     );
 
     pEqn.setReference(pRefCell, pRefValue);
 
     pEqn.solve();
 
     if (pimple.finalNonOrthogonalIter())
     {
         phi = phiHbyA - pEqn.flux();
     }
 }
 
 #include "continuityErrs.H"
 
 p.relax();
 
 // Momentum corrector
 U = HbyA - rAU*fvc::grad(p);
 //U.correctBoundaryConditions();


  
}
