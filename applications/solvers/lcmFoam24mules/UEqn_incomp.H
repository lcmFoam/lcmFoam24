
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
        + turbulence.divDevTau(rho, U)
//     ==
//        fvOptions(rho, U)
    );

    
    UEqn.relax();
    pZones.addResistance(UEqn);
    fvConstraints.constrain(UEqn);


        
        if (pimple.momentumPredictor())
        {

            //pZones.addResistance(UEqn);
            //
            //fvConstraints.constrain(UEqn);

            solve(UEqn == -fvc::grad(p));

            fvConstraints.constrain(U);
    
        }


