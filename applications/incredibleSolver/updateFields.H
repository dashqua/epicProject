{
    // Compute U
    U = rhoU/rho;
    U.correctBoundaryConditions();

    // Calculate enthalpy from rhoE
    const volScalarField Cp = thermo->Cp();
    const volScalarField Cv = thermo->Cv();

//     h = rhoE/rho - 0.5*magSqr(U) + p/rho;
    // Alternative formulation, Felipe Alves Portela TU Delft
    h = Cp/Cv*(rhoE/rho - 0.5*magSqr(U));
    h.correctBoundaryConditions();

    // Bound enthalpy
    dimensionedScalar CpMin = min(Cp);
    dimensionedScalar CpMax = max(Cp);

    dimensionedScalar hMin_ = CpMin*TMin; dimensionedScalar& hMin = hMin_;
    dimensionedScalar hMax_ = CpMax*TMax; dimensionedScalar& hMax = hMax_;

    //bound(h, hMin);
    reBoundMinMax(h, hMin, hMax);
        
    // Correct thermodynamics
    thermo->correct();

    // Bound density
    //bound(rho, rhoMin);
    reBoundMinMax(rho, rhoMin, rhoMax);

    // Compute p from rho
    p = rho*(thermo->Cp() - thermo->Cv())*T;
    p.correctBoundaryConditions();

    // Bound pressure
    //bound(p, pMin);
    reBoundMinMax(p, pMin, pMax);
}
