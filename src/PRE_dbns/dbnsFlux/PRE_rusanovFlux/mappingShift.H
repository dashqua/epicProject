    // the mapping is arbitrary chosen beforehand and the following are
    // everything that composes the coefficient that shifts eigenvalues
    const scalar x = xyz[0];
    const scalar y = xyz[1];
    const scalar z = xyz[2];
    const scalar Tper = 2;
    const scalar pii = Foam::mathematicalConstant::pi;

// Get Jw
    scalar jw = ( 1 + pii/5*Foam::cos(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(2*pii*t/Tper)  ) * \
      ( 1 + pii/5*Foam::sin(pii*x/10)*Foam::cos(2*pii*y/15)*Foam::sin(4*pii*t/Tper)  )              - \
      ( 3*pii/20*Foam::cos(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(4*pii*t/Tper)   )              * \
      ( 4*pii/15*Foam::sin(pii*x/10)*Foam::cos(2*pii*y/15)*Foam::sin(2*pii*t/Tper)   )              ;

// Get Lambdaw
    scalar lambdaw = magSf / magSfOld;

// Get Cchi
    scalar cchi    = Foam::dot( UTilde / jw , Sf / magSf) / jw;
    
    //Step 6c: modification of eigenvalues with coefficient of step 6b
    lambda1 = (lambdaw/jw) * (lambda1 - cchi);
    lambda2 = (lambdaw/jw) * (lambda2 - cchi);
    lambda3 = (lambdaw/jw) * (lambda3 - cchi);

