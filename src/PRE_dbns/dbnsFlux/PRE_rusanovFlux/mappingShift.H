    // the mapping is arbitrary chosen beforehand and the following are
    // everything that composes the coefficient that shifts eigenvalues

    const scalar x = xyz[0];
    const scalar y = xyz[1];
    const scalar z = xyz[2];
    const scalar Tper = 2;
    const scalar pii = Foam::mathematicalConstant::pi;


    scalar detfw = ( 1 + pii/5*Foam::cos(pii*x/10)*Foam::sin(2*pii*y/15)*Foam::sin(2*pii*t/Tper)  );

    scalar lambdaw = 1;
    scalar jw      = detfw;
    scalar cchi    = 0;
    
    //Step 6c: modification of eigenvalues with coefficient of step 6b
    lambda1 = (lambdaw/jw) * (lambda1 - cchi);
    lambda2 = (lambdaw/jw) * (lambda2 - cchi);
    lambda3 = (lambdaw/jw) * (lambda3 - cchi);

