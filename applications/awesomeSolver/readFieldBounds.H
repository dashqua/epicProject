// Read field bounds
//dictionary fieldBounds = mesh.solutionDict().subDict("fieldBounds");

// Pressure bounds
dimensionedScalar pMin("pMin", p.dimensions(), SMALL);
dimensionedScalar pMax("pMax", p.dimensions(), ROOTVGREAT);
//fieldBounds.lookup(p.name() >> pMin.value() >> pMax.value());

//Temperature bounds
dimensionedScalar TMin("TMin", T.dimensions(), SMALL);
dimensionedScalar TMax("TMin", T.dimensions(), ROOTVGREAT);
//fieldBounds.lookup(T.name() >> TMin.value() >> TMax.value());
