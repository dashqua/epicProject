xlen = 1;
lc = 0.1;
Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, 0.1, lc};
Point(3) = {0, 0.1, 0.1, lc};
Point(4) = {0, 0.1, 0, lc};

Point(5) = {xlen, 0, 0, lc};
Point(6) = {xlen, 0, 0.1, lc};
Point(7) = {xlen, 0.1, 0.1, lc};
Point(8) = {xlen, 0.1, 0, lc};

Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line(5) = {8, 7};
Line(6) = {7, 6};
Line(7) = {6, 5};
Line(8) = {5, 8};
Line(9) = {8, 4};
Line(10) = {7, 3};
Line(11) = {6, 2};
Line(12) = {5, 1};

Curve Loop(1) = {2, 3, 4, 1};
Plane Surface(1) = {1};
Curve Loop(2) = {6, 7, 8, 5};
Plane Surface(2) = {2};
Curve Loop(3) = {10, 2, -11, -6};
Plane Surface(3) = {3};
Curve Loop(4) = {11, 3, -12, -7};
Plane Surface(4) = {4};
Curve Loop(5) = {12, 4, -9, -8};
Plane Surface(5) = {5};
Curve Loop(6) = {9, 1, -10, -5};
Plane Surface(6) = {6};

Transfinite Curve {1,2,3,4} = xlen/lc;
Transfinite Surface{1};
Recombine Surface{1};

Transfinite Curve {5,6,7,8} = xlen/lc;
Transfinite Surface{2};
Recombine Surface{2};

Transfinite Curve {9,10,11,12} = 10*xlen/lc;
Transfinite Surface{3};
Recombine Surface{3};

Transfinite Curve {13,14,15,16} = xlen/lc;
Transfinite Surface{4};
Recombine Surface{4};

Transfinite Curve {17,18,19,20} = xlen/lc;
Transfinite Surface{5};
Recombine Surface{5};

Transfinite Curve {21,22,23,24} = xlen/lc;
Transfinite Surface{6};
Recombine Surface{6};


Mesh.RecombineAll = 1;
Coherence ;