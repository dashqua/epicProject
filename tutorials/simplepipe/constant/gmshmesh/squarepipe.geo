//+ increase lc parameter for coarser mesh
lc = 4.;


//+
//+
//+
//+ UNCOMMENT THE FOLLOWING AT YOUR OWN RISK
//Mesh.Algorithm = 8;
//Mesh.Algorithm3D = 6;
//Mesh.RecombineAll = 1;
//Mesh.Recombine3DAll = 1;
//Mesh.SubdivisionAlgorithm = 2;
//+


x = 10;
y = 1;
z = .1;



Point(0) = {x/2, y/2, 0, 1.0};
Point(1) = {0, 0, 0, 1.0};
Point(2) = {x, 0, 0, 1.0};
Point(3) = {x, y, 0, 1.0};
Point(4) = {0, y, 0, 1.0};

Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

Point(5) = {0, 0, z, 1.0};
Point(6) = {x, 0, z, 1.0};
Point(7) = {x, y, z, 1.0};
Point(8) = {0, y, z, 1.0};

Line(5) = {4, 8};
Line(6) = {3, 7};
Line(7) = {2, 6};
Line(8) = {1, 5};
Line(9) = {5, 8};
Line(10) = {8, 7};
Line(11) = {7, 6};
Line(12) = {6, 5};


Line Loop(1) = {5, -9, -8, -1};
Plane Surface(1) = {1};
Line Loop(2) = {4, 5, 10, -6};
Plane Surface(2) = {2};
Line Loop(3) = {6, 11, -7, 3};
Plane Surface(3) = {3};
Line Loop(4) = {8, -12, -7, -2};
Plane Surface(4) = {4};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(5) = {5};
Line Loop(6) = {12, 9, 10, 11};
Plane Surface(6) = {6};

//Physical Surface("side") = {1, 2, 3, 4};
//Physical Surface("inbound") = {5};
//Physical Surface("outbound") = {6};
//Surface Loop(1) = {1, 2, 5, 4, 6, 3};
//Volume(1) = {1};
//Physical Volume("pipe") = {1};
//+
Physical Surface("inbound") = {1};
//+
Physical Surface("side") = {6, 2, 5, 4};
//+
Physical Surface("outbound") = {3};
