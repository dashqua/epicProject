//+ increase lc parameter for coarser mesh
lc = 1.;

x = 10;
y = 7.5;
z = 1;

Point(1) = {-x, -y,   0, lc};
Point(2) = { x, -y,   0, lc};
Point(3) = { x,  y,   0, lc};
Point(4) = {-x,  y,   0, lc};
Point(5) = {-x, -y,   z, lc};
Point(6) = { x, -y,   z, lc};
Point(7) = { x,  y,   z, lc};
Point(8) = {-x,  y,   z, lc};

//+
Line(1) = {1, 5};
//+
Line(2) = {2, 6};
//+
Line(3) = {3, 7};
//+
Line(4) = {4, 8};
//+
Line(5) = {1, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 4};
//+
Line(8) = {4, 1};
//+
Line(9) = {5, 6};
//+
Line(10) = {6, 7};
//+
Line(11) = {7, 8};
//+
Line(12) = {8, 5};
//+
Line Loop(1) = {12, -1, -8, 4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {2, 10, -3, -6};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {11, -4, -7, 3};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {1, 9, -2, -5};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {6, 7, 8, 5};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {11, 12, 9, 10};
//+
Plane Surface(6) = {6};
//+

//+
Transfinite Line {4, 1, 3, 2} = 2 Using Progression 1;
//+
Transfinite Surface {6};
//+
Transfinite Surface {2};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6} = {8, 5, 6, 7};
//+
Transfinite Surface {5} = {4, 3, 2, 1};
//+
Transfinite Surface {1} = {4, 8, 5, 1};
//+
Transfinite Surface {4} = {6, 2, 1, 5};
//+
//Physical Surface("fixedWall") = {5, 6, 3, 4};
//+
Physical Surface("inlet") = {1};
//+
Physical Surface("outlet") = {2};
//+
//+
Surface Loop(1) = {6, 3, 1, 4, 2, 5};
//+
Volume(1) = {1};
//+
Physical Volume("box") = {1};
//+
Transfinite Surface {3};
//+
//+
Recombine Surface{1,2,3,4,5,6};
//+
//+
//+ UNCOMMENT THE FOLLOWING AT YOUR OWN RISK
//Mesh.Algorithm = 8;
//Mesh.Algorithm3D = 6;
//Mesh.RecombineAll = 1;
//Mesh.Recombine3DAll = 1;
//Mesh.SubdivisionAlgorithm = 2;
//+
Physical Surface("top") = {3};
//+
Physical Surface("bottom") = {4};
//+
Physical Surface("left") = {5};
//+
Physical Surface("right") = {6};
