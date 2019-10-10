// Gmsh project created on Mon Oct 07 10:42:15 2019
x=10;
y=7.5;
//+
Point(1) = {-x, y, 0, 1.0};
//+
Point(2) = {-x, -y, 0, 1.0};
//+
Point(3) = {x, -y, 0, 1.0};
//+
Point(4) = {x, y, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1} = {1, 2, 3, 4};
//+
Transfinite Line {1, 3} = 16 Using Progression 1;
//+
Transfinite Line {4, 2} = 21 Using Progression 1;
//+
Mesh.RecombineAll = 1;//+
Extrude {0, 0, 1} {
  Line{4}; Line{1}; Line{2}; Line{3}; Surface{1}; Layers{1}; Recombine;
}
//+
Physical Volume("volume") = {1};
