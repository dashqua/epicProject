xmin=-10;
xmax=10;
ymin=-7.5;
ymax=7.5;
z=1;
//+
Point(1) = {xmin, ymin, 0, 1.0};
//+
Point(2) = {xmax, ymin, 0, 1.0};
//+
Point(3) = {xmin, ymax, 0, 1.0};
//+
Point(4) = {xmax, ymax, 0, 1.0};
//+
Line(1) = {3, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 4};
//+
Line(4) = {4, 3};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1} = {3, 1, 2, 4};
//+
Transfinite Line {1, 4, 3, 2} = 41 Using Progression 1;
//+
Extrude {0, 0, z} {
  Surface{1}; Line{4}; Line{1}; Line{2}; Line{3}; Layers{1}; Recombine;
}
//+
Mesh.RecombineAll = 1 ;
//+
Physical Surface("inlet") = {17};
//+
Physical Surface("outlet") = {25};
//+
Physical Surface("top") = {13};
//+
Physical Surface("bottom") = {21};
//+
Physical Surface("emptyPatches") = {26, 1};
//+
Physical Volume("domain") = {1};
