x=10;
y=7.5;
z=0.5;
//+
Point(1) = {-x, y, 0, 1.0};
//+
Point(2) = {x, y, 0, 1.0};
//+
Point(3) = {x, -y, 0, 1.0};
//+
Point(4) = {-x, -y, 0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Extrude {0, 0, z} {
  Line{1}; Line{4}; Line{3}; Line{2}; 
}
//+
Line Loop(1) = {5, 17, 13, 9};
//+
Plane Surface(21) = {1};
//+
Line Loop(2) = {4, 1, 2, 3};
//+
Plane Surface(22) = {2};
//+
Physical Surface("inlet") = {8};
//+
Physical Surface("outlet") = {16};
//+
Physical Surface("top") = {12};
//+
Physical Surface("bottom") = {20};
//+
Physical Surface("sides") = {22, 21};
//+
Surface Loop(1) = {21, 8, 22, 12, 16, 20};
//+
Volume(1) = {1};
//+
Physical Volume("volume") = {1};
//+
Transfinite Surface {21} = {5, 6, 8, 7};
//+
Transfinite Surface {22} = {2, 1, 4, 3};
