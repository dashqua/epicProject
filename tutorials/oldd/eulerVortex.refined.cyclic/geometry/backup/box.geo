x=10;
y=7.5;
extr=(x+y)/20;
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
Extrude {0, 0, extr} {
  Line{4}; Line{1}; Line{2}; Line{3}; 
}
//+
Line Loop(1) = {5, 9, 13, 17};
//+
Plane Surface(21) = {1};
//+
Line Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(22) = {2};
//+
Transfinite Surface {21} = {6, 8, 10, 5};
//+
Transfinite Surface {22} = {1, 2, 3, 4};
//+
Transfinite Line {9, 1, 17, 3} = 16 Using Progression 1;
//+
Transfinite Line {4, 5, 2, 13} = 21 Using Progression 1;
//+
Transfinite Line {7, 11, 15, 6} = 2 Using Progression 1;
//+
Transfinite Surface {12} = {1, 2, 8, 6};
//+
Transfinite Surface {8} = {1, 6, 5, 4};
//+
Transfinite Surface {20} = {5, 4, 3, 10};
//+
Transfinite Surface {16} = {8, 2, 3, 10};
//+
Transfinite Line {7, 6, 11, 15} = 2 Using Progression 1;

//+
Mesh.RecombineAll = 1;

//+
Physical Surface("emptyPatches") = {21, 22};
//+
Physical Surface("inlet") = {12};
//+
Physical Surface("outlet") = {20};
//+
Physical Surface("top") = {8};
//+
Physical Surface("bottom") = {16};
//+
Surface Loop(1) = {21, 8, 22, 20, 16, 12};
//+
Volume(1) = {1};
//+
Physical Volume("volume") = {1};
//+

