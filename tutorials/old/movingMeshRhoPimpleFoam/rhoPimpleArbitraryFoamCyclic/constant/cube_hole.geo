x = 10;
y = 7.5;
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
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, .5} {
  Surface{1}; 
}
//+
Transfinite Surface {26} = {6, 10, 14, 5};
//+
Transfinite Surface {1} = {4, 3, 2, 1};
//+
Transfinite Surface {13} = {4, 1, 6, 5};
//+
Transfinite Surface {21} = {3, 2, 10, 14};
//+
Transfinite Surface {25} = {3, 14, 5, 4};
//+
Transfinite Surface {17} = {1, 6, 10, 2};
//+
Transfinite Line {12, 11, 16, 20} = 2 Using Progression 1;
//+
Transfinite Line {6, 4, 8, 2} = 21 Using Progression 1;
Transfinite Line {1, 7, 3, 9} = 16 Using Progression 1;
//+
Mesh.RecombineAll=1;
//+
Physical Surface("inlet") = {25};
//+
Physical Surface("outlet") = {17};
//+
Physical Surface("top") = {13};
//+
Physical Surface("bottom") = {21};
//+
Physical Surface("sides") = {1, 26};
//+
Physical Volume("volume") = {1};

