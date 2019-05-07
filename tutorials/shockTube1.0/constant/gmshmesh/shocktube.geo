thick=.5;

lf=1;
lc=lf/thick;

Point(1) = {-1, -1, 0, 1.0};
Point(2) = {1, -1, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {-1, 1, 0, 1.0};
Point(5) = {-1, -1, thick, 1.0};
Point(6) = {1, -1, thick, 1.0};
Point(7) = {1, 1, thick, 1.0};
Point(8) = {-1, 1, thick, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {8, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 8};
//+
Line(9) = {4, 8};
//+
Line(10) = {1, 5};
//+
Line(11) = {2, 6};
//+
Line(12) = {3, 7};
//+
Curve Loop(1) = {3, 9, -8, -12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 5, -10, -4};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, 11, -6, -10};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {11, 7, -12, -2};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {8, 5, 6, 7};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {3, 4, 1, 2};
//+
Plane Surface(6) = {6};
//+
Physical Surface("top") = {1};
//+
Physical Surface("bottom") = {3};
//+
Physical Surface("inlet") = {2};
//+
Physical Surface("oulet") = {4};
//+
Physical Surface("above") = {5};
//+
Physical Surface("below") = {6};
//+
Surface Loop(1) = {5, 1, 6, 2, 3, 4};
//+
Volume(1) = {1};
//+
Physical Volume("box") = {1};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6};
//+
Transfinite Surface {2};
//+
Transfinite Surface {1};
//+
Transfinite Surface {4};
//+
Transfinite Surface {3};
