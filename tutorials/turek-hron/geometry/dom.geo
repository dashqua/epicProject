// Gmsh project created on Thu Oct 10 14:42:11 2019
//+ dim of dom
xmin=0;
ymin=0;
xmax=5;
ymax=1;
zmax=.1;
//+ big box
Point(1) = {xmin, ymax, 0, 1.0};
//+
Point(2) = {xmin, ymin, 0, 1.0};
//+
Point(3) = {xmax, ymin, -0, 1.0};
//+
Point(4) = {xmax, ymax, 0, 1.0};
//+
Point(20) = {(xmax-xmin)/4, ymax, -0, 1.0};
//+
Point(21) = {(xmax-xmin)/4, ymin, -0, 1.0};
//+
Point(22) = {3*(xmax-xmin)/4, ymin, 0, 1.0};
//+
Point(23) = {3*(xmax-xmin)/4, ymax, 0, 1.0};
//+ position et taille languette
tl=.2;
xl=(xmax-xmin)/3.;
yl=(ymax-ymin)/2.;
ll=(xmax-xmin)/4.;
//+ cercle
r=.15;
Point(5) = {xl, yl, 0, 1.0};
//+
Point(6) = {xl, yl+tl/2., 0, 1.0};
//+
Point(7) = {xl, yl-tl/2., 0, 1.0};
//+
Point(8) = {xl-r, yl+tl/2., 0, 1.0};
//+
Point(9) = {xl+r, yl+tl/2., 0, 1.0};
//+
Point(10) = {xl+r, yl-tl/2., 0, 1.0};
//+
Point(11) = {xl-r, yl-tl/2., 0, 1.0};
//+
Point(12) = {xl+r+ll, yl-tl/2., 0, 1.0};
//+
Point(13) = {xl+r+ll, yl+tl/2., 0, 1.0};
//+
Point(24) = {xl, yl+tl/2.+r, 0, 1.0};
//+
Point(25) = {xl, yl-tl/2.-r, 0, 1.0};
//+
Point(26) = {(xmax-xmin)/4, 0.8*ymax, 0.0, 1.0};
//+
Point(27) = {(xmax-xmin)/4, 0.2*ymax, 0.0, 1.0};
//+
Point(28) = {3*(xmax-xmin)/4, 0.8*ymax, 0.0, 1.0};
//+
Point(29) = {3*(xmax-xmin)/4, 0.2*ymax, 0.0, 1.0};
//+
Point(30) = {xmin, 0.8*ymax, 0.0, 1.0};
//+
Point(31) = {xmax, 0.8*ymax, 0.0, 1.0};
//+
Point(32) = {xmin, 0.2*ymax, 0.0, 1.0};
//+
Point(33) = {xmax, 0.2*ymax, 0.0, 1.0};
//+
Line(1) = {1, 20};
//+
Line(3) = {26, 27};
//+
Line(5) = {21, 2};
//+
Line(7) = {20, 23};
//+
Line(9) = {28, 29};
//+
Line(11) = {22, 21};
//+
Line(12) = {22, 3};
//+
Line(14) = {4, 23};
//+
Line(15) = {28, 26};
//+
Line(16) = {27, 29};
//+
Line(17) = {10, 9};
//+
Line(18) = {9, 13};
//+
Line(19) = {13, 12};
//+
Line(20) = {12, 10};
//+
Line(21) = {8, 11};
//+
Circle(22) = {8, 6, 24};
//+
Circle(23) = {24, 6, 9};
//+
Circle(24) = {10, 7, 25};
//+
Circle(25) = {25, 7, 11};
//+
Line(26) = {1, 30};
//+
Line(27) = {30, 32};
//+
Line(28) = {2, 32};
//+
Line(29) = {3, 33};
//+
Line(30) = {33, 31};
//+
Line(31) = {31, 4};
//+
Line(32) = {31, 28};
//+
Line(33) = {29, 33};
//+
Line(34) = {29, 22};
//+
Line(35) = {28, 23};
//+
Line(36) = {26, 20};
//+
Line(37) = {26, 30};
//+
Line(38) = {32, 27};
//+
Line(39) = {27, 21};
//+
Transfinite Line {26, 36, 35, 31, 29, 34, 39, 28} = 3 Using Progression 1;
//+
Transfinite Line {27, 3, 9, 30} = 7 Using Progression 1;
//+
Transfinite Line {1, 37, 38, 5, 14, 32, 33, 12} = 16 Using Progression 1;
//+
Transfinite Line {7, 15, 16, 11} = 31 Using Progression 1;
//+
Line Loop(1) = {1, -36, 37, -26};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {37, 27, 38, -3};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {7, -35, 15, 36};
//+
Plane Surface(3) = {3};
//+
Transfinite Surface {1} = {1, 30, 26, 20};
//+
Transfinite Surface {2} = {30, 32, 27, 26};
//+
Transfinite Surface {3} = {26, 28, 23, 20};
//+
Line Loop(4) = {15, 3, 16, -9};
//+
Line Loop(5) = {23, 18, 19, 20, 24, 25, -21, 22};
//+
Plane Surface(4) = {4, 5};
//+
Line Loop(6) = {38, 39, 5, 28};
//+
Plane Surface(5) = {6};
//+
Transfinite Surface {5} = {32, 2, 21, 27};
//+
Line Loop(7) = {16, 34, 11, -39};
//+
Plane Surface(6) = {7};
//+
Line Loop(8) = {33, -29, -12, -34};
//+
Plane Surface(7) = {8};
//+
Line Loop(9) = {30, 32, 9, 33};
//+
Plane Surface(8) = {9};
//+
Line Loop(10) = {31, 14, -35, -32};
//+
Plane Surface(9) = {10};
//+
Transfinite Surface {6} = {27, 21, 22, 29};
//+
Transfinite Surface {7} = {29, 22, 3, 33};
//+
Transfinite Surface {8} = {28, 29, 33, 31};
//+
Transfinite Surface {9} = {23, 28, 31, 4};
//+
Transfinite Line {18, 20} = 16 Using Progression 1;
//+
Transfinite Line {21, 19} = 3 Using Progression 1;
//+
Extrude {0, 0, .5} {
  Surface{1}; Surface{2}; Surface{4}; Surface{3}; Surface{6}; Surface{5}; Surface{7}; Surface{8}; Surface{9}; Line{17}; Layers{1}; Recombine;
}
//+
Physical Surface("emptyPatches") = {83, 61, 211, 189, 145, 167, 277, 255, 233, 1, 2, 5, 6, 4, 3, 9, 8, 7};
//+
Physical Surface("inlet") = {60, 74, 210};
//+
Physical Surface("outlet") = {264, 242, 224};
//+
Physical Surface("top") = {48, 154, 268};
//+
Physical Surface("bottom") = {206, 184, 228};
//+
Physical Surface("fixed_hole") = {116, 120, 124, 144, 128};
//+
Physical Surface("bar") = {140, 136, 132};
//+
Physical Volume("domain") = {4, 1, 2, 6, 5, 8, 9, 7, 3};
//+
Line Loop(11) = {95, -278, 93, 94};
//+
Plane Surface(282) = {11};
//+
Line Loop(12) = {20, 17, 18, 19};
//+
Plane Surface(283) = {12};
//+
Surface Loop(1) = {281, 283, 282, 136, 140, 132};
//+
Volume(10) = {1};
//+
Physical Volume("domBar") = {10};
//+
//+
Mesh 2;
RefineMesh;
