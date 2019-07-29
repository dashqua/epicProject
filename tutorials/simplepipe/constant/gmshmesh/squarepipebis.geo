// tuto found in pipermail
//+ increase lc parameter for coarser mesh
lc = 4;

x = 10;
y = 7.5;
z = 1.;


//Point(0) = {0, 0, 0, lc};
Point(1) = {-x, -y, 0, lc};
Point(2) = {+x, -y, 0, lc};
Point(3) = {+x, +y, 0, lc};
Point(4) = {-x, +y, 0, lc};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Point(5) = {-x, -y, z, lc};
Point(6) = {+x, -y, z, lc};
Point(7) = {+x, +y, z, lc};
Point(8) = {-x, +y, z, lc};
Line(5) = {8, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
//+
Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};
//+
Line Loop(14) = {1, 2, 3 ,4};
Plane Surface(15) = {14};
Line Loop(16) = {5, 6, 7, 8};
Plane Surface(17) = {16};
Line Loop(18) = {4, 12, -8, -11};
Plane Surface(19) = {18};
Line Loop(20) = {2, 10, -6, -9};
Plane Surface(21) = {20};
Line Loop(22) = {3, 11, -7, -10};
Plane Surface(23) = {22};
Line Loop(24) = {1, 9, -5, -12};
Plane Surface(25) = {24};
//+
Transfinite Line{1:4,5:8} = 10;
Transfinite Line{9:12}    = 1;
Transfinite Surface{15} = {1,2,3,4};
Transfinite Surface{17} = {5,6,7,8};
Transfinite Surface{19} = {3,4,8,7};
Transfinite Surface{21} = {2,1,5,6};
Transfinite Surface{23} = {2,3,7,6};
Transfinite Surface{25} = {1,4,8,5};

Recombine Surface{19,21,23,25} = 45;

Physical Surface("inbound") = {25};
Physical Surface("right") = {17};
Physical Surface("left") = {15};
Physical Surface("outbound") = {23};
Physical Surface("top") = {19};
Physical Surface("bottom") = {21};
Surface Loop(1) = {17, 25, 15, 21, 23, 19};
Volume(1) = {1};
Physical Volume("domain") = {1};



Mesh 3;