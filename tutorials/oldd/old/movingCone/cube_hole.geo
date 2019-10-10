//+
lc1=0.1;
//+
Point(1) = {-1, 1, 0, lc1};
//+
Point(2) = {-1, -1, 0, lc1};
//+
Point(3) = {1, -1, 0, lc1};
//+
Point(4) = {1, 1, 0, lc1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
lc2=.1;
//+
Point(5) = {-0.1, -0.1, 0, lc2};
//+
Point(6) = {0.1, -0.1, -0, lc2};
//+
Point(7) = {0.1, 0.1, -0, lc2};
//+
Point(8) = {-0.1, 0.1, 0, lc2};
//+
Line(5) = {8, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 8};

//+
Line Loop(1) = {5, 6, 7, 8};
//+
Line Loop(2) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1, 2};
//+
Extrude {0, 0, 0.1} {
  Line{5}; Line{8}; Line{7}; Line{6}; Line{2}; Line{1}; Line{4}; Line{3}; 
}
//+
Line Loop(3) = {13, 9, 21, 17};
//+
Line Loop(4) = {33, 29, 25, 37};
//+
Plane Surface(41) = {3, 4};
//+
Transfinite Surface {32};
//+
Transfinite Surface {36};
//+
Transfinite Surface {40};
//+
Transfinite Surface {28};

//+
Transfinite Surface {16};
//+
Transfinite Surface {12};
//+
Transfinite Surface {24};
//+
Transfinite Surface {20};
