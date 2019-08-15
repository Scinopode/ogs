// Grid number
GridNum = 41;

// Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.75, 0, 0, 1.0};
Point(3) = {2, 0, 0, 1.0};
Point(4) = {2, 3, 0, 1.0};
Point(5) = {2, 5, 0, 1.0};
Point(6) = {0.75, 5, 0, 1.0};
Point(7) = {0.75, 3, 0, 1.0};
Point(8) = {0, 3, 0, 1.0};

// Lines
Line(1) = {1, 2}; Transfinite Line {1} = 16 Using Progression 1;
Line(2) = {2, 7}; Transfinite Line {2} = 61 Using Progression 1;
Line(3) = {7, 8}; Transfinite Line {3} = 16 Using Progression 1;
Line(4) = {8, 1}; Transfinite Line {4} = 61 Using Progression 1;
Line(5) = {2, 3}; Transfinite Line {5} = 26 Using Progression 1;
Line(6) = {3, 4}; Transfinite Line {6} = 61 Using Progression 1;
Line(7) = {4, 7}; Transfinite Line {7} = 26 Using Progression 1;
Line(8) = {4, 5}; Transfinite Line {8} = 41 Using Progression 1;
Line(9) = {5, 6}; Transfinite Line {9} = 26 Using Progression 1;
Line(10) = {6, 7}; Transfinite Line {10} = 41 Using Progression 1;

// Surfaces
Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {11};
Line Loop(13) = {5, 6, 7, -2};
Plane Surface(14) = {13};
Line Loop(15) = {10, -7, 8, 9};
Plane Surface(16) = {15};

// For structured mapping
Transfinite Surface {12};
Transfinite Surface {14};
Transfinite Surface {16};
Recombine Surface {12};
Recombine Surface {14};
Recombine Surface {16};

// Physical groups
Physical Surface(47) = {12, 14, 16};
