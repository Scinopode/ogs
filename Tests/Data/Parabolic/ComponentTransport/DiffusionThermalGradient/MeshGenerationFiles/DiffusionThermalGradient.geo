// Grid number
GridNum = 41;

// Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.038, 0, 0, 1.0};
Point(3) = {0.038, 0.5, 0, 1.0};
Point(4) = {0, 0.5, 0, 1.0};
Point(5) = {0.2, 0, 0, 1.0};
Point(6) = {0.2, 0.5, 0, 1.0};
Point(7) = {0.2, 5.5, 0, 1.0};
Point(8) = {0.038, 5.5, 0, 1.0};
Point(9) = {0.2, 6, 0, 1.0};
Point(10) = {3, 0, 0, 1.0};
Point(11) = {3, 0.5, 0, 1.0};
Point(12) = {3, 5.5, 0, 1.0};
Point(13) = {3, 6, 0, 1.0};

// Lines
Line(1) = {1, 2}; Transfinite Line {1} = 6 Using Progression 1;
Line(2) = {2, 3}; Transfinite Line {2} = 51 Using Progression 1;
Line(3) = {3, 4}; Transfinite Line {3} = 6 Using Progression 1;
Line(4) = {4, 1}; Transfinite Line {4} = 51 Using Progression 1;
Line(5) = {2, 5}; Transfinite Line {5} = 21 Using Progression 1;
Line(6) = {5, 6}; Transfinite Line {6} = 51 Using Progression 1;
Line(7) = {6, 3}; Transfinite Line {7} = 21 Using Progression 1;
Line(8) = {6, 7}; Transfinite Line {8} = 501 Using Progression 1;
Line(9) = {7, 8}; Transfinite Line {9} = 21 Using Progression 1;
Line(10) = {8, 3}; Transfinite Line {10} = 501 Using Progression 1;
Line(11) = {5, 10}; Transfinite Line {11} = 281 Using Progression 1;
Line(12) = {10, 11}; Transfinite Line {12} = 51 Using Progression 1;
Line(13) = {11, 6}; Transfinite Line {13} = 281 Using Progression 1;
Line(14) = {11, 12}; Transfinite Line {14} = 501 Using Progression 1;
Line(15) = {12, 7}; Transfinite Line {15} = 281 Using Progression 1;
Line(16) = {12, 13}; Transfinite Line {16} = 51 Using Progression 1;
Line(17) = {13, 9}; Transfinite Line {17} = 281 Using Progression 1;
Line(18) = {9, 7}; Transfinite Line {18} = 51 Using Progression 1;

// Surfaces
Line Loop(19) = {1, 2, 3, 4};
Plane Surface(20) = {19};
Line Loop(21) = {5, 6, 7, -2};
Plane Surface(22) = {21};
Line Loop(23) = {11, 12, 13, -6};
Plane Surface(24) = {23};
Line Loop(25) = {8, 9, 10, -7};
Plane Surface(26) = {25};
Line Loop(27) = {14, 15, -8, -13};
Plane Surface(28) = {27};
Line Loop(29) = {16, 17, 18, -15};
Plane Surface(30) = {29};

// For structured mapping
Transfinite Surface {20};
Transfinite Surface {22};
Transfinite Surface {24};
Transfinite Surface {26};
Transfinite Surface {28};
Transfinite Surface {30};
Recombine Surface {20};
Recombine Surface {22};
Recombine Surface {24};
Recombine Surface {26};
Recombine Surface {28};
Recombine Surface {30};

// Physical groups
Physical Surface(47) = {20, 22, 24, 26, 28, 30};

