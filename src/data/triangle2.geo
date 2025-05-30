// triangle2.geo
// Single right‐angle triangle [0,0]–[1,0]–[0,1] with characteristic length 0.1

h=0.1;
// Define corner points
Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {1.0, 0.0, 0.0, h};
Point(3) = {0.0, 1.0, 0.0, h};

// Define edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};

// Create loop and surface
Line Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};

// Physical groups
Physical Surface("Triangle") = {1};
Physical Line("Dirichlet")    = {1,2,3};

