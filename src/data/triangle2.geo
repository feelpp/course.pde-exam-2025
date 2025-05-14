// triangle2.geo
// Single right‐angle triangle [0,0]–[1,0]–[0,1] with characteristic length 0.1

// Define corner points
Point(1) = {0.0, 0.0, 0.0, 0.1};
Point(2) = {1.0, 0.0, 0.0, 0.1};
Point(3) = {0.0, 1.0, 0.0, 0.1};

// Define edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};

// Create loop and surface
Line Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};

// Physical groups
Physical Surface("Triangle") = {1};
Physical Line("Edge12")    = {1};
Physical Line("Edge23")    = {2};
Physical Line("Edge31")    = {3};