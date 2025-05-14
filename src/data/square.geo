// square.geo
// Unit square [0,1]×[0,1] with characteristic length 0.1
h=1;
// Define corner points
Point(1) = {0.0, 0.0, 0.0, h};  // bottom‐left
Point(2) = {1.0, 0.0, 0.0, h};  // bottom‐right
Point(3) = {1.0, 1.0, 0.0, h};  // top‐right
Point(4) = {0.0, 1.0, 0.0, h};  // top‐left

// Define edges
Line(1) = {1, 2};  // bottom
Line(2) = {2, 3};  // right
Line(3) = {3, 4};  // top
Line(4) = {4, 1};  // left

// Create loop and surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Physical groups for boundary conditions
Physical Surface("Square") = {1};
Physical Line("Dirichlet") = {1,2,3,4};
