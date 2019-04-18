dx=0.2;

Point(1)={0.0,0.0,0.0,dx};
Point(2)={1.0,0.0,0.0,dx};
Point(3)={1.0,1.0,0.0,dx};
Point(4)={0.0,1.0,0.0,dx};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Line("left") = {4};
//+
Physical Line("right") = {2};
//+
Physical Line("bottom") = {1};
//+
Physical Line("top") = {3};
//+
Physical Surface("matrix") = {1};
