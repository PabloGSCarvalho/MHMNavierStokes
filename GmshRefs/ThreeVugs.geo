IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Definição de parâmtros

h = 2;
L = 2;

refp = 16;
refc1 = 48;
refc5 = 60;
nx1 = refp;
nx2 = refp;
nx3 = refp;
nx4 = refp;

pr = 1;

// Coordenadas dos pontos

  //Domínio Omega - Meio poroso
  Point(1) = {0, 0, 0};
  Point(2) = {L, 0, 0};
  Point(3) = {L, h, 0};
  Point(4) = {0, h, 0};

// Fronteiras

  //Domínio Omega  
  Line(1) = {1,2};
  Line(2) = {2,3};
  Line(3) = {3,4};
  Line(4) = {4,1};

  Transfinite Line{1,2} = nx1 Using Progression pr;
  Transfinite Line{2,3} = nx2 Using Progression pr;
  Transfinite Line{3,4} = nx3 Using Progression pr;
  Transfinite Line{4,1} = nx4 Using Progression pr;
  
  Line Loop(1) = {1,2,3,4};

// Vug 1
 cl = 0.0; 
 co01 = 0.5; co02 = 0.5;  //center point
 h01 = 0.2; h02 = 0.5;    //dimensions
 Point(5) = {co01,co02,0,cl}; 
 Point(6) = {co01-h01/2.,co02,0,cl};
 Point(7) = {co01+h01/2,co02,0,cl};
 Point(8) = {co01,co02+h02/2.,0,cl};
 Point(9) = {co01,co02-h02/2.,0,cl};
 
 Ellipse(5) = {8,5,6};
 Ellipse(6) = {7,5,8};
 Ellipse(7) = {6,5,9};
 Ellipse(8) = {9,5,7};

Transfinite Line{5,6,7,8} = refc1 Using Progression pr;
Curve Loop(2) = {5,7,8,6};

// Vug 2
 cl = 0.0; 
 co01 = 1.5; co02 = 1.65;  //center point
 h01 = 0.6; h02 = 0.3;    //dimensions
 Point(10) = {co01,co02,0,cl}; 
 Point(11) = {co01-h01/2.,co02,0,cl};
 Point(12) = {co01+h01/2,co02,0,cl};
 Point(13) = {co01,co02+h02/2.,0,cl};
 Point(14) = {co01,co02-h02/2.,0,cl};
 
 Ellipse(9) = {13,10,11};
 Ellipse(10) = {12,10,13};
 Ellipse(11) = {11,10,14};
 Ellipse(12) = {14,10,12};

Transfinite Line{9,10,11,12} = refc1 Using Progression pr;
Curve Loop(3) = {9,11,12,10};

// Vug 3
 cl = 0.0; 
 co01 = 1.25; co02 = 0.9;  //center point
 h01 = 0.5; h02 = 0.5;    //dimensions
 Point(15) = {co01,co02,0,cl}; 
 Point(16) = {co01-h01/2.,co02,0,cl};
 Point(17) = {co01+h01/2,co02,0,cl};
 Point(18) = {co01,co02+h02/2.,0,cl};
 Point(19) = {co01,co02-h02/2.,0,cl};
 
 Ellipse(13) = {18,15,16};
 Ellipse(14) = {17,15,18};
 Ellipse(15) = {16,15,19};
 Ellipse(16) = {19,15,17};

Transfinite Line{13,14,15,16} = refc1 Using Progression pr;
Curve Loop(4) = {13,15,16,14};

// Vug 4
 cl = 0.0; 
 co01 = 0.6; co02 = 1.3;  //center point
 h01 = 0.5; h02 = 0.5;    //dimensions
// Point(20) = {co01,co02,0,cl}; 
 Point(20) = {co01-h01/4.,co02,0,cl};
 Point(21) = {co01+h01/4,co02,0,cl};
 Point(22) = {co01,co02+h02/2.,0,cl};
 Point(23) = {co01,co02-h02/4,0,cl};
 
 Spline(17) = {20,22,21,23,20};

Transfinite Line{17} = refc5 Using Progression pr;
Curve Loop(5) = {17};

// Definição dos domínios (superfícies) 

  Plane Surface(2) = {1,2,3,4,5};
  Plane Surface(1) = {2};
  Plane Surface(3) = {3};
  Plane Surface(4) = {4};
  Plane Surface(5) = {5};

  If(IsquadQ)

  Recombine Surface {1,2};

  EndIf

  Physical Surface("Omega") = {1,3,4,5};
  Physical Surface("Omega2") = {2};
  Physical Line("bottom") = {1};
  Physical Line("right") = {2};
  Physical Line("top") = {3};
  Physical Line("left") = {4};
  
  Coherence Mesh;


