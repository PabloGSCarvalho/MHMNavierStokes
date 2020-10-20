IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Definição de parâmtros

h = 2;
L = 3;

refp = 10;
refc1 = 20;
refc5 = 20;
nx1 = refp;
nx2 = refp;
nx3 = refp;
nx4 = refp;
lc = .1;
pr = 1;

// Coordenadas dos pontos

  //Domínio Omega - Meio poroso
  Point(1) = {0, 0, 0, lc};
  Point(2) = {L, 0, 0, lc};
  Point(3) = {L, h, 0, lc};
  Point(4) = {0, h, 0, lc};

// Fronteiras

  //Domínio Omega  
  Line(1) = {1,2};
  Line(2) = {2,3};
  Line(3) = {3,4};
  Line(4) = {4,1};

//  Transfinite Line{1,2} = nx1 Using Progression pr;
//  Transfinite Line{2,3} = nx2 Using Progression pr;
//  Transfinite Line{3,4} = nx3 Using Progression pr;
//  Transfinite Line{4,1} = nx4 Using Progression pr;
  
  Line Loop(1) = {1,2,3,4};

// Vug 1
 co01 = 0.5; co02 = 0.5;  //center point
 h01 = 0.2; h02 = 0.5;    //dimensions
 Point(5) = {co01,co02,0,lc}; 
 Point(6) = {co01-h01/2.,co02,0,lc};
 Point(7) = {co01+h01/2,co02,0,lc};
 Point(8) = {co01,co02+h02/2.,0,lc};
 Point(9) = {co01,co02-h02/2.,0,lc};
 
 Ellipse(5) = {8,5,6};
 Ellipse(6) = {7,5,8};
 Ellipse(7) = {6,5,9};
 Ellipse(8) = {9,5,7};

//Transfinite Line{5,6,7,8} = refc1 Using Progression pr;
Curve Loop(2) = {5,7,8,6};

// Vug 2
 co01 = 1.5; co02 = 1.65;  //center point
 h01 = 0.6; h02 = 0.3;    //dimensions
 Point(10) = {co01,co02,0,lc}; 
 Point(11) = {co01-h01/2.,co02,0,lc};
 Point(12) = {co01+h01/2,co02,0,lc};
 Point(13) = {co01,co02+h02/2.,0,lc};
 Point(14) = {co01,co02-h02/2.,0,lc};
 
 Ellipse(9) = {13,10,11};
 Ellipse(10) = {12,10,13};
 Ellipse(11) = {11,10,14};
 Ellipse(12) = {14,10,12};

//Transfinite Line{9,10,11,12} = refc1 Using Progression pr;
Curve Loop(3) = {9,11,12,10};

// Vug 3

 co01 = 1.25; co02 = 0.9;  //center point
 h01 = 0.5; h02 = 0.5;    //dimensions
 Point(15) = {co01,co02,0,lc}; 
 Point(16) = {co01-h01/2.,co02,0,lc};
 Point(17) = {co01+h01/2,co02,0,lc};
 Point(18) = {co01,co02+h02/2.,0,lc};
 Point(19) = {co01,co02-h02/2.,0,lc};
 
 Ellipse(13) = {18,15,16};
 Ellipse(14) = {17,15,18};
 Ellipse(15) = {16,15,19};
 Ellipse(16) = {19,15,17};

//Transfinite Line{13,14,15,16} = refc1 Using Progression pr;
Curve Loop(4) = {13,15,16,14};

// Vug 4

 co01 = 0.6; co02 = 1.3;  //center point
 h01 = 0.5; h02 = 0.5;    //dimensions
// Point(20) = {co01,co02,0,lc}; 
 Point(20) = {co01-h01/4.,co02,0,lc};
 Point(21) = {co01+h01/4,co02,0,lc};
 Point(22) = {co01,co02+h02/2.,0,lc};
 Point(23) = {co01,co02-h02/4,0,lc};
 
 Spline(17) = {20,22,21,23,20};

//Transfinite Line{17} = refc5 Using Progression pr;
Curve Loop(5) = {17};

// Vug 5

 co01 = 2.45; co02 = 1.5;  //center point
 h01 = 0.3; h02 = 0.6; rot = Pi/4;   //dimensions

 Point(24) = {co01,co02,0,lc}; 
 Point(25) = {co01+(-h01/2.)*Cos(rot)+(0)*Sin(rot),co02-(-h01/2.)*Sin(rot)+(0)*Cos(rot),0,lc};
 Point(26) = {co01+(h01/2.)*Cos(rot)+(0)*Sin(rot),co02-(h01/2.)*Sin(rot)+(0)*Cos(rot),0,lc};
 Point(27) = {co01+(0)*Cos(rot)+(h02/2.)*Sin(rot),co02-(0)*Sin(rot)+(h02/2.)*Cos(rot),0,lc};
 Point(28) = {co01+(0)*Cos(rot)+(-h02/2.)*Sin(rot),co02-(0)*Sin(rot)+(-h02/2.)*Cos(rot),0,lc};
 
 Ellipse(18) = {27,24,25};
 Ellipse(19) = {26,24,27};
 Ellipse(20) = {25,24,28};
 Ellipse(21) = {28,24,26};

 Curve Loop(6) = {18,20,21,19};

// Vug 6

 co01 = 1.9; co02 = 0.9;  //center point
 h01 = 0.2; h02 = 0.2;    //dimensions
 Point(29) = {co01-h01/2.,co02,0,lc};
 Point(30) = {co01+h01/2,co02,0,lc};
 Point(31) = {co01,co02+h02/1.2,0,lc};
 Point(32) = {co01,co02-h02/1.2,0,lc};
 
 Spline(22) = {32,30,31,29,32};
 Curve Loop(7) = {22};

// Vug 7

 co01 = 2.6; co02 = 0.75;  //center point
 h01 = 0.3; h02 = 0.3;    //dimensions
 Point(33) = {co01-h01/2.,co02,0,lc};
 Point(34) = {co01+h01/2,co02,0,lc};
 Point(35) = {co01,co02+h02/1.5,0,lc};
 Point(36) = {co01,co02-h02/1.1,0,lc};
 
 Spline(23) = {36,34,35,33,36};
 Curve Loop(8) = {23};

// Vug 8

 Point(37) = {2.12-0.4,0.5,0,lc};
 Point(38) = {2.5-0.4,0.35,0,lc};
 Point(39) = {2.8-0.4,0.3,0,lc};
 Point(40) = {2.7-0.4,0.14,0,lc};
 Point(41) = {2.2-0.4,0.23,0,lc}; 
 Point(42) = {1.9-0.4,0.13,0,lc};
 Point(43) = {1.9-0.4,0.32,0,lc};

 Spline(24) = {37,38,39,40,41,42,43,37};
 Curve Loop(9) = {24};


// Definição dos domínios (superfícies) 

  Plane Surface(2) = {1,2,3,4,5,6,7,8,9};
  Plane Surface(1) = {2};
  Plane Surface(3) = {3};
  Plane Surface(4) = {4};
  Plane Surface(5) = {5};
  Plane Surface(6) = {6};
  Plane Surface(7) = {7};
  Plane Surface(8) = {8};
  Plane Surface(9) = {9};

Field[1] = Distance;
Field[1].NNodesByEdge = 100;
Field[1].EdgesList = {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc/7;
Field[2].LcMax = lc;
Field[2].DistMin = 0.05;
Field[2].DistMax = 0.3;

Field[7] = Min;
Field[7].FieldsList = {2};
Background Field = 7;


  If(IsquadQ)

  Recombine Surface {1,2,3,4,5,6,7,8,9};

  EndIf

  Physical Surface("Omega") = {1,3,4,5,6,7,8,9};
  Physical Surface("Omega2") = {2};
  Physical Line("bottom") = {1};
  Physical Line("right") = {2};
  Physical Line("top") = {3};
  Physical Line("left") = {4};
  
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;

  Coherence Mesh;


