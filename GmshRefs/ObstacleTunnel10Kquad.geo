IsquadQ = 1;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Definição de parâmtros

h = 4;
L = 12;

refp = 32;
refc1 = 24*2;
nx1 = refp*2*1.2;
nx2 = refp*2;
nx3 = refp;
nx4 = refp;
lc = 0.09;
pr = 1.0;
pr2 = 1.;

// Coordenadas dos pontos

  //Domínio Omega - Meio poroso
  Point(1) = {0, 0, 0,lc};
  Point(2) = {L, 0, 0,lc};
  Point(3) = {L, h, 0,lc};
  Point(4) = {0, h, 0,lc};

// Fronteiras

  //Domínio Omega  
  Line(1) = {1,2};
  Line(2) = {2,3};
  Line(3) = {3,4};
  Line(4) = {4,1};

//  Transfinite Line{1} = nx1 Using Progression pr2;
//  Transfinite Line{3} = nx1 Using Progression pr2;
  Transfinite Line{2,4} = nx2 Using Progression pr;
//  Transfinite Line{3,4} = nx3 Using Progression pr;
//  Transfinite Line{4,1} = nx4 Using Progression pr;
  
  Line Loop(1) = {1,2,3,4};

// Obstacle
 cl = 0.5; 
 co01 = 2.0; co02 = 2.0;  //center point
 h01 = 1.0; h02 = 1.0;    //dimensions
 Point(5) = {co01,co02,0,cl}; 
 Point(6) = {co01-h01/2.,co02,0,cl};
 Point(7) = {co01+h01/2,co02,0,cl};
 Point(8) = {co01,co02+h02/2.,0,cl};
 Point(9) = {co01,co02-h02/2.,0,cl};
 
 Ellipse(5) = {8,5,6};
 Ellipse(6) = {7,5,8};
 Ellipse(7) = {6,5,9};
 Ellipse(8) = {9,5,7};

//Transfinite Line{5,6,7,8} = refc1 Using Progression pr;
Curve Loop(2) = {5,6,7,8};

// Definição dos domínios (superfícies) 

  Plane Surface(1) = {1,2};

Field[1] = Distance;
Field[1].EdgesList = {1,3};
Field[1].NNodesByEdge = 100;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc/1.5;
Field[2].LcMax = lc;
Field[2].DistMin = 0.0;
Field[2].DistMax = 0.4;

Field[3] = Distance;
Field[3].EdgesList = {5,6,7,8};
Field[3].NNodesByEdge = 100;

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = lc/4;
Field[4].LcMax = lc;
Field[4].DistMin = 0.1;
Field[4].DistMax = 1.5;


Field[7] = Min;
Field[7].FieldsList = {2,4};
Background Field = 7;



  If(IsquadQ)

  Recombine Surface {1,2};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {1};
  Physical Line("right") = {2};
  Physical Line("top") = {3};
  Physical Line("left") = {4};
  Physical Line("obstacle") = {5,6,7,8};

  Coherence Mesh;


Mesh.CharacteristicLengthExtendFromBoundary = 0;


