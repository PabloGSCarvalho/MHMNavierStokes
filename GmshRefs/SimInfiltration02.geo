IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Definição de parâmtros

h = 7;
hp = 3;
hp1 = 5;
L = 64;
L1 = 20;
L2 = 24;


refp = 10;
refc1 = 20;
refc5 = 20;
nfl = 8;
npr = 6;
nx = 32;
lc = .6;
pr = 1;

// Coordenadas dos pontos

  //Domínio Omega - Meio poroso
  Point(1) = {0, 0, 0, lc};
  Point(2) = {L, 0, 0, lc};
  Point(3) = {L, h, 0, lc};
  Point(4) = {0, h, 0, lc};

  Point(5) = {0, hp, 0, lc};
  Point(6) = {L, hp, 0, lc};

  Point(7) = {L1, hp, 0, lc};
  Point(8) = {L1, hp1, 0, lc};
  Point(9) = {L2, hp1, 0, lc};
  Point(10) = {L2, hp, 0, lc};


// Fronteiras

  //Domínio Omega  
  Line(1) = {1,2};
  Line(2) = {2,6};
  Line(3) = {6,3};
  Line(4) = {3,4};
  Line(5) = {4,5};
  Line(6) = {5,1};

  Line(7) = {5,7};
  Line(8) = {7,8};
  Line(9) = {8,9};
  Line(10) = {9,10};
  Line(11) = {10,6};

//  Transfinite Line{3,5} = nfl Using Progression pr;
//  Transfinite Line{2,6} = npr Using Progression pr;
//  Transfinite Line{1,4,7,8,9,10,11} = nx Using Progression pr;
 
  Line Loop(1) = {1,2,-11,-10,-9,-8,-7,6};
  Line Loop(2) = {7,8,9,10,11,3,4,5};

// Definição dos domínios (superfícies) 

  Plane Surface(2) = {1};
  Plane Surface(1) = {2};

Field[1] = Distance;
Field[1].NNodesByEdge = 100;
Field[1].EdgesList = {8,9,10};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc/5;
Field[2].LcMax = lc;
Field[2].DistMin = 0.1;
Field[2].DistMax = 2.5;

Field[7] = Min;
Field[7].FieldsList = {2};
Background Field = 7;


 // Transfinite Surface {1,2};
  If(IsquadQ)

  Recombine Surface {1,2};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Surface("Omega2") = {2};
  Physical Line("bottom") = {1};
  Physical Line("right") = {3};
  Physical Line("right2") = {2};
  Physical Line("top") = {4};
  Physical Line("left") = {5};
  Physical Line("left2") = {6};
  
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;

  Coherence Mesh;


