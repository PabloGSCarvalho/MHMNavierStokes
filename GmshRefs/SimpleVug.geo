IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Definição de parâmtros

h = 2;
L = 2;

refp = 16;
nx1 = refp;
nx2 = refp;
nx3 = refp;
nx4 = refp;

pr = 1.0;

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

  Point(5) = {1.0, 0.5, 0.0};
  Point(6) = {1.5, 1.0, 0.0};
  Point(7) = {1.0, 1.5, 0.0};
  Point(8) = {0.5, 1.0, 0.0};
  Spline(5) = {5, 6, 7, 8, 5};
  Transfinite Line{5} = 60 Using Progression pr;
  Line Loop(2)=5;

// Definição dos domínios (superfícies) 

  Plane Surface(1) = {2};
  Plane Surface(2) = {1,2};

  If(IsquadQ)

//  Recombine Surface {1,2};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Surface("Omega2") = {2};
  Physical Line("bottom") = {1};
  Physical Line("right") = {2};
  Physical Line("top") = {3};
  Physical Line("left") = {4};

  Coherence Mesh;



