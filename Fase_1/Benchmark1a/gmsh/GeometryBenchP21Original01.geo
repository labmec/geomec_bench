

IsquadQ = 1;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

// Definição de parâmtros

a = 5; 
b = 5;
h = 10;
L = 200;
Lf = 200;

n_bc = 2;
nx = 8;
nx1 = 4;
nx2 = 3;
nx3 = 4;
nx4 = 5;

ny = 3;
pr = 1;

nf = 4;

// Coordenadas dos pontos

  //Domínio Omega
  Point(1) = {0, 0, 0};
  Point(2) = {L, 0, 0};
  Point(3) = {L, h, 0};
  Point(4) = {0, h, 0};

  //Fraturas

  Point(5) = {59.56005859375, 0, 0}; 
  Point(6) = {57.26318359375, 10, 0};

  Point(7) = {76.92333984375, 5.39453125, 0};
  Point(8) = {76.5068359375, 10, 0};

  Point(9) = {89.01513671875, 0, 0};
  Point(10) = {87.35107421875, 10, 0};


// Fronteiras

  //Domínio Omega  
  Line(1) = {1,5};
  Line(2) = {5,9};
  Line(3) = {9,2};
  Line(4) = {2,3};
  Line(5) = {3,10};
  Line(6) = {10,8};
  Line(7) = {8,6};
  Line(8) = {6,4};
  Line(9) = {4,1};

  //Fratura

  Line(10) = {5,6};
  Line(11) = {7,8};
  Line(12) = {9,10};

  Transfinite Line{1,8} = nx1 Using Progression pr;
  Transfinite Line{6,7} = nx2 Using Progression pr;
  Transfinite Line{2} = nx3 Using Progression pr;
  Transfinite Line{3,5} = nx4 Using Progression pr;
  Transfinite Line{4,9} = ny Using Progression pr;

  Transfinite Line{10} = nf Using Progression pr;
  Transfinite Line{11} = nf Using Progression pr;
  Transfinite Line{12} = nf Using Progression pr;    


// Definição da superfície 

  Line Loop(1) = {1,2,3,4,5,6,7,8,9};
  Plane Surface(1) = {1};

  Line{10} In Surface{1};
  Line{11} In Surface{1};
  Line{12} In Surface{1};

  Point{5,6} In Surface {1};
  Point{7,8} In Surface {1};
  Point{9,10} In Surface {1};

  //Transfinite Surface {1};

  If(IsquadQ)

  Recombine Surface {1};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {1,2,3};
  Physical Line("top") = {5,6,7,8};
  Physical Line("right") = {4};
  Physical Line("left") = {9};
  
  Physical Line("f1") = {10};
  Physical Line("f2") = {11};
  Physical Line("f3") = {12};
  

  Physical Point("PointRight0") = {5};
  Physical Point("PointLeft0") = {6};

  Physical Point("PointRight1") = {7};
  Physical Point("PointLeft1") = {8};

  Physical Point("PointRight2") = {9};
  Physical Point("PointLeft2") = {10};



  Coherence Mesh;


