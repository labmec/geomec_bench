

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
nx = 15;
ny = 3;
pr = 1;

nf = 6;
pf = 1;

// Coordenadas dos pontos

  //Domínio Omega
  Point(1) = {0, 0, 0};
  Point(2) = {L, 0, 0};
  Point(3) = {L, h, 0};
  Point(4) = {0, h, 0};

  //Fratura 1
  Point(5) = {25, 5, 0};
  Point(6) = {75, 5, 0};  

  //Fratura 2
  Point(7) = {125, 2.5, 0};
  Point(8) = {175, 7.5, 0};  

// Fronteiras

  //Domínio Omega  
  Line(1) = {1, 2};
  Line(2) = {2, 3};
  Line(3) = {3, 4};
  Line(4) = {4, 1};  

  //Fratura
  Line(5) = {5, 6};
  Line(6) = {7, 8};

  Transfinite Line{1,3} = nx Using Progression pr;
  Transfinite Line{2,4} = ny Using Progression pr;
  Transfinite Line{5} = nf Using Progression pf;
  Transfinite Line{6} = nf Using Progression pf;
  //Transfinite Line{1} = nx Using Progression pr;


// Definição da superfície 
  Line Loop(1) = {1, 2, 3, 4};
  Plane Surface(1) = {1};
  Line{5} In Surface {1};
  Line{6} In Surface {1};
  Point{5,6} In Surface {1};
  Point{7,8} In Surface {1};

 // Transfinite Surface {1};

  If(IsquadQ)
    Recombine Surface {1};
  EndIf


  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {1};
  Physical Line("right") = {2};
  Physical Line("top") = {3};
  Physical Line("left") = {4};

  Physical Line("frac") = {5};
  Physical Line("frac2") = {6};

  Physical Point("PointLeft0") = {5};
  Physical Point("PointRight0") = {6};

  Physical Point("PointLeft1") = {7};
  Physical Point("PointRight1") = {8};

  Coherence Mesh;


