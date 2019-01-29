

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

nf = 5;
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
  Point(7) = {125, 0.1, 0};
  Point(8) = {175, 9.9, 0};  

  //Fratura 2
  Point(9) = {180, 0.1, 0};
  Point(10) = {188, 10, 0};  


// Fronteiras

  //Domínio Omega  
  Line(1) = {1, 2};
  Line(2) = {2, 3};
  Line(3) = {3, 10};
  Line(4) = {10, 4};
  Line(5) = {4, 1};  

  //Fratura
  Line(6) = {5, 6};
  Line(7) = {7, 8};
  Line(8) = {9, 10};

  Transfinite Line{1,4} = nx Using Progression pr;
  Transfinite Line{2,5} = ny Using Progression pr;
  Transfinite Line{3} = ny Using Progression pf;
  Transfinite Line{6} = nf Using Progression pf;
  Transfinite Line{7} = nf Using Progression pf;
  Transfinite Line{8} = ny Using Progression pf;
  //Transfinite Line{1} = nx Using Progression pr;


// Definição da superfície 
  Line Loop(1) = {1, 2, 3, 4, 5};
  Plane Surface(1) = {1};
  Line{6} In Surface {1};
  Line{7} In Surface {1};
  Line{8} In Surface {1};
  Point{5,6} In Surface {1};
  Point{7,8} In Surface {1};
  Point{9,10} In Surface {1};

 // Transfinite Surface {1};

  If(IsquadQ)
    Recombine Surface {1};
  EndIf


  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {1};
  Physical Line("right") = {2};
  Physical Line("top") = {3,4};
  Physical Line("left") = {5};

  Physical Line("f1") = {6};
  Physical Line("f2") = {7};
  Physical Line("f3") = {8};

  Physical Point("PointLeft0") = {5};
  Physical Point("PointRight0") = {6};

  Physical Point("PointLeft1") = {7};
  Physical Point("PointRight1") = {8};

  Physical Point("PointLeft2") = {9};
  Physical Point("PointRight2") = {10};

  Coherence Mesh;


