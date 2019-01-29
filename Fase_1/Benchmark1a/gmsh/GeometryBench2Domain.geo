

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
ny = 5;
pr = 1;

nf = 5;
pf = 1;

// Coordenadas dos pontos

  //Domínio Omega
  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {L, 0, 0, 1e+22};
  Point(3) = {L, h, 0, 1e+22};
  Point(4) = {0, h, 0, 1e+22};

  //Fratura 1
  Point(5) = {100, 0, 0, 1e+22};
  Point(6) = {100, 10, 0, 1e+22};  


// Fronteiras

  //Domínio Omega  
  Line(1) = {1, 5};
  Line(2) = {5, 2};
  Line(3) = {2, 3};
  Line(4) = {3, 6};
  Line(5) = {6, 4};
  Line(6) = {4, 1};

  //Fratura
  Line(7) = {5, 6};

  Transfinite Line{1,2,4,5} = nx Using Progression pr;
  Transfinite Line{3,6} = ny Using Progression pr;
  Transfinite Line{7} = nf Using Progression pr;

// Definição da superfície 
  Line Loop(1) = {1, 7, 5, 6};
  Line Loop(2) = {2, 3, 4, -7};
  Plane Surface(1) = {1};
  Plane Surface(2) = {2};

  Transfinite Surface {1};
  Transfinite Surface {2};

  If(IsquadQ)
    Recombine Surface {1};
    Recombine Surface {2};
  EndIf


  Physical Surface("Omega") = {1,2};
  Physical Line("bottom") = {1,2};
  Physical Line("top") = {4,5};
  Physical Line("right") = {3};
  Physical Line("left") = {6};

  Physical Line("f1") = {7};

  Physical Point("PointLeft0") = {5};
  Physical Point("PointRight0") = {6};


  Coherence Mesh;


