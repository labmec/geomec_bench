

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
nx = 5;
ny = 3;
pr = 1;

nf = 3;
pf = 1;

// Coordenadas dos pontos

  //Domínio Omega
  Point(1) = {0, 0, 0};
  Point(2) = {L, 0, 0};
  Point(3) = {L, h, 0};
  Point(4) = {0, h, 0};

  //Fratura 1
  Point(5) = {66, 0, 0};
  Point(6) = {66, 10, 0};  

  //Fratura 2
  Point(7) = {137, 0, 0};
  Point(8) = {137, 10, 0};  


// Fronteiras

  //Domínio Omega  
  Line(1) = {1, 5};
  Line(2) = {5, 7};
  Line(3) = {7, 2};
  Line(4) = {2, 3};
  Line(5) = {3, 8};
  Line(6) = {8, 6};
  Line(7) = {6, 4};
  Line(8) = {4, 1};

  //Fratura
  Line(9) = {5, 6};
  Line(10) = {7, 8};


  Transfinite Line{1,2,3,5,6,7} = nx Using Progression pr;
  Transfinite Line{4,8} = ny Using Progression pr;
  Transfinite Line{9,10} = nf Using Progression pr;

// Definição da superfície 
  Line Loop(1) = {1, 9, 7, 8};
  Line Loop(2) = {2, 10, 6, -9};
  Line Loop(3) = {3, 4, 5, -10};
  Plane Surface(1) = {1};
  Plane Surface(2) = {2};
  Plane Surface(3) = {3};

//  Transfinite Surface {1};
//  Transfinite Surface {2};
//  Transfinite Surface {3};

  If(IsquadQ)
    Recombine Surface {1};
    Recombine Surface {2};
    Recombine Surface {3};
  EndIf


  Physical Surface("Omega") = {1,2,3};
  Physical Line("bottom") = {1,2,3};
  Physical Line("right") = {4};
  Physical Line("top") = {5,6,7};
  Physical Line("left") = {8};

  Physical Line("f1") = {9};
  Physical Line("f2") = {10};

  Physical Point("PointLeft0") = {5};
  Physical Point("PointRight0") = {6};

  Physical Point("PointLeft1") = {7};
  Physical Point("PointRight1") = {8};

  Coherence Mesh;


