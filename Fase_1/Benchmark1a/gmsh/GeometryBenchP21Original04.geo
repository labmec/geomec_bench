

IsquadQ = 0;
 
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

// Coordenadas dos pontos

  //Domínio Omega
  Point(1) = {0, 0, 0};
  Point(2) = {L, 0, 0};
  Point(3) = {L, h, 0};
  Point(4) = {0, h, 0};

  //Fraturas

  Point(5) = {59.56005859375, 0, 0}; 
  Point(6) = {57.26318359375, 10, 0};

  Point(7) = {89.01513671875, 0, 0};
  Point(8) = {87.35107421875, 10, 0};

  Point(9) = {97.42724609375,  0, 0};
  Point(10) = {90.37890625, 10, 0};

  Point(11) = {127.0439453125,  0, 0};
  Point(12) = {121.64697265625, 10, 0};

  Point(13) = {191.47509765625, 0, 0};
  Point(14) = {188.8212890625,  10, 0};

  Point(15) = {184.17529296875, 0, 0};
  Point(16) = {182.79931640625, 10, 0};

  Point(17) = {131.27392578125, 0, 0};
  Point(18) = {125.04931640625, 10, 0};


// Fronteiras

  //Domínio Omega  
  Line(1) = {1,5};
  Line(2) = {5,7};
  Line(3) = {7,9};
  Line(4) = {9,11};
  Line(5) = {11,17};  
  Line(6) = {17,15};
  Line(7) = {15,13};  
  Line(8) = {13,2};
  Line(9) = {2,3};
  Line(10) = {3,14};
  Line(11) = {14,16};
  Line(12) = {16,18};
  Line(13) = {18,12};  
  Line(14) = {12,10};
  Line(15) = {10,8};
  Line(16) = {8,6};
  Line(17) = {6,4};
  Line(18) = {4,1};

  //Fratura

  Line(19) = {5,6};
  Line(20) = {7,8};
  Line(21) = {9,10};
  Line(22) = {11,12};
  Line(23) = {13,14};
  Line(24) = {15,16};
  Line(25) = {17,18};

  Transfinite Line{1} = 8 Using Progression pr;
  Transfinite Line{2} = 5 Using Progression pr;
  Transfinite Line{3} = 3 Using Progression pr;
  Transfinite Line{4} = 6 Using Progression pr;
  Transfinite Line{5} = 3 Using Progression pr;
  Transfinite Line{6} = 7 Using Progression pr;
  Transfinite Line{7} = 3 Using Progression pr;
  Transfinite Line{8} = 3 Using Progression pr;
  Transfinite Line{9} = 3 Using Progression pr;
  Transfinite Line{10} = 3 Using Progression pr;
  Transfinite Line{11} = 3 Using Progression pr;
  Transfinite Line{12} = 7 Using Progression pr;
  Transfinite Line{13} = 3 Using Progression pr;
  Transfinite Line{14} = 6 Using Progression pr;
  Transfinite Line{15} = 3 Using Progression pr;
  Transfinite Line{16} = 5 Using Progression pr;
  Transfinite Line{17} = 8 Using Progression pr;
  Transfinite Line{18} = 3 Using Progression pr;


  Transfinite Line{19} = nf Using Progression pr;
  Transfinite Line{20} = nf Using Progression pr;
  Transfinite Line{21} = nf Using Progression pr;
  Transfinite Line{22} = nf Using Progression pr;
  Transfinite Line{23} = nf Using Progression pr;
  Transfinite Line{24} = nf Using Progression pr;
  Transfinite Line{25} = nf Using Progression pr;

// Definição da superfície 

  Line Loop(1) = {1,19,17,18};
  Line Loop(2) = {2,20,16,-19};
  Line Loop(3) = {3,21,15,-20};
  Line Loop(4) = {4,22,14,-21};
  Line Loop(5) = {5,25,13,-22};
  Line Loop(6) = {6,24,12,-25};
  Line Loop(7) = {7,23,11,-24};
  Line Loop(8) = {8,9,10,-23};

  Plane Surface(1) = {1};
  Plane Surface(2) = {2};
  Plane Surface(3) = {3};
  Plane Surface(4) = {4};
  Plane Surface(5) = {5};
  Plane Surface(6) = {6};
  Plane Surface(7) = {7};
  Plane Surface(8) = {8};

  Transfinite Surface {1};
  Transfinite Surface {2};
  Transfinite Surface {3};
  Transfinite Surface {4};
  Transfinite Surface {5};
  Transfinite Surface {6};
  Transfinite Surface {7};
  Transfinite Surface {8};

  If(IsquadQ)

  Recombine Surface {1};
  Recombine Surface {2};
  Recombine Surface {3};
  Recombine Surface {4};
  Recombine Surface {5};
  Recombine Surface {6};
  Recombine Surface {7};
  Recombine Surface {8};

  EndIf

  Physical Surface("Omega") = {1,2,3,4,5,6,7,8};
  Physical Line("bottom") = {1,2,3,4,5,6,7,8};
  Physical Line("right") = {9};
  Physical Line("top") = {9,10,11,12,13,14,15,16,17};
  Physical Line("left") = {18};
  
  Physical Line("f1") = {19};
  Physical Line("f2") = {20};
  Physical Line("f3") = {21};
  Physical Line("f4") = {22};
  Physical Line("f5") = {23};
  Physical Line("f6") = {24};
  Physical Line("f7") = {25};


  Physical Point("PointRight0") = {5};
  Physical Point("PointLeft0") = {6};

  Physical Point("PointRight1") = {7};
  Physical Point("PointLeft1") = {8};

  Physical Point("PointRight2") = {9};
  Physical Point("PointLeft2") = {10};

  Physical Point("PointRight3") = {11};
  Physical Point("PointLeft3") = {12};

  Physical Point("PointRight4") = {17};
  Physical Point("PointLeft4") = {18};

  Physical Point("PointRight5") = {15};
  Physical Point("PointLeft5") = {16};

  Physical Point("PointRight6") = {13};
  Physical Point("PointLeft6") = {14};

  Coherence Mesh;


