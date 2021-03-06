

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

  Point(11) = {97.42724609375,  0, 0};
  Point(12) = {90.37890625, 10, 0};

  Point(13) = {102.81201171875, 10, 0};
  Point(14) = {102.45849609375, 2.873046875, 0};

  Point(15) = {127.0439453125,  0, 0};
  Point(16) = {121.64697265625, 10, 0};

  Point(17) = {191.47509765625, 0, 0};
  Point(18) = {188.8212890625,  10, 0};


// Fronteiras

  //Domínio Omega  
  Line(1) = {1,5};
  Line(2) = {5,9};
  Line(3) = {9,11};
  Line(4) = {11,15};
  Line(5) = {15,17};  
  Line(6) = {17,2};
  Line(7) = {2,3};
  Line(8) = {3,18};
  Line(9) = {18,16};
  Line(10) = {16,13};
  Line(11) = {13,12};
  Line(12) = {12,10};
  Line(13) = {10,8};
  Line(14) = {8,6};
  Line(15) = {6,4};
  Line(16) = {4,1};

  //Fratura

  Line(17) = {5,6};
  Line(18) = {7,8};
  Line(19) = {9,10};
  Line(20) = {11,12};
  Line(21) = {13,14};
  Line(22) = {15,16};
  Line(23) = {17,18};


  Transfinite Line{1} = 8 Using Progression pr;
  Transfinite Line{2} = 6 Using Progression pr;
  Transfinite Line{3} = 5 Using Progression pr;
  Transfinite Line{4} = 9 Using Progression pr;
  Transfinite Line{5} = 8 Using Progression pr;
  Transfinite Line{6} = 6 Using Progression pr;
  Transfinite Line{7} = 5 Using Progression pr;
  Transfinite Line{8} = 6 Using Progression pr;
  Transfinite Line{9} = 8 Using Progression pr;
  Transfinite Line{10} = 9 Using Progression pr;
  Transfinite Line{11} = 3 Using Progression pr;
  Transfinite Line{12} = 5 Using Progression pr;
  Transfinite Line{13} = 4 Using Progression pr;
  Transfinite Line{14} = 5 Using Progression pr;
  Transfinite Line{15} = 8 Using Progression pr;
  Transfinite Line{16} = 5 Using Progression pr;

  Transfinite Line{17} = 5 Using Progression pr;
  Transfinite Line{18} = 3 Using Progression pr;
  Transfinite Line{19} = 5 Using Progression pr;
  Transfinite Line{20} = 5 Using Progression pr;
  Transfinite Line{21} = 4 Using Progression pr;
  Transfinite Line{22} = 5 Using Progression pr;
  Transfinite Line{23} = 5 Using Progression pr;


// Definição da superfície 

  Line Loop(1) = {1,17,15,16};
  Line Loop(2) = {2,19,13,14,-17};
  Line Loop(3) = {3,20,12,-19};
  Line Loop(4) = {4,22,10,11,-20};
  Line Loop(5) = {5,23,9,-22};
  Line Loop(6) = {6,7,8,-23};


  Plane Surface(1) = {1};
  Plane Surface(2) = {2};
  Plane Surface(3) = {3};
  Plane Surface(4) = {4};
  Plane Surface(5) = {5};
  Plane Surface(6) = {6};

  Line{18} In Surface{2};
  Line{21} In Surface{4};

  Point{7,8} In Surface {2};
  Point{13,14} In Surface {4};

  Transfinite Surface {1};
  Transfinite Surface {3};
  Transfinite Surface {5};
  Transfinite Surface {6};

  If(IsquadQ)

  Recombine Surface {1};
  Recombine Surface {2};
  Recombine Surface {3};
  Recombine Surface {4};
  Recombine Surface {5};
  Recombine Surface {6};

  EndIf

  Physical Surface("Omega") = {1,2,3,4,5,6};
  Physical Line("bottom") = {1,2,3,4,5,6,7,8,9,10,11,12};
  Physical Line("top") = {8,9,10,11,12,13,14,15};
  Physical Line("right") = {7};
  Physical Line("left") = {16};
  
  Physical Line("f1") = {17}; 
  Physical Line("f2") = {18};
  Physical Line("f3") = {19};
  Physical Line("f4") = {20};
  Physical Line("f5") = {21};
  Physical Line("f6") = {22};
  Physical Line("f7") = {23};
      

  Physical Point("PointRight0") = {5};
  Physical Point("PointLeft0") = {6};

  Physical Point("PointRight1") = {7};
  Physical Point("PointLeft1") = {8};

  Physical Point("PointRight2") = {9};
  Physical Point("PointLeft2") = {10};

  Physical Point("PointRight3") = {11};
  Physical Point("PointLeft3") = {12};

  Physical Point("PointRight4") = {13};
  Physical Point("PointLeft4") = {14};

  Physical Point("PointRight5") = {15};
  Physical Point("PointLeft5") = {16};

  Physical Point("PointRight6") = {17};
  Physical Point("PointLeft6") = {18};


  Coherence Mesh;

