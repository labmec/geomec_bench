

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

//  Point(19) = {157.5068359375,  1.970703125, 0};
//  Point(19) = {157.3876953125,  1, 0};

  Point(20) = {157.3876953125,  0, 0};



  Point(21) = {177.08203125,  0, 0};
  Point(22) = {173.21044921875, 7.1572265625, 0};

  Point(23) = {184.17529296875, 0, 0};
  Point(24) = {182.79931640625, 10, 0};

  Point(25) = {131.27392578125, 0, 0};
  Point(26) = {125.04931640625, 10, 0};

  Point(27) = {192.71142578125, 0, 0};
  Point(28) = {192.60400390625, 1.470703125, 0};

  Point(29) = {193.7646484375,  3.07763671875, 0};
  Point(30) = {192.798828125, 10, 0};

  Point(31) = {195.78076171875, 0, 0};
  Point(32) = {195.13818359375, 1.36083984375, 0};

// Fronteiras

  //Domínio Omega  
  Line(1) = {1,5};
  Line(2) = {5,9};
  Line(3) = {9,11};
  Line(4) = {11,15};
  Line(5) = {15,25};  
  Line(6) = {25,20};
  Line(7) = {20,21};
  Line(8) = {21,23};
  Line(9) = {23,17};
  Line(10) = {17,27};
  Line(11) = {27,31};
  Line(12) = {31,2};
  Line(13) = {2,3};
  Line(14) = {3,30};
  Line(15) = {30,18};
  Line(16) = {18,24};
  Line(17) = {24,26};
  Line(18) = {26,16};
  Line(19) = {16,13};
  Line(20) = {13,12};
  Line(21) = {12,10};
  Line(22) = {10,8};
  Line(23) = {8,6};
  Line(24) = {6,4};
  Line(25) = {4,1};

  //Fratura

  Line(26) = {5,6};
  Line(27) = {7,8};
  Line(28) = {9,10};
  Line(29) = {11,12};
  Line(30) = {13,14};
  Line(31) = {15,16};
  Line(32) = {17,18};
//  Line(33) = {19,20};
  Line(34) = {21,22};
  Line(35) = {23,24};
  Line(36) = {25,26};
  Line(37) = {27,28};
  Line(38) = {29,30};
  Line(39) = {31,32};

  Transfinite Line{1} = 5 Using Progression pr;
  Transfinite Line{2} = 6 Using Progression pr;
  Transfinite Line{3} = 4 Using Progression pr;
  Transfinite Line{4} = 11 Using Progression pr;
  Transfinite Line{5} = 4 Using Progression pr;
  Transfinite Line{6} = 8 Using Progression pr;
  Transfinite Line{7} = 8 Using Progression pr;
  Transfinite Line{8} = 4 Using Progression pr;
  Transfinite Line{9} = 5 Using Progression pr;
  Transfinite Line{10} = 5 Using Progression pr;
  Transfinite Line{11} = 4 Using Progression pr;
  Transfinite Line{12} = 4 Using Progression pr;
  Transfinite Line{13} = 4 Using Progression pr;
  Transfinite Line{14} = 4 Using Progression pr;
  Transfinite Line{15} = 4 Using Progression pr;
  Transfinite Line{16} = 5 Using Progression pr;
  Transfinite Line{17} = 15 Using Progression pr;
  Transfinite Line{18} = 4 Using Progression pr;
  Transfinite Line{19} = 7 Using Progression pr;
  Transfinite Line{20} = 5 Using Progression pr;
  Transfinite Line{21} = 4 Using Progression pr;
  Transfinite Line{22} = 5 Using Progression pr;
  Transfinite Line{23} = 5 Using Progression pr;
  Transfinite Line{24} = 5 Using Progression pr;
  Transfinite Line{25} = 4 Using Progression pr;


  Transfinite Line{26} = 4 Using Progression pr;
  Transfinite Line{27} = 4 Using Progression pr;
  Transfinite Line{28} = nf Using Progression pr;    
  Transfinite Line{29} = nf Using Progression pr;
  Transfinite Line{30} = nf Using Progression pr;
  Transfinite Line{31} = 5 Using Progression pr;
  Transfinite Line{32} = 5 Using Progression pr;
//  Transfinite Line{33} = 3 Using Progression pr;  
  Transfinite Line{34} = nf Using Progression pr;    
  Transfinite Line{35} = 5 Using Progression pr;
  Transfinite Line{36} = 5 Using Progression pr;
  Transfinite Line{37} = 3 Using Progression pr;
  Transfinite Line{38} = nf Using Progression pr;
  Transfinite Line{39} = 3 Using Progression pr;


// Definição da superfície 

  Line Loop(1) = {1,26,24,25};
  Line Loop(2) = {2,28,22,23,-26};
  Line Loop(3) = {3,29,21,-28};
  Line Loop(4) = {4,31,19,20,-29};
  Line Loop(5) = {5,36,18,-31};
  Line Loop(6) = {6,7,8,35,17,-36};
  Line Loop(7) = {9,32,16,-35};
  Line Loop(8) = {10,11,12,13,14,15,-32};

  Plane Surface(1) = {1};
  Plane Surface(2) = {2};
  Plane Surface(3) = {3};
  Plane Surface(4) = {4};
  Plane Surface(5) = {5};
  Plane Surface(6) = {6};
  Plane Surface(7) = {7};
  Plane Surface(8) = {8};

  Line{27} In Surface{2};
  Line{30} In Surface{4};
//  Line{33} In Surface{6};
  Line{34} In Surface{6};
  Line{37} In Surface{8};
  Line{38} In Surface{8};
  Line{39} In Surface{8};

  Point{7,8} In Surface {2};
  Point{13,14} In Surface {4};
//  Point{19,20} In Surface {6};
  Point{21,22} In Surface {6};
  Point{27,28} In Surface {8};
  Point{29,30} In Surface {8};
  Point{31,32} In Surface {8};

  Transfinite Surface {1};
  Transfinite Surface {3};
  Transfinite Surface {5};
  Transfinite Surface {7};

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
  Physical Line("bottom") = {1,2,3,4,5,6,7,8,9,10,11,12};
  Physical Line("top") = {14,15,16,17,18,19,20,21,22,23,24};
  Physical Line("right") = {13};
  Physical Line("left") = {25};
  
  Physical Line("f1") = {26}; 
  Physical Line("f2") = {27};
  Physical Line("f3") = {28};
  Physical Line("f4") = {29};
  Physical Line("f5") = {30};
  Physical Line("f6") = {31};
  Physical Line("f7") = {32};
//  Physical Line("f8") = {33};
  Physical Line("f8") = {34};
  Physical Line("f9") = {35}; 
  Physical Line("f10") = {36};
  Physical Line("f11") = {37};
  Physical Line("f12") = {38};
  Physical Line("f13") = {39};      
  

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

//  Physical Point("PointRight7") = {19};
//  Physical Point("PointLeft7") = {20};

  Physical Point("PointRight7") = {21};
  Physical Point("PointLeft7") = {22};

  Physical Point("PointRight8") = {23};
  Physical Point("PointLeft8") = {24};

  Physical Point("PointRight9") = {25};
  Physical Point("PointLeft9") = {26};

  Physical Point("PointRight10") = {27};
  Physical Point("PointLeft10") = {28};

  Physical Point("PointRight11") = {29};
  Physical Point("PointLeft11") = {30};

  Physical Point("PointRight12") = {31};
  Physical Point("PointLeft12") = {32};



  Coherence Mesh;


