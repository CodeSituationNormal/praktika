SetFactory("OpenCASCADE");

// Параметры
L = 2;
h = 1;
n = L / h + 1; // 3 узла по каждой координате

// Куб
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};
Point(5) = {0, 0, L, h};
Point(6) = {L, 0, L, h};
Point(7) = {L, L, L, h};
Point(8) = {0, L, L, h};

// Рёбра
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Поверхности
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};

Line Loop(15) = {5, 6, 7, 8};
Plane Surface(16) = {15};

Line Loop(17) = {1, 10, -5, -9};
Plane Surface(18) = {17};

Line Loop(19) = {2, 11, -6, -10};
Plane Surface(20) = {19};

Line Loop(21) = {3, 12, -7, -11};
Plane Surface(22) = {21};

Line Loop(23) = {4, 9, -8, -12};
Plane Surface(24) = {23};
 
// Объём
Surface Loop(25) = {14, 16, 18, 20, 22, 24};
Volume(26) = {25};

// Трансфинитная сетка — это сетка, построенная по заранее заданному числу узлов/элементов вдоль краевых линий поверхности, с регулярным (упорядоченным) расположением узлов внутри.
// разбить каждую из линий от 1 до 12 на n узлов.
Transfinite Line {1:12} = n;
Transfinite Surface {14, 16, 18, 20, 22, 24};
Transfinite Volume {26};
// Определение физических поверхностей с именами (краевые условия)
Physical Surface("bc1", 30) = {14, 16,18,20,22,24};  // 
//Physical Surface("bc1", 30) = {22};  
//Physical Surface("bc2", 31) = {20, 24};  
//Physical Surface("bc3", 32) = {14, 16, 18};  




Recombine Volume {26}; // превращает гексаэдры в "рекомбинированные" гексаэдры (для разбиения)
Physical Volume("Cube", 33) = {26};
//Mesh.Recombine3DAll = 1; // рекомбинация всех 3D объёмов

Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
//Mesh 3;

Mesh.Algorithm3D = 4;
Mesh.RecombinationAlgorithm = 1;
Mesh.RecombineAll = 1;
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;