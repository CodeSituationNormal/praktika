// Параметры цилиндра
radius = 1.0;
height = 3.0;

// 1. Создаем базовые точки (нижнее основание)
Point(1) = {0, 0, 0, 1.0};  // Центр
Point(2) = {radius, 0, 0, 1.0};  // Ось X
Point(3) = {0, radius, 0, 1.0};  // Ось Y
Point(4) = {-radius, 0, 0, 1.0};
Point(5) = {0, -radius, 0, 1.0};

// 2. Нижняя окружность (4 сегмента)
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// 3. Верхние точки
Point(6) = {radius, 0, height, 1.0};
Point(7) = {0, radius, height, 1.0};
Point(8) = {-radius, 0, height, 1.0};
Point(9) = {0, -radius, height, 1.0};
Point(10) = {0, 0, height, 1.0};  // Центр верхнего основания

// 4. Вертикальные линии
Line(5) = {2, 6};
Line(6) = {3, 7};
Line(7) = {4, 8};
Line(8) = {5, 9};

// 5. Верхняя окружность
Circle(9) = {6, 10, 7};
Circle(10) = {7, 10, 8};
Circle(11) = {8, 10, 9};
Circle(12) = {9, 10, 6};

// 6. Боковые поверхности
Line Loop(1) = {1, 6, -9, -5};
Surface(1) = {1};
Line Loop(2) = {2, 7, -10, -6};
Surface(2) = {2};
Line Loop(3) = {3, 8, -11, -7};
Surface(3) = {3};
Line Loop(4) = {4, 5, -12, -8};
Surface(4) = {4};

// 7. Основания
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(5) = {5};
Line Loop(6) = {9, 10, 11, 12};
Plane Surface(6) = {6};

// 8. Объем
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// 9. Трансфинитное разбиение
Transfinite Curve {1, 2, 3, 4, 9, 10, 11, 12} = 17;  // Окружности
Transfinite Curve {5, 6, 7, 8} = 17;  // Вертикальные линии

// 10. Поверхности
Transfinite Surface {1, 2, 3, 4};
Recombine Surface {1, 2, 3, 4};
Transfinite Surface {5, 6};
Recombine Surface {5, 6};

// 11. Объем (правильный синтаксис)
Transfinite Volume{1} = {2, 3, 4, 5, 6, 7, 8, 9};
Recombine Volume{1};

// 12. Физические группы
Physical Surface("Bottom") = {5};
Physical Surface("Top") = {6};
Physical Surface("Sides") = {1, 2, 3, 4};
Physical Volume("Cylinder") = {1};

// 13. Настройки сетки
Mesh.Algorithm3D = 4;
Mesh.RecombinationAlgorithm = 1;
Mesh.RecombineAll = 1;
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;