#!/bin/bash

cd build

gmshh=true

if [ "$gmshh" = true ]; then
   gmsh ../input/meshh.geo -3 -format msh2 -o ../input/meshh.msh
   if [ $? -ne 0 ]; then
      echo "Ошибка при запуске Gmsh."
      exit 1
   fi
   read -p "Нажмите любую клавишу для продолжения..."
else
   ./generate.exe
fi

./solve.exe
cd ..
if [ "$gmshh" = false ]; then
   python 3DGrid.py
fi