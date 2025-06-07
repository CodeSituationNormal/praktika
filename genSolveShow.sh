#!/bin/bash

cd build

gmshh=0
llt=0
format=2

if [ -f ../input/meshh.geo ]; then
   echo "Файл meshh.geo найден."
   read -p "Хотите ли вы использовать Gmsh для генерации сетки? (y/n): " answer
   if [[ "$answer" == "y" || "$answer" == "Y" ]]; then
      gmshh=1
      read -p "Хотите ли вы использовать format 1 для генерации сетки? (по умолчанию 2) (y/n): " answer2
      if [[ "$answer2" == "y" || "$answer2" == "Y" ]]; then
         format=1
      fi
   fi
   read -p "Хотите ли вы использовать LLt для решения СЛАУ? (по умолчанию МСГ) (y/n): " answer3
   if [[ "$answer3" == "y" || "$answer3" == "Y" ]]; then
      llt=1
   fi
else
   echo "Файл meshh.geo не найден. Использование Gmsh невозможно."
fi

echo $gmshh $llt $format > ../input/settings.txt

if [ "$gmshh" = 1 ]; then
   if [ "$format" -eq 1 ]; then
      echo "Использование формата msh1."
      gmsh ../input/meshh.geo -3 -format msh1 -o ../input/meshh.msh
   else
      echo "Использование формата msh2."
      gmsh ../input/meshh.geo -3 -format msh2 -o ../input/meshh.msh
   fi
   if [ $? -ne 0 ]; then
      echo "Ошибка при запуске Gmsh."
      exit 1
   fi
   read -p "Нажмите любую клавишу для продолжения..."
   ./generate.exe
else
   ./generate.exe
fi

./solve.exe
cd ..
if [ "$gmshh" = false ]; then
   python 3DGrid.py
fi