#include "common_includes.h"

double u_a(int i, double t) {
   u[i] = sin(nodes[i].x); // modify manually if needed
   return u[i];
}

double u_c(double x, double y, double z, double t) {
   return sin(x); // modify manually if needed
}

double f_auto(double x, double y, double z, double t) {
   return sin(x); // modify manually if needed
}

void dif_u() {
   dif.resize(nodes_c);
   u.resize(nodes_c);
   ofstream difFile("../dif.txt");
   for (int i = 0; i < nodes_c; i++) {
      dif[i] = u_a(i, 0) - q[i];
      difFile << scientific << setprecision(10) << dif[i] << endl;
   }
   difFile.close();
}

void print_u() {
   ofstream uFile("../u.txt");
   // cout << endl << "u ";
   for (int i = 0; i < nodes_c; i++) {
      // cout << u[i] << " ";
      uFile << scientific << setprecision(10) << u[i] << endl;
   }
   // cout << endl;
   uFile.close();
}