#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
using namespace std;

struct node {
   int number;
   double x, y, z;
   double f;
};
struct EL {
   int number;
   int physical_tag;
   int elementary_tag;
   int node_n[8]{};
   double hx = 0, hy = 0, hz = 0;
};


extern int nx, ny, nz, nt, n;
extern double x_min, x_max, y_min, y_max, z_min, z_max, t_min, t_max, kx, ky, kz, kt, hx, hy, hz, ht;
extern vector<int> bc1;


extern int nodes_c, el_c, face_c, physical_names_c, maxiter;
extern double lam, sig, alpha, beta_beta, norma_pr, eps;
extern vector <node> nodes;
extern vector <EL> el;
extern vector <int> ig, jg;
extern vector<double> u, dif, di, gg, r, b, q, x, y, z, val, Az, Ar, Mr, b_loc;
extern vector<vector<double>> A_loc, M_loc, G_loc;

extern vector<double> q1, q2, q3;
extern vector<double> t; 
extern int times_c;
extern double current_t, t_1, t_2, t_3, c_0, c_1, c_2, c_3, eta_0, eta_1, eta_2, eta_3;

extern vector <pair<int, int>> faces; // Node num, tag
extern vector <pair<int, string>> physical_names; // Physical tag, name

extern int format;

#endif // STRUCTURES_H