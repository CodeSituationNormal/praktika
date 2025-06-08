#include "common_includes.h"

int nx, ny, nz, nt, n;
double x_min, x_max, y_min, y_max, z_min, z_max, t_min, t_max, kx, ky, kz, kt, hx, hy, hz, ht;
vector<int> bc1;

int nodes_c, el_c, face_c, physical_names_c, maxiter = 10000;
double lam = 1, sig = 1, alpha, beta_beta, norma_pr, eps = 1e-18;
vector <node> nodes;
vector <EL> el;
vector <int> ig, jg;
vector<double> u, dif, di, gg, r, b, q, x, y, z, val, Az, Ar, Mr, b_loc;
vector<vector<double>> A_loc, M_loc, G_loc;

vector<double> q1, q2, q3;
vector<double> t;
int times_c;
double current_t, t_1, t_2, t_3, c_0, c_1, c_2, c_3, eta_0, eta_1, eta_2, eta_3; 

vector <pair<int, int>> faces; // Node num, physical tag
vector <pair<int, string>> physical_names;

int format;