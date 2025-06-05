#include "common_includes.h"

int t_it = 3;
ofstream outf, outdif;

double sigma() {
   return 1;
}

double lambda() {
   return 1;
}

vector<vector<double>> M_quad, Mz, M_hex, A_hex;  
vector<vector<double>> G_quad, Gz, G_hex;  
vector<double> b_hex;
double s, e;
vector<double> N, dN_ds, dN_de;

void shape_functions() {
   N.resize(4), dN_ds.resize(4), dN_de.resize(4);
   N[0] = (1.0 - s) * (1.0 - e);
   N[1] = s * (1.0 - e);
   N[2] = s * e;
   N[3] = (1.0 - s) * e;

   dN_ds[0] = -(1.0 - e);
   dN_ds[1] =  (1.0 - e);
   dN_ds[2] =  e;
   dN_ds[3] = -e;

   dN_de[0] = -(1.0 - s);
   dN_de[1] = -s;
   dN_de[2] =  s;
   dN_de[3] =  (1.0 - s);
}

void loc_quad(int f_el_n) {
   vector<double> x_lc, y_lc;
   x_lc.resize(4), y_lc.resize(4);
   for (int i = 0; i < 4; i++) {
      x_lc[i] = nodes[el[f_el_n].node_n[i]].x;
      y_lc[i] = nodes[el[f_el_n].node_n[i]].y;
   }

   M_quad = vector<vector<double>>(4, vector<double>(4, 0.0));
   G_quad = vector<vector<double>>(4, vector<double>(4, 0.0));

   const int nGauss = 2;
   const double sqrt3 = 1.0 / sqrt(3.0);
   double gaussPoints[nGauss] = { 0.5 * (1.0 - sqrt3), 0.5 * (1.0 + sqrt3) };
   double gaussWeights[nGauss] = {0.5, 0.5};

   for (int iq = 0; iq < nGauss; iq++) {
      s = gaussPoints[iq];
      double ws = gaussWeights[iq];
      for (int jq = 0; jq < nGauss; jq++) {
         e = gaussPoints[jq];
         double we = gaussWeights[jq];
         double weight = ws * we;

         shape_functions();

         double dx_ds = 0.0, dx_de = 0.0;
         double dy_ds = 0.0, dy_de = 0.0;
         for (int k = 0; k < 4; k++) {
               dx_ds += x_lc[k] * dN_ds[k];
               dx_de += x_lc[k] * dN_de[k];
               dy_ds += y_lc[k] * dN_ds[k];
               dy_de += y_lc[k] * dN_de[k];
         }

         double detJ = dx_ds * dy_de - dx_de * dy_ds;
         if (detJ <= 0.0) {
            detJ = detJ * (-1.0);
         }

         double invJ[2][2] = {
               {  dy_de / detJ, -dy_ds / detJ },
               { -dx_de / detJ,  dx_ds / detJ }
         };

         double dN_dx[4], dN_dy[4];
         for (int k = 0; k < 4; k++) {
               dN_dx[k] = invJ[0][0] * dN_ds[k] + invJ[0][1] * dN_de[k];
               dN_dy[k] = invJ[1][0] * dN_ds[k] + invJ[1][1] * dN_de[k];
         }

         for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
               M_quad[i][j] += weight * N[i] * N[j] * detJ;
               G_quad[i][j] += weight * (dN_dx[i] * dN_dx[j] + dN_dy[i] * dN_dy[j]) * detJ;
            }
      }
   }
}

int mu(int index) {
   return index % 4;
}

int nu(int index) {
   int res = index / 4;
   return res;
}

void loc_z(int f_el_n) {
   int node_min_index = el[f_el_n].node_n[0];
   int node_max_index = el[f_el_n].node_n[7];

   el[f_el_n].hz = nodes[node_max_index].z - nodes[node_min_index].z;   

   double coef_g = 1 / el[f_el_n].hz;
   
   Gz.resize(2,vector<double>(2));
   Gz[0][0] = Gz[1][1] = coef_g;
   Gz[1][0] = Gz[0][1] = -coef_g;

   double coef_m = el[f_el_n].hz / 6;
   
   Mz.resize(2,vector<double>(2));
   Mz[0][0] = Mz[1][1] = 2 * coef_m;
   Mz[1][0] = Mz[0][1] = coef_m;
}

void loc_hex(int f_el_n) {
   loc_quad(f_el_n);
   loc_z(f_el_n);

   M_hex.resize(8, vector<double>(8,0));
   G_hex.resize(8, vector<double>(8,0));
   A_hex.resize(8, vector<double>(8,0));
   b_hex.clear();
   b_hex.resize(8,0);

   for (int i = 0; i < 8; i++) {
      double Mq1 = 0, Mq2 = 0, Mq3 = 0;
      for (int j = 0; j < 8; j++) {
         M_hex[i][j] = M_quad[mu(i)][mu(j)] * Mz[nu(i)][nu(j)];
         G_hex[i][j] = G_quad[mu(i)][mu(j)] * Mz[nu(i)][nu(j)] + M_quad[mu(i)][mu(j)] * Gz[nu(i)][nu(j)];

         A_hex[i][j] = sigma() * M_hex[i][j] * c_0 + lambda() * G_hex[i][j];
			b_hex[i] += nodes[el[f_el_n].node_n[j]].f * M_hex[i][j]; 

         Mq1 += M_hex[i][j] * q1[el[f_el_n].node_n[j]];
         Mq2 += M_hex[i][j] * q2[el[f_el_n].node_n[j]];
         Mq3 += M_hex[i][j] * q3[el[f_el_n].node_n[j]];
      }
      b_hex[i] -= sigma() * (c_1 * Mq1 + c_2 * Mq2 + c_3 * Mq3);
   }

   M_quad.clear(), G_quad.clear(), Mz.clear(), Gz.clear();
}

void generate_q123() {
   q1.resize(nodes_c);
   q2.resize(nodes_c);
   q3.resize(nodes_c);
   for (int i = 0; i < nodes_c; i++) {
      q1[i] = u_c(nodes[i].x, nodes[i].y, nodes[i].z, t[2]);
      q2[i] = u_c(nodes[i].x, nodes[i].y, nodes[i].z, t[1]);
      q3[i] = u_c(nodes[i].x, nodes[i].y, nodes[i].z, t[0]);
   }
   current_t = t[3];
   t_1 = t[2];
   t_2 = t[1];
   t_3 = t[0];
}

void generate_bc1() {
   val.resize(face_c);
   for (int i = 0; i < face_c; i++) {
      val[i] = u_c(nodes[faces[i].first].x, nodes[faces[i].first].y, nodes[faces[i].first].z, current_t);
   }
}
void generate_f() {
   for (int i = 0; i < nodes_c; i++) {
      nodes[i].f = f_auto(nodes[i].x, nodes[i].y, nodes[i].z, current_t);
   }
}

void portrait() {
   map<int, set<int>> list;
   for (int i = 0; i < el_c; ++i) {
      for (int j = 0; j < 8; ++j) {
         int node_index = el[i].node_n[j];
         list[node_index].insert(node_index);
         for (int k = 0; k < 8; ++k) {
            int neighbor_index = el[i].node_n[k];
            if (node_index != neighbor_index) {
               list[node_index].insert(neighbor_index);
            }
         }
      }
   }
   
   // for (size_t i = 0; i < list.size(); ++i) {
   //    cout << "Node " << i << ": ";
   //    for (int node : list[i]) 
   //       cout << node << " ";
   //    cout << endl;
   // }
   int k = 0;
   ig.push_back(0);
   for (int i = 0; i < nodes_c; i++) {
      k = 0;
      for (int node : list[i]) {
         if (node < i) {
            jg.push_back(node);
            k++;
         }
      }
      ig.push_back(ig.back() + k);
   }
   // cout << "jg " << ": ";
   // for (int node : jg) {
   //    cout << node << " ";
   // }
   // cout << endl;
   // cout << "ig " << ": ";
   // for (int node : ig) {
   //    cout << node << " ";
   // }
   // cout << endl;
}
void calc_c() {
   double t0 = current_t, t1 = t[t_it - 1], t2 = t[t_it - 2], t3 = t[t_it - 3], d01, d02, d03, d12, d13, d23, t0p2;

   d01 = t0 - t1;
   d02 = t0 - t2;
   d03 = t0 - t3;
   d12 = t1 - t2;
   d13 = t1 - t3;
   d23 = t2 - t3;

   t0p2 = t0 * t0;

   // cout << " times " << t0 << " " << t1 << " " << t2 << " " << t3 << endl;

   c_0 = (3 * t0p2 - 2 * t0 * (t1 + t2 + t3) + (t1 * t2) + (t1 * t3) + (t2 * t3)) / (d01 * d02 * d03);
   // cout << "nu0: " << c_0 << endl;

   c_1 = (t0p2 - t0 * t2 - t3 * t0 + t3 * t2) / (d13 * d12 * (-d01));
   // cout <<"nu1: " << c_1 << endl;

   c_2 = (t0p2 - t1 * t0 - t3 * t0 + t3 * t1) / (d23 * (-d12) * (-d02));
   // cout << "nu2: " << c_2 << endl;

   c_3 = (t0p2 - t1 * t0 - t2 * t0 + t2 * t1) / ((-d23) * (-d13) * (-d03));
   // cout << "nu3: " << c_3 << endl;
}
int find_el_pos(int i, int j) {
   int s = ig[i], e = ig[i + 1];
   int cur = 0;
   bool find = false;
   for (int p = s; p < e && !find; p++) {
      if (jg[p] == j) {
         cur = p;
         find = true;
      }
   }
   return cur;
}
void global_A_hex() {
   int i_gg, glob_i, glob_j;

   A_hex.resize(8, vector<double>(8, 0));
   b_hex.resize(8, 0);

   di.resize(nodes_c, 0);
   gg.resize(ig[nodes_c], 0);
   b.resize(nodes_c, 0);

   // Should be generated every time step
   generate_bc1();
   generate_f();

   calc_c();

   for (int k = 0; k < el_c; k++) {
      loc_hex(k);
      for (int i = 0; i < 8; i++) {
         glob_i = el[k].node_n[i];
         for (int j = 0; j < 8; j++) {
            glob_j = el[k].node_n[j];
            if (glob_i > glob_j) {
               i_gg = find_el_pos(glob_i, glob_j);
               gg[i_gg] += A_hex[i][j];
            }
         }
         di[glob_i] += A_hex[i][i];
         b[glob_i] += b_hex[i];
      }
      A_hex.clear();
      M_hex.clear();
      G_hex.clear();
      b_hex.clear();
   }
   for (int j = 0; j < face_c; j++) {
      glob_i = faces[j].first;
      di[glob_i] = 1;
      b[glob_i] = val[j];
      for (int i = ig[glob_i]; i < ig[glob_i + 1]; i++) {
         b[jg[i]] -= gg[i] * val[j];
         gg[i] = 0;
      }
      for (int p = glob_i + 1; p < nodes_c; p++)
         for (int i = ig[p]; i < ig[p + 1]; i++) {
            if (jg[i] == glob_i) {
               b[p] -= gg[i] * val[j];
               gg[i] = 0;
            }
         }
   }
   // cout << endl;
   // for (int j = 0; j < face_c; j++) cout << faces[j] << endl;
   // cout << "di " << ": ";
   // for (double node : di)
   //    cout << node << " ";
   // cout << "gg " << ": ";
   // for (double node : gg)
   //    cout << node << " ";
   // cout << endl;
   // cout << "b " << ": ";
   // for (double node : b)
   //    cout << node << " ";
   // cout << endl;

   val.clear();
}

void trans_Gauss(vector<vector<double>> A, int i, int n) {
   int line = i;
   for (int j = i + 1; j < n; j++)
      if (fabs(A[j][i]) > fabs(A[line][i]))
         line = j;
   if (line != i) 
      for (int j = 0; j < 2 * n; j++)
         swap(A[i][j], A[line][j]);
}
void GaussJordan(vector<vector<double>> A, vector<vector<double>> alpha, int n) {
   vector<vector<double>> matrix;
   matrix.resize(n);
   for (int i = 0; i < n; i++) {
      matrix[i].resize(2 * n, 0.0);
      for (int j = 0; j < n; j++) 
         matrix[i][j] = A[i][j];
      for (int j = 0; j < n; j++) 
         matrix[i][j + n] = (j == i) ? 1.0 : 0.0;
   }

   for (int i = 0; i < n; i++) {
      trans_Gauss(matrix, i, n);
      double m_d = matrix[i][i];
      for (int j = 0; j < 2 * n; j++)
         matrix[i][j] /= m_d;
      for (int j = 0; j < n; j++) {
         if (j != i) {
            double m_j = matrix[j][i];
            for (int k = 0; k < 2 * n; k++)
               matrix[j][k] -= matrix[i][k] * m_j;
         }
      }
   }
   for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) 
         alpha[i][j] = matrix[i][j + n]; 
   }
   for (int i = 0; i < n; i++) matrix[i].clear();
   matrix.clear();
}


static void calc_Av(vector<double>& v, vector<double>& res) {
   for (int i = 0; i < nodes_c; i++) {
      res[i] = di[i] * v[i];
      for (int k = ig[i]; k < ig[i + 1]; k++) {
         int j = jg[k];
         res[i] += gg[k] * v[j];
         res[j] += gg[k] * v[i];
      }
   }
}
static void calc_r0() {
   for (int i = 0; i < nodes_c; i++) {
      r[i] = b[i] - di[i] * q[i];
      for (int k = ig[i]; k < ig[i + 1]; k++) {
         int j = jg[k];
         r[i] -= gg[k] * q[j];
         r[j] -= gg[k] * q[i];
      }
   }
}
static void calc_x() {
   for (int i = 0; i < nodes_c; i++) q[i] += alpha * z[i];
}
static void calc_r(vector<double>& x) {
   for (int i = 0; i < nodes_c; i++) r[i] -= alpha * x[i];
}
static void calc_z(vector<double>& x) {
   for (int i = 0; i < nodes_c; i++) z[i] = x[i] + beta_beta * z[i];
}
static double scMult(vector<double>& x, vector<double>& y) {
   double res = 0;
   for (int i = 0; i < nodes_c; i++) res += x[i] * y[i];
   return res;
}
static vector<double> vecMult(vector<double>& x, vector<double>& y, vector<double>& res) {
   for (int i = 0; i < nodes_c; i++) res[i] = x[i] * y[i];
   return res;
}
static void CGM() {
   double nev = 1;
   q.resize(nodes_c), r.resize(nodes_c), z.resize(nodes_c), Az.resize(nodes_c), Ar.resize(nodes_c), Mr.resize(nodes_c);
   for (int i = 0; i < nodes_c; i++) q[i] = 0;  
 
   norma_pr = sqrt(scMult(b, b));
   calc_r0();
   double r_scMult = scMult(r, r);
   z = r;
   for (int k = 0; k < maxiter && nev > eps; k++) {
      double r_scMult_old = r_scMult;
      calc_Av(z, Az);
      double zAz_scMult = scMult(Az, z);
      alpha = r_scMult_old / zAz_scMult;
      calc_x();
      calc_r(Az);
      r_scMult = scMult(r, r);
      nev = sqrt(r_scMult) / norma_pr;
      beta_beta = r_scMult / r_scMult_old;
      calc_z(r);
   }
   // ofstream qFile("../q.txt");
  
   // cout << "q ";
   // for (int i = 0; i < nodes_c; i++) {
   //    // cout << q[i] << " ";
   //    qFile << scientific << setprecision(10) << q[i] << endl;
   // }

   dif.resize(nodes_c);
   u.resize(nodes_c);
   for (int i = 0; i < nodes_c; i++) {
      dif[i] = u_a(i, current_t) - q[i];

      // if((nodes[i].x == 0) && (nodes[i].y == 0) && (nodes[i].z == 1.5)) 
      //    outdif << current_t << " " << scientific << setprecision(10) << dif[i] << endl;

      // if(i == 111) 
         // outdif << current_t << " " << scientific << setprecision(10) << dif[i] << endl;
   }

   double norma_dif = sqrt(scMult(dif, dif));

   // outdif << current_t << " " << scientific << setprecision(10) << dif[13] << endl;

   cout << norma_dif << endl;

   // outf << norma_dif << " "; //????
   dif.clear();
   u.clear();

   // qFile.close();
   r.clear(); z.clear(); Az.clear();
}

void fourth_order_temporal_scheme() {
   generate_q123();

   outf << scientific << setprecision(10) << t[0] << " ";
   for (int i = 0; i < nodes_c; i++)
      outf << scientific << setprecision(10) << q3[i] << " ";
   outf << endl;
   outf << t[1] << " ";
   for (int i = 0; i < nodes_c; i++)
      outf << scientific << setprecision(10) << q2[i] << " ";
   outf << endl;
   outf << t[2] << " ";
   for (int i = 0; i < nodes_c; i++)
      outf << scientific << setprecision(10) << q1[i] << " ";
   outf << endl;
   
   while(t_it != times_c){
      global_A_hex();
      CGM();

      outf << current_t << " ";
      for (int i = 0; i < nodes_c; i++){
         outf << scientific << setprecision(10) << q[i] << " ";
         q3[i] = q2[i]; 
         q2[i] = q1[i];
         q1[i] = q[i];
      }
      outf << endl;

      t_it++;
      current_t = t[t_it];
      gg.clear();
      di.clear();
      b.clear();
   }
}

void input_all(bool is_gmsh) {
   if (is_gmsh) {
      if (format == 1) input_gmsh_1("../input/meshh.msh");
      else if (format == 2) input_gmsh_2("../input/meshh.msh");
      else {
         cerr << "Unknown mesh format" << endl;
         exit(1);
      }
      input_t();
   } 
   else {
      input_nodes();
      input_el();
      input_faces();
      input_el_coef();
      input_t();
   }
}

void call_functions(bool is_static) {
   if (is_static) {
      global_A_hex();
      CGM();
      dif_u();
      print_u();
   }
   else {
      fourth_order_temporal_scheme();
   }
}

int main() {
   outf.open("../output/output.txt");
   outdif.open("../output/dif.txt");

   format = 2;

   input_all(false); // Change to true for GMSH input and false for built-in input
   portrait();

   call_functions(false); // Change to true for static case and false for dynamic case

   outf.close();
   outdif.close();

   return 0;
}

//LU LLt
//t