#include "common_includes.h"

int t_it = 3;
ofstream outf, outdif;

double sigma() {
   return 1;
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
      val[i] = u_c(nodes[faces[i]].x, nodes[faces[i]].y, nodes[faces[i]].z, current_t);
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

static void calc_h() {
	for (int i = 0; i < el_c; i++) {
		int node_min_index = el[i].node_n[0];
		int node_max_index = el[i].node_n[7];

		el[i].hx = nodes[node_max_index].x - nodes[node_min_index].x;
		el[i].hy = nodes[node_max_index].y - nodes[node_min_index].y;
		el[i].hz = nodes[node_max_index].z - nodes[node_min_index].z;
	}
}

static void local_el(int f_el_n) {
	double hx = el[f_el_n].hx;
	double hy = el[f_el_n].hy;
	double hz = el[f_el_n].hz;

	double dx2 = 1 / hx;
	double x2 = hx / 3;

	double dy2 = 1 / hy;
	double y2 = hy / 3;

	double dz2 = 1 / hz;
	double z2 = hz / 3;

	double** gx = nullptr;
	double** gy = nullptr;
	double** gz = nullptr;

   A_loc.resize(8, vector<double>(8));
   b_loc.clear();
   b_loc.resize(8, 0);
   M_loc.resize(8, vector<double>(8));
   G_loc.resize(8, vector<double>(8));
	gx = new double* [8] ();
	gy = new double* [8] ();
	gz = new double* [8] ();

	for (int i = 0; i < 8; i++) {
		gx[i] = new double[8] ();
		gy[i] = new double[8] ();
		gz[i] = new double[8] ();  
	}
	double coefXd = 0, coefYd = 0, coefZd = 0;
	double coefX = 0, coefY = 0, coefZ = 0;
	for (int i = 0; i < 8; i++) {
      double Mq1 = 0, Mq2 = 0, Mq3 = 0;
		for (int j = 0; j < 8; j++) {
			if ((i % 2 != 0 && j % 2 != 0) || (i % 2 == 0 && j % 2 == 0)) coefXd = dx2, coefX = x2;
			else coefXd = -dx2, coefX = x2 / 2;
			if ((i < 4 && j < 4) || (i >= 4 && j >= 4)) coefZd = dz2, coefZ = z2;
			else coefZd = -dz2, coefZ = z2 / 2;
			if (((i == 1 || i == 0 || i == 4 || i == 5) && (j == 1 || j == 0 || j == 4 || j == 5)) || ((i != 1 && i != 0 && i != 4 && i != 5) && (j != 1 && j != 0 && j != 4 && j != 5))) coefYd = dy2, coefY = y2;
			else coefYd = -dy2, coefY = y2 / 2;

			gx[i][j] = coefXd * coefY * coefZ;
			gy[i][j] = coefX * coefYd * coefZ;
			gz[i][j] = coefX * coefY * coefZd;

			M_loc[i][j] = (hx * hy * hz) * coefX * coefY * coefZ; 
			G_loc[i][j] = (hx * hy * hz) * (gx[i][j] / (hx * hx) + gy[i][j] / (hy * hy) + gz[i][j] / (hz * hz));

			A_loc[i][j] = lam * G_loc[i][j] + sigma() * M_loc[i][j] * c_0;
			b_loc[i] += nodes[el[f_el_n].node_n[j]].f * M_loc[i][j]; 

         Mq1 += M_loc[i][j] * q1[el[f_el_n].node_n[j]];
         Mq2 += M_loc[i][j] * q2[el[f_el_n].node_n[j]];
         Mq3 += M_loc[i][j] * q3[el[f_el_n].node_n[j]];
      }
      b_loc[i] -= sigma() * (c_1 * Mq1 + c_2 * Mq2 + c_3 * Mq3);
	}
	for (int i = 0; i < 8; i++)
	{
		delete[] gx[i];
		delete[] gy[i];
		delete[] gz[i];
	}
	delete[] gx;
	delete[] gy;
	delete[] gz;
   M_loc.clear();
   G_loc.clear();
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
void global_A() {
   int i_gg, glob_i, glob_j;

   A_loc.resize(8, vector<double>(8, 0));
   b_loc.resize(8, 0);

   di.resize(nodes_c, 0);
   gg.resize(ig[nodes_c], 0);
   b.resize(nodes_c, 0);

   // Should be generated every time step
   generate_bc1();
   generate_f();
   calc_c();

   for (int k = 0; k < el_c; k++) {
      local_el(k);
      for (int i = 0; i < 8; i++) {
         glob_i = el[k].node_n[i];
         for (int j = 0; j < 8; j++) {
            glob_j = el[k].node_n[j];
            if (glob_i > glob_j) {
               i_gg = find_el_pos(glob_i, glob_j);
               gg[i_gg] += A_loc[i][j];
            }
         }
         di[glob_i] += A_loc[i][i];
         b[glob_i] += b_loc[i];
      }
      A_loc.clear();
      b_loc.clear();
   }
   for (int j = 0; j < face_c; j++) {
      glob_i = faces[j];
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
   }

   // double norma_dif = sqrt(scMult(dif, dif));
   outdif << current_t << " " << scientific << setprecision(10) << dif[13] << endl;
   // cout << norma_dif << endl;
   // outf << norma_dif << " "; //????
   dif.clear();
   u.clear();

   // qFile.close();
   r.clear(); z.clear(); Az.clear();
}

void fourth_order_temporal_scheme() {
   generate_q123();

   // cout << "t = " << t[0] << endl;
   // for (int i = 0; i < nodes_c; i++)
   //    cout << q3[i] << " ";
   // cout << endl;
   // cout << "t = " << t[1] << endl;
   // for (int i = 0; i < nodes_c; i++)
   //    cout << q2[i] << " ";
   // cout << endl;
   // cout << "t = " << t[2] << endl;
   // for (int i = 0; i < nodes_c; i++)
   //    cout << q1[i] << " ";
   // cout << endl;

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
      global_A();
      CGM();

      // cout << "t = " << current_t << endl;
      outf << current_t << " ";
      for (int i = 0; i < nodes_c; i++){
         // cout << q[i] << " ";
         outf << scientific << setprecision(10) << q[i] << " ";
         q3[i] = q2[i]; 
         q2[i] = q1[i];
         q1[i] = q[i];
      }
      // cout << endl;
      outf << endl;

      t_it++;
      current_t = t[t_it];
      gg.clear();
      di.clear();
      b.clear();
   }
}

int main() {
   // eps = 1e-14;   
   // maxiter = 10000;
   outf.open("../output.txt");
   outdif.open("../dif.txt");
   int testNumber = 0;

   // cout << "Enter the test number: ";
   // cin >> testNumber;

   input_nodes(testNumber);
   input_el(testNumber);
   input_faces(testNumber);
   input_el_coef(testNumber);
   // input_f(testNumber);
   input_t(testNumber);
   
   portrait();
   calc_h();

   fourth_order_temporal_scheme();

   // global_A();
   // CGM();
   // dif_u();
   // print_u();

   outf.close();
   outdif.close();

   return 0;
}