#include "common_includes.h"

void input_nodes(int testNumber) {
   // char filename[50];
   // sprintf_s(filename, "nodes%d.txt", testNumber);
   string filename = "../nodes_out.txt";

   ifstream nodes_f(filename);
   if (!nodes_f.is_open()) {
      perror("No such file");
      return;
   }

   nodes.clear();
   while (!nodes_f.eof()) {
      node n;
      nodes_f >> n.x >> n.y >> n.z;
      if (nodes_f.fail()) break; 
      n.number = nodes.size();
      nodes.push_back(n);
   }
   
   nodes_f.close();
   nodes_c = nodes.size(); 
}

void input_el(int testNumber) {
   // char filename[50];
   // sprintf_s(filename, "el%d.txt", testNumber);
   string filename = "../elements_out.txt";

   ifstream el_f(filename);
   if (!el_f.is_open()) {
      perror("No such file");
      return;
   }

   el.clear();
   while (!el_f.eof()) {
      EL e;
      for (int i = 0; i < 8; i++) {
         el_f >> e.node_n[i];
         if (el_f.fail()) break; 
      }
      if (el_f.fail()) break; 
      e.number = el.size();
      el.push_back(e);
   }
   el_f.close();
   el_c = el.size();
}

void input_faces(int testNumber) {
   // char filename[50];
   // sprintf_s(filename, "face%d.txt", testNumber);
   string filename = "../faces.txt";

   ifstream face_f(filename);
   if (!face_f.is_open()) {
      perror("No such file");
      return;
   }
   
   faces.clear();
   while(!face_f.eof()) {
      int face_node;
      // double value;
      // face_f >> face_node >> value;
      face_f >> face_node;
      if (face_f.fail()) break; 
      faces.push_back(face_node);
      // val.push_back(value);
   }

   face_f.close();
   face_c = faces.size();
}

void input_el_coef(int testNumber) {
   // char filename[50];
   // sprintf_s(filename, "coef%d.txt", testNumber); 

   string filename = "../coef.txt";

   ifstream el_coef_f(filename);
   if (!el_coef_f.is_open()) {
      perror("No such file");
      return;
   }
   el_coef_f >> lam >> sig;
   if (el_coef_f.fail()) {
      perror("Error reading file");
      return;
   }
   el_coef_f.close();
}

void input_f(int testNumber) {
   // char filename[50];
   // sprintf_s(filename, "f%d.txt", testNumber);

   string filename = "../f.txt";

   ifstream f_f(filename);
   if (!f_f.is_open()) {
      perror("No such file");
      return;
   }
   
   for (int i = 0; i < nodes_c; i++) 
      f_f >> nodes[i].f;
   if (f_f.fail()) {
      perror("Error reading file");
      return;
   }
   f_f.close();
} 

void input_t(int testNumber) {
   // char filename[50];
   // sprintf_s(filename, "f%d.txt", testNumber);

   string filename = "../time.txt";

   ifstream t_t(filename);
   if (!t_t.is_open()) {
      perror("No such file");
      return;
   }
   
   while(!t_t.eof()) {
      double time_value;
      t_t >> time_value;
      if (t_t.fail()) break; 
      t.push_back(time_value);
   }

   times_c = t.size();
   t_min = t[0];
   t_max = t[times_c - 1];
   ht = t[1] - t[0];
   
   t_t.close();
} 
