#include "common_includes.h"

void input_nodes() {
   string filename = "../processFiles/nodes_out.txt";

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

void input_el() {
   string filename = "../processFiles/elements_out.txt";

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

void input_faces() {
   string filename = "../processFiles/faces.txt";

   ifstream face_f(filename);
   if (!face_f.is_open()) {
      perror("No such file");
      return;
   }
   
   faces.clear();
   while(!face_f.eof()) {
      int face_node, face_tag;
      face_f >> face_node >> face_tag;
      pair<int, int> face_node_tag;
      face_node_tag.first = face_node;
      face_node_tag.second = face_tag;
      if (face_f.fail()) break; 
      faces.push_back(face_node_tag);
   }

   face_f.close();
   face_c = faces.size();
}

void input_el_coef() {
   string filename = "../input/coef.txt";

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

void input_f() {
   string filename = "../processFiles/f.txt";

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

void input_t() {
   string filename = "../processFiles/time.txt";

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

void input_gmsh_1(const string& filename) {
   ifstream file(filename);
   if (!file.is_open()) {
      cerr << "Error opening file: " << filename << endl;
      return;
   }

   string line;
   bool in_nodes = false;
   bool in_elements = false;

   while (getline(file, line)) {
      // Удаляем лишние пробелы в начале/конце строки
      line.erase(0, line.find_first_not_of(" \t"));
      line.erase(line.find_last_not_of(" \t") + 1);

      if (line.empty()) continue;

      if (line == "$NOD") {
         in_nodes = true;
         in_elements = false;
         
         // Читаем количество узлов
         getline(file, line);
         int num_nodes = stoi(line);
         nodes.resize(num_nodes);
         
         // Читаем узлы
         for (int i = 0; i < num_nodes; ++i) {
            getline(file, line);
            istringstream iss(line);
            int id;
            double x, y, z;
            iss >> id >> x >> y >> z;
            nodes[id - 1] = {id - 1, x, y, z}; // преобразуем к 0-based индексации
         }
      }
      else if (line == "$ENDNOD") {
         in_nodes = false;
      }
      else if (line == "$ELM") {
         in_elements = true;
         in_nodes = false;
         
         // Читаем количество элементов
         getline(file, line);
         int num_elements = stoi(line);
         
         // Читаем элементы
         for (int i = 0; i < num_elements; ++i) {
            getline(file, line);
            istringstream iss(line);
            
            EL elem;
            int num_nodes;
            
            // Формат: номер тип физич_тег элем_тег количество_узлов узлы...
            int type;
            iss >> elem.number >> type >> elem.physical_tag >> elem.elementary_tag >> num_nodes;
            
            // Читаем узлы элемента
            for (int j = 0; j < num_nodes; ++j) {
               int node_id;
               iss >> node_id;
               if (j < 8) { // сохраняем только первые 8 узлов (для гексов)
                  elem.node_n[j] = node_id - 1; // преобразуем к 0-based
               }
            }
               
            if(type == 5) el.push_back(elem);
            
            // Если это поверхностный элемент (тип 3), добавляем в faces
            if (type == 3) {
               for (int i = 0; i < 4; i++) {
                  pair<int, int> face;
                  face.first = elem.node_n[i];
                  face.second = elem.physical_tag;
                  faces.push_back(face);
               }
            }
         }
      }
      else if (line == "$ENDELM") {
         in_elements = false;
      }
   }

   file.close();

   nodes_c = nodes.size();
   el_c = el.size();
   face_c = faces.size();
}

void input_gmsh_2(const string& filename) {
   ifstream file_g2(filename);
   if (!file_g2.is_open()) {
      cerr << "Error opening file: " << filename << endl;
      return;
   }

   string line;

   while (getline(file_g2, line)) {
      // Удаляем лишние пробелы в начале/конце строки
      line.erase(0, line.find_first_not_of(" \t"));
      line.erase(line.find_last_not_of(" \t") + 1);

      if (line.empty()) continue;

      if (line == "$MeshFormat") {
         // Пропускаем строку с форматом
         getline(file_g2, line);
         continue;
      } 
      else if (line == "$EndMeshFormat") {
         // Конец формата
         continue;
      } 
      else if (line == "$PhysicalNames") {
         getline(file_g2, line);
         physical_names_c = stoi(line);
         for (int i = 0; i < physical_names_c; ++i) {
            getline(file_g2, line);
            istringstream iss(line);
            int tag, dim;
            string name;
            iss >> dim >> tag >> name;
            if (iss.fail()) {
               cerr << "Error reading physical name data: " << line << endl;
               continue;
            } 
            pair<int, string> physical_name;
            physical_name.first = tag;
            physical_name.second = name;
            physical_names.push_back(physical_name);
         }
         continue;
      } 
      else if (line == "$EndPhysicalNames") {
         // Конец физических тегов
         continue;
      }
      else if (line == "$Nodes") {
         // Читаем количество узлов
         getline(file_g2, line);
         nodes_c = stoi(line);
         nodes.resize(nodes_c);
         for (int i = 0; i < nodes_c; ++i) {
            getline(file_g2, line);
            istringstream iss(line);
            int node_id;
            double x, y, z;
            iss >> node_id >> x >> y >> z;
            if (iss.fail()) {
               cerr << "Error reading node data: " << line << endl;
               continue;
            }
            nodes[node_id - 1] = {node_id - 1, x, y, z}; // преобразуем к 0-based индексации
         }
         continue;
      } 
      else if (line == "$EndNodes") {
         // Конец узлов
         continue;
         // break;
      } 
      else if (line == "$Elements") {
         // Читаем количество элементов
         getline(file_g2, line);
         el_c = stoi(line);
         for (int i = 0; i < el_c; ++i) {
            getline(file_g2, line);
            istringstream iss(line);
            EL elem;
            int type, tags_num;
            iss >> elem.number >> type >> tags_num >> elem.physical_tag >> elem.elementary_tag;

            if (iss.fail()) {
               cerr << "Error reading element data: " << line << endl;
               continue;
            }

            // Читаем узлы элемента
            for (int j = 0; j < 8; ++j) {
               int node_id;
               iss >> node_id;
               elem.node_n[j] = node_id - 1; // преобразуем к 0-based
            }

            if(type == 5) el.push_back(elem);
            
            // Если это поверхностный элемент (тип 3), добавляем в faces номер и физический тег
            if (type == 3) {
               for (int i = 0; i < 4; i++) {
                  pair<int, int> face;
                  face.first = elem.node_n[i];
                  face.second =  elem.physical_tag;
                  faces.push_back(face);
               }
            }
         }
         continue;
      } 
      else if (line == "$EndElements") {
         // Конец элементов
         break;
      }
   }

   file_g2.close();

   nodes_c = nodes.size();
   el_c = el.size();
   face_c = faces.size();
}