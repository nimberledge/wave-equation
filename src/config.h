#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>

using std::cerr;
using std::cout;
using std::endl;
using std::string;

struct initial_conditions_t {
  double xmax, ymax, dx, dy, output_dt, t_max, c;
  bool periodic, neumann_bc, benchmark;
  int ic_func_choice;
};

void default_conditions(struct initial_conditions_t &ics) {
  ics.xmax = 10.0;
  ics.ymax = 10.0;
  ics.dx = 0.05;
  ics.dy = 0.05;
  ics.output_dt = 0.5;
  ics.t_max = 30.0;
  ics.c = 1.0;
  ics.periodic = false;
  ics.neumann_bc = true;
  ics.benchmark = false;
  ics.ic_func_choice = 4;
}

void print_initial_conditions(struct initial_conditions_t ics) {
  string func_names[] = {"sin2pix", "sin2pixy", "sinpixy", "sin2piy",
                         "point_disturbance"};
  cout << "xmax: " << ics.xmax << ", ymax: " << ics.ymax << endl;
  cout << "dx: " << ics.dx << ", dy: " << ics.dy
       << ", output_dt: " << ics.output_dt << endl;
  cout << "t_max: " << ics.t_max << ", c (wave speed): " << ics.c << endl;
  cout << "Periodic BCs: " << ics.periodic
       << ", Neumann BCs: " << ics.neumann_bc << endl;
  cout << "Benchmark run: " << ics.benchmark << endl;
  cout << "ic_func_choice: " << ics.ic_func_choice
       << ", ic_func name: " << func_names[ics.ic_func_choice] << endl;
}

void read_conditions_file(string filename, struct initial_conditions_t &ics) {
  default_conditions(ics);
  std::ifstream in_file(filename.c_str());
  if (!in_file.good()) {
    cerr << "Bad initial conditions file, setting values to defaults" << endl;
    return;
  }
  string line;
  int entry_cnt = 0;
  while (std::getline(in_file, line)) {
    // pass
    if (line[0] == '#') {
      continue;
    }
    if (entry_cnt == 0) {
      ics.xmax = std::atof(line.c_str());
    } else if (entry_cnt == 1) {
      ics.ymax = std::atof(line.c_str());
    } else if (entry_cnt == 2) {
      ics.dx = std::atof(line.c_str());
    } else if (entry_cnt == 3) {
      ics.dy = std::atof(line.c_str());
    } else if (entry_cnt == 4) {
      ics.output_dt = std::atof(line.c_str());
    } else if (entry_cnt == 5) {
      ics.t_max = std::atof(line.c_str());
    } else if (entry_cnt == 6) {
      ics.c = std::atof(line.c_str());
    } else if (entry_cnt == 7) {
      ics.periodic = (bool)std::atoi(line.c_str());
    } else if (entry_cnt == 8) {
      ics.neumann_bc = (bool)std::atoi(line.c_str());
    } else if (entry_cnt == 9) {
      ics.benchmark = (bool)std::atoi(line.c_str());
    } else if (entry_cnt == 10) {
      ics.ic_func_choice = std::atoi(line.c_str());
    }
    entry_cnt++;
  }
}

MPI_Datatype get_ic_type(struct initial_conditions_t ics) {
  MPI_Aint addr_start;
  MPI_Get_address(&ics, &addr_start);
  int block_len[11];
  MPI_Datatype types[11];
  MPI_Aint displacements[11];

  MPI_Get_address(&ics.xmax, &displacements[0]);
  MPI_Get_address(&ics.ymax, &displacements[1]);
  MPI_Get_address(&ics.dx, &displacements[2]);
  MPI_Get_address(&ics.dy, &displacements[3]);
  MPI_Get_address(&ics.output_dt, &displacements[4]);
  MPI_Get_address(&ics.t_max, &displacements[5]);
  MPI_Get_address(&ics.c, &displacements[6]);
  MPI_Get_address(&ics.periodic, &displacements[7]);
  MPI_Get_address(&ics.neumann_bc, &displacements[8]);
  MPI_Get_address(&ics.benchmark, &displacements[9]);
  MPI_Get_address(&ics.ic_func_choice, &displacements[10]);
  for (unsigned int i = 0; i < 11; i++) {
    block_len[i] = 1;
    if (i < 7) {
      types[i] = MPI_DOUBLE;
    } else if (i < 10) {
      types[i] = MPI_CXX_BOOL;
    } else {
      types[i] = MPI_INT;
    }
    displacements[i] = displacements[i] - addr_start;
  }
  MPI_Datatype ic_type;
  MPI_Type_create_struct(11, block_len, displacements, types, &ic_type);
  MPI_Type_commit(&ic_type);
  return ic_type;
}

// int main() {
//   string test_filename = "src/conf_file.txt";
//   cout << "Trying to read from: " << test_filename << endl;
//   initial_conditions_t ic = read_conditions_file(test_filename);
//   print_initial_conditions(ic);
//   return 0;
// }