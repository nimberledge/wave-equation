#include "Subdomain.h"
#include "config.h"
#include "initialconditions.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>

using std::cerr;
using std::cout;
using std::endl;
using std::string;

int id, num_proc;
// Define the grid for our rectangular domain
int imax, jmax;
double xmax, ymax;
double dx, dy;
double dt;
double t_max;
double output_dt;
double c;
double (*ic_func)(double, double);
int ic_func_choice;
bool periodic, neumann_bc;
int id_i, id_j;

// Flag to determine if we're running a benchmark.
bool benchmark;

void get_grid_size() {
  // Really strange fix but I think because of the case of 2
  // where the loop finds 2 as a divisor of 2, I have to deal with
  // that case separately. All other primes are fine.
  if (num_proc == 2) {
    imax = 1;
    jmax = 2;
    return;
  }
  imax = 1;
  for (int i = (int)sqrt(num_proc) + 1; i > 1; i--) {
    if (num_proc % i == 0) {
      imax = i;
      break;
    }
  }
  jmax = num_proc / imax;
}

void setup_iteration(struct initial_conditions_t ics) {
  xmax = ics.xmax;
  ymax = ics.ymax;
  dx = ics.dx;
  dy = ics.dy;
  output_dt = ics.output_dt;
  t_max = ics.t_max;
  c = ics.c;
  dt = 0.1 * std::min(dx, dy) / c;
  ic_func_choice = ics.ic_func_choice;
  ic_func = func_choice[ic_func_choice];
  periodic = ics.periodic;
  neumann_bc = ics.neumann_bc;
  benchmark = ics.benchmark;
}

bool internal_bd_square(double x, double y) {
  if (x > 5 && x < 7) {
    if (y > 5 && y < 7) {
      return true;
    }
  }
  return false;
}

bool internal_bd_circle(double x, double y) {
  double dist = pow(x - 7, 2) + pow(y - 2, 2);
  if (dist < 1) {
    return true;
  }
  return false;
}

bool test_compound_internal_bd(double x, double y) {
  return internal_bd_square(x, y) || internal_bd_circle(x, y);
}

bool no_internal_bd(double x, double y) { return false; }

bool tall_thin_internal_bd(double x, double y) {
  if (5 < x && x < 6) {
    if (0.1 < y && y < 7.5) {
      return true;
    }
  }
  return false;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  struct initial_conditions_t ics;
  MPI_Datatype ic_type = get_ic_type(ics);
  if (id == 0) {
    if (argc > 1) {
      char *file_char_star = argv[1];
      string ic_filename(file_char_star);
      read_conditions_file(ic_filename, ics);
    } else {
      default_conditions(ics);
    }
    // cout << "Trying to broadcast ic_type" << endl;
    MPI_Bcast(&ics, 1, ic_type, 0, MPI_COMM_WORLD);
    print_initial_conditions(ics);
  } else {
    // cout << id << " trying to receive the ic_type" << endl;
    MPI_Bcast(&ics, 1, ic_type, 0, MPI_COMM_WORLD);
  }

  get_grid_size();
  id_i = id / jmax;
  id_j = id % jmax;
  setup_iteration(ics);
  if (id == 0) {
    cout << "iteration dt: " << dt << endl;
    cout << "num_proc: " << num_proc << endl;
  }

  Subdomain sd(id_i, id_j, imax, jmax, periodic, xmax, ymax, dx, dy, c,
               neumann_bc, &tall_thin_internal_bd);
  sd.set_initial_condition(ic_func);
  char fname_buffer[100];
  double t_current = 0.0;
  double t_interval = 0.0;
  int iter = 0;

  // Write the initial condition - only if not benchmarking
  if (!benchmark) {
    sprintf(fname_buffer, "data/iter_%d_thread_%d.txt", iter, id);
    string it_file(fname_buffer);
    if (!sd.write_grid(it_file)) {
      cerr << "id: " << id << " failed to write file: " << it_file << endl;
    }
  }
  double start_time = MPI_Wtime();
  while (t_current < t_max) {
    if (!benchmark) {
      if (t_interval >= output_dt) {
        sprintf(fname_buffer, "data/iter_%d_thread_%d.txt", ++iter, id);
        string it_file(fname_buffer);
        if (!sd.write_grid(it_file)) {
          cerr << "id: " << id << " failed to write file: " << it_file << endl;
        }
        t_interval = 0.0;
      }
    }
    sd.do_update(dt);
    t_current += dt;
    t_interval += dt;
  }
  double taken = MPI_Wtime() - start_time;
  MPI_Barrier(MPI_COMM_WORLD);
  if (id == 0) {
    cout << "Program completed in " << taken << " seconds." << endl;
  }
  // Just write out the last iteration's grid, in case we missed it
  if (t_interval > 0 && !benchmark) {
    sprintf(fname_buffer, "data/iter_%d_thread_%d.txt", ++iter, id);
    string it_file(fname_buffer);
    if (!sd.write_grid(it_file)) {
      cerr << "id: " << id << " failed to write file: " << it_file << endl;
    }
  }
  MPI_Finalize();
  return 0;
}
