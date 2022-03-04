#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

class Subdomain {
public:
  // Use neumann_bc = 0 for Homogenous Dirichlet BCs, neumann_bc = 1 for
  // homogenous Neumann conditions
  Subdomain(int id_i, int id_j, int imax, int jmax, bool periodic, double xmax,
            double ymax, double dx, double dy, double c, bool neumann_bc,
            bool (*internal_bd_func)(double, double));
  ~Subdomain();
  double *old_grid, *current_grid, *new_grid;

  // Pass time-step here
  void do_update(double dt);

  MPI_Datatype Left_Type, Right_Type, Top_Type, Bottom_Type;
  // The following variable is just used to make communication
  // loops cleaner. Order is (L, R, T, B)
  MPI_Datatype MPI_TypeList[4];
  void create_MPI_Types();
  void set_nbrs();
  // Pass a function f(x, y) as the initial condition
  void set_initial_condition(double (*func_ptr)(double, double));
  bool write_grid(string filename);

private:
  int imax, jmax, id_i, id_j;
  int rows, cols;
  // Boundary condition code. 0 = Dirichlet, 1 = Neumann
  bool neumann_bc;
  // periodic = true skips either of the above BCs
  bool periodic;
  // Bools to determine if we have to apply BCs
  bool bd_left, bd_right, bd_top, bd_bottom;
  double dx, dy, xmax, ymax, xstart, ystart, c;
  // List of neighbours, -1 if neighbour doesn't exist
  // order is (L, R, T, B)
  int nbrs[4];
  int num_nbrs;
  bool (*internal_bd)(double, double);
};