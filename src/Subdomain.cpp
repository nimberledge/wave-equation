#include "Subdomain.h"

Subdomain::Subdomain(int id_i, int id_j, int imax, int jmax, bool periodic,
                     double xmax, double ymax, double dx, double dy, double c,
                     bool neumann_bc, bool (*internal_bd_func)(double, double))
    : id_i(id_i), id_j(id_j), imax(imax), jmax(jmax), periodic(periodic),
      xmax(xmax), ymax(ymax), dx(dx), dy(dy), c(c), neumann_bc(neumann_bc),
      internal_bd(internal_bd_func) {
  double width = xmax / jmax;
  double height = ymax / imax;
  // Add one to ensure we have a point on the boundary
  cols = (int)(width / dx) + 1;
  rows = (int)(height / dy) + 1;
  xstart = id_j * cols * dx;
  ystart = id_i * rows * dy;
  old_grid = new double[rows * cols];
  current_grid = new double[rows * cols];
  new_grid = new double[rows * cols];
  bd_bottom = (id_i == 0);
  bd_top = (id_i == imax - 1);
  bd_left = (id_j == 0);
  bd_right = (id_j == jmax - 1);
  num_nbrs = 0;
  this->set_nbrs();
  this->create_MPI_Types();
}

Subdomain::~Subdomain() {
  delete[] old_grid;
  delete[] current_grid;
  delete[] new_grid;
}

void Subdomain::set_nbrs() {
  // Left neighbour - only -1 if id_j is 0 and no periodic BC
  // the neighbour is the process that would be on the same row
  // and one column left (id_i * jmax) + id_j - 1
  // TODO: Ensure we do not communicate with ourselves.
  // Now doing that TODO - if bd_left and bd_right, then we
  // don't have to send stuff across the boundary
  // same condition for RHS boundary
  if (id_j == 0) {
    if (periodic && !(bd_left && bd_right)) {
      int nbr_id = (id_i * jmax) + jmax - 1;
      nbrs[0] = nbr_id;
      num_nbrs++;
    } else {
      nbrs[0] = -1;
    }
  } else {
    int nbr_id = (id_i * jmax) + id_j - 1;
    nbrs[0] = nbr_id;
    num_nbrs++;
  }
  // Right neighbour - only -1 if id_j = jmax-1 and no periodic BC
  // the neighbour is the process that would be on the same row
  // and one column right
  if (id_j == jmax - 1) {
    if (periodic && !(bd_left && bd_right)) {
      int nbr_id = (id_i * jmax);
      nbrs[1] = nbr_id;
      num_nbrs++;
    } else {
      nbrs[1] = -1;
    }
  } else {
    int nbr_id = id_i * jmax + id_j + 1;
    nbrs[1] = nbr_id;
    num_nbrs++;
  }
  // Top neighbour - only -1 if id_i = imax-1 and no periodic BC
  // the neighbour is the process that would be one row up and
  // on the same column
  // We only communicate with ourselves if bd_top and bd_bottom
  // are both true and it's periodic
  if (id_i == imax - 1) {
    if (periodic && !(bd_top && bd_bottom)) {
      int nbr_id = id_j;
      nbrs[2] = nbr_id;
      num_nbrs++;
    } else {
      nbrs[2] = -1;
    }
  } else {
    int nbr_id = (id_i + 1) * jmax + id_j;
    nbrs[2] = nbr_id;
    num_nbrs++;
  }
  // Bottom neighbour - only -1 if id_i = 0 and no periodic BC
  // the neighbour is the process that would be one row down and
  // on the same column
  if (id_i == 0) {
    if (periodic && !(bd_top && bd_bottom)) {
      int nbr_id = (imax - 1) * jmax + id_j;
      nbrs[3] = nbr_id;
      num_nbrs++;
    } else {
      nbrs[3] = -1;
    }
  } else {
    int nbr_id = (id_i - 1) * jmax + id_j;
    nbrs[3] = nbr_id;
    num_nbrs++;
  }
}

void Subdomain::create_MPI_Types() {
  MPI_Aint addr_start;
  int block_len[rows];
  MPI_Datatype type_arr[rows];
  MPI_Aint displacements[rows];

  MPI_Get_address(this, &addr_start);
  // Left Type
  for (unsigned int i = 0; i < rows; i++) {
    type_arr[i] = MPI_DOUBLE;
    block_len[i] = 1;
    MPI_Get_address(&current_grid[i * cols], &displacements[i]);
    displacements[i] = displacements[i] - addr_start;
  }
  MPI_Type_create_struct(rows, block_len, displacements, type_arr, &Left_Type);
  MPI_Type_commit(&Left_Type);
  // Right Type
  for (unsigned int i = 0; i < rows; i++) {
    MPI_Get_address(&current_grid[i * cols + cols - 1], &displacements[i]);
    displacements[i] = displacements[i] - addr_start;
  }
  MPI_Type_create_struct(rows, block_len, displacements, type_arr, &Right_Type);
  MPI_Type_commit(&Right_Type);

  // For top and bottom types we can just specify the starting address and block
  // size, since MPI will just read that contiguous memory and put it in a type
  // nicely
  // Top Type
  block_len[0] = cols;
  MPI_Get_address(&current_grid[(rows - 1) * (cols)], &displacements[0]);
  displacements[0] = displacements[0] - addr_start;
  MPI_Type_create_struct(1, block_len, displacements, type_arr, &Top_Type);
  MPI_Type_commit(&Top_Type);

  // Bottom Type
  MPI_Get_address(&current_grid[0], &displacements[0]);
  displacements[0] = displacements[0] - addr_start;
  MPI_Type_create_struct(1, block_len, displacements, type_arr, &Bottom_Type);
  MPI_Type_commit(&Bottom_Type);

  MPI_TypeList[0] = Left_Type;
  MPI_TypeList[1] = Right_Type;
  MPI_TypeList[2] = Top_Type;
  MPI_TypeList[3] = Bottom_Type;
}

void Subdomain::set_initial_condition(double (*func_ptr)(double, double)) {
  for (unsigned int i = 0; i < rows; i++) {
    for (unsigned int j = 0; j < cols; j++) {
      double x = xstart + j * dx;
      double y = ystart + i * dy;
      if (internal_bd(x, y)) {
        old_grid[i * cols + j] = 0.0;
        current_grid[i * cols + j] = 0.0;
        continue;
      }
      double fval = func_ptr(x, y);
      old_grid[i * cols + j] = fval;
      current_grid[i * cols + j] = fval;
    }
  }
  // If we have a periodic boundary, we don't have to do anything else here
  if (periodic) {
    return;
  }
  // Check for Dirichlet BCs
  if (!neumann_bc) {
    if (bd_left) {
      for (unsigned int i = 0; i < rows; i++) {
        old_grid[i * cols] = 0.0;
        current_grid[i * cols] = 0.0;
      }
    }
    if (bd_right) {
      for (unsigned int i = 0; i < rows; i++) {
        old_grid[i * cols + cols - 1] = 0.0;
        current_grid[i * cols + cols - 1] = 0.0;
      }
    }
    if (bd_top) {
      for (unsigned int j = 0; j < cols; j++) {
        old_grid[(rows - 1) * cols + j] = 0.0;
        current_grid[(rows - 1) * cols + j] = 0.0;
      }
    }
    if (bd_bottom) {
      for (unsigned int j = 0; j < cols; j++) {
        old_grid[j] = 0.0;
        current_grid[j] = 0.0;
      }
    }
    return;
  }
  // If we're here, we have Neumann BCs
  // Left boundary
  if (bd_left) {
    for (unsigned int i = 0; i < rows; i++) {
      double fval = old_grid[i * cols + 1];
      old_grid[i * cols] = fval;
      current_grid[i * cols] = fval;
    }
  }
  // Bottom boundary
  if (bd_bottom) {
    for (unsigned int j = 0; j < cols; j++) {
      double fval = old_grid[cols + j];
      old_grid[j] = fval;
      current_grid[j] = fval;
    }
  }
  // Right boundary
  if (bd_right) {
    for (unsigned int i = 0; i < rows; i++) {
      double fval = old_grid[i * cols + cols - 2];
      old_grid[i * cols + cols - 1] = fval;
      current_grid[i * cols + cols - 1] = fval;
    }
  }
  // Top boundary
  if (bd_top) {
    for (unsigned int j = 0; j < cols; j++) {
      double fval = old_grid[(rows - 2) * cols + j];
      old_grid[(rows - 1) * cols + j] = fval;
      current_grid[(rows - 1) * cols + j] = fval;
    }
  }
}

bool Subdomain::write_grid(string filename) {
  std::ofstream out_file(filename);
  if (!out_file.good()) {
    cerr << "Bad filename, aborting write" << endl;
    return false;
  }
  for (unsigned int i = 0; i < rows; i++) {
    for (unsigned int j = 0; j < cols; j++) {
      out_file << current_grid[i * cols + j] << "\t";
    }
    out_file << "\n";
  }
  out_file.close();
  return true;
}

void Subdomain::do_update(double dt) {
  // First do comms
  double left_bdvals[rows], right_bdvals[rows], top_bdvals[cols],
      bottom_bdvals[cols];
  double *bd_data[] = {left_bdvals, right_bdvals, top_bdvals, bottom_bdvals};
  MPI_Request requests[2 * num_nbrs];
  int req_cnt = 0;
  for (int i = 0; i < 4; i++) {
    if (nbrs[i] >= 0) {
      // MPI send
      // Send L, R, T, B as    0, 1, 2, 3
      // Receive L, R, T, B as 1, 0, 3, 2
      int send_tag_num = i;
      int recv_tag_num;
      if (send_tag_num == 0) {
        recv_tag_num = 1;
      } else if (send_tag_num == 1) {
        recv_tag_num = 0;
      } else if (send_tag_num == 2) {
        recv_tag_num = 3;
      } else {
        recv_tag_num = 2;
      }

      MPI_Isend(this, 1, this->MPI_TypeList[i], nbrs[i], send_tag_num,
                MPI_COMM_WORLD, &requests[req_cnt++]);
      int recv_cnt = rows;
      if (i >= 2) {
        recv_cnt = cols;
      }
      MPI_Irecv(bd_data[i], recv_cnt, MPI_DOUBLE, nbrs[i], recv_tag_num,
                MPI_COMM_WORLD, &requests[req_cnt++]);
    }
  }
  // Precompute some values to simplify calculations
  double cdt2 = pow(dt * c, 2);
  double dx2 = pow(dx, 2);
  double dy2 = pow(dy, 2);
  // Now do the inner iteration
  for (unsigned int i = 1; i < rows - 1; i++) {
    for (unsigned int j = 1; j < cols - 1; j++) {
      double x = xstart + j * dx;
      double y = ystart + i * dy;
      if (internal_bd(x, y)) {
        new_grid[i * cols + j] = 0.0;
        continue;
      }
      new_grid[i * cols + j] = cdt2 * (((current_grid[(i + 1) * cols + j] +
                                         current_grid[(i - 1) * cols + j] -
                                         2 * current_grid[i * cols + j]) /
                                        dy2) +
                                       ((current_grid[i * cols + j + 1] +
                                         current_grid[i * cols + j - 1] -
                                         2 * current_grid[i * cols + j]) /
                                        dx2)) +
                               2 * current_grid[i * cols + j] -
                               old_grid[i * cols + j];
    }
  }
  // Now that we are done with the inner loop, think about the comms
  // so that we can do the boundaries
  if (num_nbrs > 0) {
    MPI_Waitall(req_cnt, requests, MPI_STATUSES_IGNORE);
  }
  // Deal with subdomain boundaries, only loop if we have that neighbour,
  // we'll deal with BCs later
  // Left boundary
  if (nbrs[0] >= 0) {
    int col = 0;
    for (unsigned int i = 1; i < rows - 1; i++) {
      double x = xstart + col * dx;
      double y = ystart + i * dy;
      if (internal_bd(x, y)) {
        new_grid[i * cols + col] = 0.0;
        continue;
      }
      new_grid[i * cols + col] =
          cdt2 * (((current_grid[(i + 1) * cols + col] +
                    current_grid[(i - 1) * cols + col] -
                    2 * current_grid[i * cols + col]) /
                   dy2) +
                  ((current_grid[i * cols + col + 1] + left_bdvals[i] -
                    2 * current_grid[i * cols + col]) /
                   dx2)) +
          2 * current_grid[i * cols + col] - old_grid[i * cols + col];
    }
  }
  // Right boundary
  if (nbrs[1] >= 0) {
    int col = cols - 1;
    for (unsigned int i = 1; i < rows - 1; i++) {
      double x = xstart + col * dx;
      double y = ystart + i * dy;
      if (internal_bd(x, y)) {
        new_grid[i * cols + col] = 0.0;
        continue;
      }
      new_grid[i * cols + col] =
          cdt2 * (((current_grid[(i + 1) * cols + col] +
                    current_grid[(i - 1) * cols + col] -
                    2 * current_grid[i * cols + col]) /
                   dy2) +
                  ((current_grid[i * cols + col - 1] + right_bdvals[i] -
                    2 * current_grid[i * cols + col]) /
                   dx2)) +
          2 * current_grid[i * cols + col] - old_grid[i * cols + col];
    }
  }
  // Top boundary
  if (nbrs[2] >= 0) {
    int row = rows - 1;
    for (unsigned int j = 1; j < cols - 1; j++) {
      double x = xstart + j * dx;
      double y = ystart + row * dy;
      if (internal_bd(x, y)) {
        new_grid[row * cols + j] = 0.0;
        continue;
      }
      new_grid[row * cols + j] =
          cdt2 * (((current_grid[(row - 1) * cols + j] + top_bdvals[j] -
                    2 * current_grid[row * cols + j]) /
                   dy2) +
                  ((current_grid[row * cols + j - 1] +
                    current_grid[row * cols + j + 1] -
                    2 * current_grid[row * cols + j]) /
                   dx2)) +
          2 * current_grid[row * cols + j] - old_grid[row * cols + j];
    }
  }
  // Bottom boundary
  if (nbrs[3] >= 0) {
    int row = 0;
    for (unsigned int j = 1; j < cols - 1; j++) {
      double x = xstart + j * dx;
      double y = ystart + row * dy;
      if (internal_bd(x, y)) {
        new_grid[row * cols + j] = 0.0;
        continue;
      }
      new_grid[row * cols + j] =
          cdt2 * (((current_grid[(row + 1) * cols + j] + bottom_bdvals[j] -
                    2 * current_grid[row * cols + j]) /
                   dy2) +
                  ((current_grid[row * cols + j - 1] +
                    current_grid[row * cols + j + 1] -
                    2 * current_grid[row * cols + j]) /
                   dx2)) +
          2 * current_grid[row * cols + j] - old_grid[row * cols + j];
    }
  }
  // Now for the corners since that requires special treatment
  // Top left corner
  if (nbrs[0] >= 0 && nbrs[2] >= 0) {
    int row = rows - 1;
    int col = 0;
    double x = xstart + col * dx;
    double y = ystart + row * dy;
    if (internal_bd(x, y)) {
      new_grid[row * cols + col] = 0.0;
    } else {
      new_grid[row * cols + col] =
          cdt2 * (((top_bdvals[col] + current_grid[(row - 1) * cols + col] -
                    2 * current_grid[row * cols + col]) /
                   dy2) +
                  ((left_bdvals[row] + current_grid[row * cols + col + 1] -
                    2 * current_grid[row * cols + col]) /
                   dx2)) +
          2 * current_grid[row * cols + col] - old_grid[row * cols + col];
    }
  }
  // Top right corner
  if (nbrs[1] >= 0 && nbrs[2] >= 0) {
    int row = rows - 1;
    int col = cols - 1;
    double x = xstart + col * dx;
    double y = ystart + row * dy;
    if (internal_bd(x, y)) {
      new_grid[row * cols + col] = 0.0;
    } else {
      new_grid[row * cols + col] =
          cdt2 * (((top_bdvals[col] + current_grid[(row - 1) * cols + col] -
                    2 * current_grid[row * cols + col]) /
                   dy2) +
                  ((right_bdvals[row] + current_grid[row * cols + col - 1] -
                    2 * current_grid[row * cols + col]) /
                   dx2)) +
          2 * current_grid[row * cols + col] - old_grid[row * cols + col];
    }
  }
  // Bottom left corner
  if (nbrs[0] >= 0 && nbrs[3] >= 0) {
    int row = 0;
    int col = 0;
    double x = xstart + col * dx;
    double y = ystart + row * dy;
    if (internal_bd(x, y)) {
      new_grid[row * cols + col] = 0.0;
    } else {
      new_grid[row * cols + col] =
          cdt2 * (((bottom_bdvals[col] + current_grid[(row + 1) * cols + col] -
                    2 * current_grid[row * cols + col]) /
                   dy2) +
                  ((left_bdvals[row] + current_grid[row * cols + col + 1] -
                    2 * current_grid[row * cols + col]) /
                   dx2)) +
          2 * current_grid[row * cols + col] - old_grid[row * cols + col];
    }
  }
  // Bottom right corner
  if (nbrs[1] >= 0 && nbrs[3] >= 0) {
    int row = 0;
    int col = cols - 1;
    double x = xstart + col * dx;
    double y = ystart + row * dy;
    if (internal_bd(x, y)) {
      new_grid[row * cols + col] = 0.0;
    } else {
      new_grid[row * cols + col] =
          cdt2 * (((bottom_bdvals[col] + current_grid[(row + 1) * cols + col] -
                    2 * current_grid[row * cols + col]) /
                   dy2) +
                  ((right_bdvals[row] + current_grid[row * cols + col - 1] -
                    2 * current_grid[row * cols + col]) /
                   dx2)) +
          2 * current_grid[row * cols + col] - old_grid[row * cols + col];
    }
  }

  // If we have Dirichlet BCs, then we're done here
  if (!periodic && !neumann_bc) {
    std::swap(current_grid, new_grid);
    std::swap(old_grid, new_grid);
    return;
  }
  // If we have a periodic boundary, let's just check if we have to send
  // comms to ourselves.
  if (periodic) {
    // This split-up is nice because we can be assured we have
    // bd_vals on the remaining sides. Unless it's one process.
    // Damn.
    if (bd_top && bd_bottom) {
      for (unsigned int j = 1; j < cols - 1; j++) {
        // First do the top one, just replace top_bdvals with
        // the value of the bottom row at column j
        int row = rows - 1;
        new_grid[row * cols + j] =
            cdt2 * (((current_grid[(row - 1) * cols + j] +
                      current_grid[0 * cols + j] -
                      2 * current_grid[row * cols + j]) /
                     dy2) +
                    ((current_grid[row * cols + j - 1] +
                      current_grid[row * cols + j + 1] -
                      2 * current_grid[row * cols + j]) /
                     dx2)) +
            2 * current_grid[row * cols + j] - old_grid[row * cols + j];
        // Now the bottom one
        row = 0;
        new_grid[row * cols + j] =
            cdt2 * (((current_grid[(row + 1) * cols + j] +
                      current_grid[(rows - 1) * cols + j] -
                      2 * current_grid[row * cols + j]) /
                     dy2) +
                    ((current_grid[row * cols + j - 1] +
                      current_grid[row * cols + j + 1] -
                      2 * current_grid[row * cols + j]) /
                     dx2)) +
            2 * current_grid[row * cols + j] - old_grid[row * cols + j];
      }
      // Now do the corners
      // Case 1, there's only one process, and we're in charge of all 4
      // boundaries
      if (bd_left && bd_right) {
        // Bottom left corner
        int row = 0;
        int col = 0;
        new_grid[row * cols + col] =
            cdt2 * (((current_grid[cols + col] +
                      current_grid[(rows - 1) * cols + col] -
                      2 * current_grid[row * cols + col]) /
                     dy2) +
                    ((current_grid[row * cols + (cols - 1)] +
                      current_grid[row * cols + 1] -
                      2 * current_grid[row * cols + col]) /
                     dx2)) +
            2 * current_grid[row * cols + col] - old_grid[row * cols + col];
        // Bottom right corner
        col = cols - 1;
        new_grid[row * cols + col] =
            cdt2 * (((current_grid[cols + col] +
                      current_grid[(rows - 1) * cols + col] -
                      2 * current_grid[row * cols + col]) /
                     dy2) +
                    ((current_grid[row * cols + (col - 1)] +
                      current_grid[row * cols + 0] -
                      2 * current_grid[row * cols + col]) /
                     dx2)) +
            2 * current_grid[row * cols + col] - old_grid[row * cols + col];
        // Top right corner
        row = rows - 1;
        new_grid[row * cols + col] =
            cdt2 * (((current_grid[0 * cols + col] +
                      current_grid[(row - 1) * cols + col] -
                      2 * current_grid[row * cols + col]) /
                     dy2) +
                    ((current_grid[row * cols + (col - 1)] +
                      current_grid[row * cols + 0] -
                      2 * current_grid[row * cols + col]) /
                     dx2)) +
            2 * current_grid[row * cols + col] - old_grid[row * cols + col];
        // Top left corner
        col = 0;
        new_grid[row * cols + col] =
            cdt2 * (((current_grid[0 * cols + col] +
                      current_grid[(row - 1) * cols + col] -
                      2 * current_grid[row * cols + col]) /
                     dy2) +
                    ((current_grid[row * cols + (cols - 1)] +
                      current_grid[row * cols + col + 1] -
                      2 * current_grid[row * cols + col]) /
                     dx2)) +
            2 * current_grid[row * cols + col] - old_grid[row * cols + col];
      } else {
        // Case 2: We have left_bdvals and right_bdvals
        int row = 0;
        int col = 0;
        // Bottom left corner
        new_grid[row * cols + col] =
            cdt2 * (((current_grid[cols + col] +
                      current_grid[(rows - 1) * cols + col] -
                      2 * current_grid[row * cols + col]) /
                     dy2) +
                    ((left_bdvals[row] + current_grid[row * cols + col + 1] -
                      2 * current_grid[row * cols + col]) /
                     dx2)) +
            2 * current_grid[row * cols + col] - old_grid[row * cols + col];
        // Bottom right corner
        col = cols - 1;
        new_grid[row * cols + col] =
            cdt2 * (((current_grid[cols + col] +
                      current_grid[(rows - 1) * cols + col] -
                      2 * current_grid[row * cols + col]) /
                     dy2) +
                    ((right_bdvals[row] + current_grid[row * cols + col - 1] -
                      2 * current_grid[row * cols + col]) /
                     dx2)) +
            2 * current_grid[row * cols + col] - old_grid[row * cols + col];
        // Top right corner
        row = rows - 1;
        new_grid[row * cols + col] =
            cdt2 * (((current_grid[0 * cols + col] +
                      current_grid[(row - 1) * cols + col] -
                      2 * current_grid[row * cols + col]) /
                     dy2) +
                    ((current_grid[row * cols + col - 1] + right_bdvals[row] -
                      2 * current_grid[row * cols + col]) /
                     dx2)) +
            2 * current_grid[row * cols + col] - old_grid[row * cols + col];
        // Top left corner
        col = 0;
        new_grid[row * cols + col] =
            cdt2 * (((current_grid[0 * cols + col] +
                      current_grid[(row - 1) * cols + col] -
                      2 * current_grid[row * cols + col]) /
                     dy2) +
                    ((current_grid[row * cols + col + 1] + left_bdvals[row] -
                      2 * current_grid[row * cols + col]) /
                     dx2)) +
            2 * current_grid[row * cols + col] - old_grid[row * cols + col];
      }
    }
    // Now check if left and right need to communicate
    // based on how the code is structured, we'll only ever get this if
    // we have only one process. The top case is faaar more important.
    // but for completeness, we do this. If we have bd_left && bd_right, we are
    // actually assured to have bd_top && bd_bottom
    // This actually means we can also skip the corners
    if (bd_left && bd_right) {
      for (unsigned int i = 1; i < rows - 1; i++) {
        // Do left boundary first
        int col = 0;
        new_grid[i * cols + col] =
            cdt2 * (((current_grid[(i + 1) * cols + col] +
                      current_grid[(i - 1) * cols + col] -
                      2 * current_grid[i * cols + col]) /
                     dy2) +
                    ((current_grid[i * cols + col + 1] +
                      current_grid[i * cols + cols - 1] -
                      2 * current_grid[i * cols + col]) /
                     dx2)) +
            2 * current_grid[i * cols + col] - old_grid[i * cols + col];
        // Now the right boundary
        col = cols - 1;
        new_grid[i * cols + col] =
            cdt2 * (((current_grid[(i + 1) * cols + col] +
                      current_grid[(i - 1) * cols + col] -
                      2 * current_grid[i * cols + col]) /
                     dy2) +
                    ((current_grid[i * cols + col - 1] +
                      current_grid[i * cols + 0] -
                      2 * current_grid[i * cols + col]) /
                     dx2)) +
            2 * current_grid[i * cols + col] - old_grid[i * cols + col];
      }
    }
    std::swap(current_grid, new_grid);
    std::swap(old_grid, new_grid);
    return;
  }

  // If we're here it's a non-periodic Neumann BC how great
  // Left domain boundary
  if (bd_left) {
    for (unsigned int i = 0; i < rows; i++) {
      new_grid[i * cols] = new_grid[i * cols + 1];
    }
  }
  // Right domain boundary
  if (bd_right) {
    for (unsigned int i = 0; i < rows; i++) {
      new_grid[i * cols + cols - 1] = new_grid[i * cols + cols - 2];
    }
  }
  // Top domain boundary
  if (bd_top) {
    for (unsigned int j = 0; j < cols; j++) {
      new_grid[(rows - 1) * cols + j] = new_grid[(rows - 2) * cols + j];
    }
  }
  // Bottom domain boundary
  if (bd_bottom) {
    for (unsigned int j = 0; j < cols; j++) {
      new_grid[j] = new_grid[cols + j];
    }
  }

  // Swap the pointers now
  std::swap(current_grid, new_grid);
  std::swap(old_grid, new_grid);
}