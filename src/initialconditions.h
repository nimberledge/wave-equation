#define _USE_MATH_DEFINES
#include <cmath>

extern double xmax, ymax;
double sin2pix(double x, double y) { return sin(2 * M_PI * x / xmax); }

double sin2pixy(double x, double y) { return sin(2 * M_PI * x * y); }

double sinpixy(double x, double y) { return sin(M_PI * x * y); }

double sin2piy(double x, double y) { return sin(2 * M_PI * y / ymax); }

double point_disturbance(double x, double y) {
  double dist = (x - 3.) * (x - 3.) + (y - 3.) * (y - 3.);
  double r_splash = 1.0;
  if (dist < r_splash) {
    return cos(dist / r_splash * M_PI) + 1;
  }
  return 0.0;
}

double (*func_choice[])(double, double) = {&sin2pix, &sin2pixy, &sinpixy,
                                           &sin2piy, &point_disturbance};