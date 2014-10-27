#include <stdio.h>
#include <stdlib.h>
#include "duffing.h"
#include "euler.h"

// Print a line
//    time x[0] x[1] ...
// to standard out
void PrintState(double n, double t, const double *x) {
  printf("%15.8f", t);
  for (int i = 0; i < n; ++i) {
    printf("%15.8f", x[i]);
  }
  printf("\n");
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    printf("USAGE: %s <timestep> <numsteps>\n", argv[0]);
    exit(1);
  }
  const double dt = atof(argv[1]);
  const int nsteps = atoi(argv[2]);

  const double delta = 0.2;
  const double gamma = 0.3;
  const double omega = 1.0;
  Model *model = new Duffing(delta, gamma, omega);

  // Initial condition
  const int dimen = model->dimen();   // number of states
  double x[dimen];
  for (int i = 0; i < dimen; ++i) {
    x[i] = 0;
  }

  // Set up integration scheme
  Integrator *integrator = new Euler(dt, *model);

  double t = 0;
  PrintState(dimen, t, x);
  for (int i = 0; i < nsteps; ++i) {
    integrator->Step(t, x);
    t = (i+1) * dt;
    PrintState(dimen, t, x);
  }

  delete integrator;
  delete model;
  return 0;
}
