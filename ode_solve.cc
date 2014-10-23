#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "model.h"
#include "integrator.h"
#include "duffing.h"
#include "lorenz.h"
#include "linear-oscillator.h"
#include "euler.h"
#include "runge-kutta.h"
#include "adams-bashforth.h"

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


int main(int argc, char *argv[]){
  if (argc != 5) {
    printf("USAGE: %s <equation> <integrator> <timestep> <numsteps>\n", argv[0]);
    exit(1);
  }
  std::string equation    (argv[1]);
  std::string integrator  (argv[2]);
  const double dt = atof  (argv[3]);
  const int nsteps= atoi  (argv[4]);

  Model  * model;
  double * x;
  if(equation == "duffing"){
    const double delta = 0.2;
    const double gamma = 0.3;
    const double omega = 1.0;
    model = new Duffing (delta,gamma,omega);
    x     = new double[2] ();
  }
  else if(equation == "lorenz"){
    const double sigma = 10.0; 
    const double rho   = 28.0;
    const double beta  = 8.0/3.0;
    model     = new Lorenz (sigma,rho,beta);
    x         = new double[3] ();
    x[1] = 0.01;
  }
  else if(equation == "linear"){
    const double beta  = 0.1;
    const double gamma = 1.0;
    const double omega = 0.9;
    model = new LinearOscillator (beta,gamma,omega);
    x     = new double[2]();
  }
  else{
    printf("equation must be one of duffing, lorenz and linear\n");
    exit(1);
  }

  Integrator * solver;
  if(integrator == "euler"){
    solver = new Euler (dt, *model);
  }
  else if(integrator == "rk4"){
    solver = new Rk4 (dt, *model); 
  }
  else if(integrator == "ab2"){
    solver = new AdamsBashforth (dt, *model);
  }
  else{
    printf("integrator must be one of euler, rk4 and ab2\n");
    exit(1);
  }

  const int dimen = model->dimen();
  double t = 0;
  PrintState(dimen, t, x); 
  for(int i=0; i<nsteps; ++i){
    solver->Step(t,x);
    t = (i+1) * dt;
    PrintState(dimen,t,x);
  }

  delete solver;
  delete model;
  delete [] x;
  return 0;
}
