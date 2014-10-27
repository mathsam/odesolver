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
#include "math.h"

int linear_exact(double t,     double beta, 
                 double gamma, double omega,
                 const double *x_initial, double *x_out){
    static double A = (1-omega*omega)*gamma
                     /((1-omega*omega)*(1-omega*omega)+4*beta*beta*omega*omega);
    static double B = 2*beta*omega*gamma
                     /((1-omega*omega)*(1-omega*omega)+4*beta*beta*omega*omega);
    if(beta*beta>=1) {
        printf("Wrong parameters for linear oscillator: beta should be smaller than 1");
        exit(1);
    }
    static double omega_d = sqrt(1-beta*beta); 
    double c1 = x_initial[0] - A; //x[0] -> x; x[1] -> \dot x 
    double c2 = (x_initial[1] + beta*(x_initial[0]-A) - omega*B)/omega_d;
    x_out[0] = exp(-beta*t) * (c1*cos(omega_d*t) + c2*sin(omega_d*t))
               + A*cos(omega*t) + B*sin(omega*t);
    x_out[1] = exp(-beta*t) * (-beta * (c1*cos(omega_d*t) + c2*sin(omega_d*t))
                               + omega_d*(-c1*sin(omega_d*t) + c2*cos(omega_d*t)))
               - omega*(A*sin(omega*t) - B*cos(omega*t));
    return 0; 
}

/* calculate (x1[0] - x2[0])^2 + (x1[1] - x2[1])^2 + ... */
double error_square(const double *x1, const double *x2, int dimen){
  double norm = 0;
  for(int i=0;i<dimen;++i)
    norm += (x1[i] - x2[i])*(x1[i] - x2[i]);
  return norm;
}

int main(int argc, char *argv[]){
  if (argc != 4) {
    printf("USAGE: %s <integrator> <timestep> <numsteps>\n", argv[0]);
    exit(1);
  }
  std::string integrator  (argv[1]);
  const double dt = atof  (argv[2]);
  const int nsteps= atoi  (argv[3]);

  Model  * model;
  double * x;
  const double beta  = 0.1;
  const double gamma = 1.0;
  const double omega = 0.9;
  model = new LinearOscillator (beta,gamma,omega);
  x     = new double[2]();

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
  double error = 0;
  double x_exact[dimen];
  double x_initial[dimen];
  for(int i=0; i<dimen; ++i) x_initial[i] = 0;

  for(int i=0; i<nsteps; ++i){
    solver->Step(t,x);
    t = (i+1) * dt;

    linear_exact(t, beta, gamma, omega, x_initial, x_exact);
    error += error_square(x,x_exact,dimen);
  }
  
  double error_two_norm;
  error_two_norm = sqrt(error*dt);
  printf("%f\n",error_two_norm);

  delete solver;
  delete model;
  delete [] x;
  return 0;
}
