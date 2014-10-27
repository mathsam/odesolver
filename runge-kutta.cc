#include "runge-kutta.h"
#include "model.h"

Rk4::Rk4(double dt, const Model &model)
    : dimen_(model.dimen()),
      dt_(dt),
      model_(model) {}

Rk4::~Rk4(){}

/* 4th order Runge-Kutta                           */
/* for an ODE dx/dt = f(x,t), it is integrated as  */
/* x_j+1 = x_j + dt/6*(k1+2*k2+2*k3+k4), where     */
/* k1 = f(x_j,t)                                   */
/* k2 = f(x_j+0.5*dt*k1, t+0.5*dt)                 */
/* k3 = f(x_j+0.5*dt*k2, t+0.5*dt)                 */
/* k4 = f(x_j+    dt*k3, t+    dt)                 */

int Rk4::Step(const double t, double *x){
/* no idea why the commented code has bug   */
/* 1.5 times of true value after 1 step     */
/*	double k1[dimen_];
	model_.rhs(t,        x, k1);

	double x2[dimen_];
	double k2[dimen_];
	for(int i=0;i<dimen_;i++) x2[i] = x[i] + 0.5*dt_*k1[i];
	model_.rhs(t+0.5*dt_,x2,k2);

	double x3[dimen_];
	double k3[dimen_];
	for(int i=0;i<dimen_;i++) x3[i] = x[i] + 0.5*dt_*k2[i];
	model_.rhs(t+0.5*dt_,x3,k3);

	double x4[dimen_];
	double k4[dimen_];
	for(int i=0;i<dimen_;i++) x4[i] = x[i] +     dt_*k3[i];
	model_.rhs(t+    dt_,x4,k4);

    for(int i=0;i<dimen_;i++)
    	x[i] += dt_/6.0*(k1[i]+k2[i]+k3[i]+k4[i]);

	return 0;*/
        double  k[4][dimen_]; // k[0][] points to k_1, k[1][] points to k_2 and so on
        double xk[3][dimen_]; // corresponds to x2, x3, x4

        int status = 0;
        status = model_.rhs(t,x,k[0]);

        if(status != 0) return 1;
        for(int i=0; i<dimen_; i++){
            xk[0][i] = x[i] + 0.5*dt_*k[0][i];
        }
        status = model_.rhs(t+0.5*dt_,xk[0],k[1]);
        if(status != 0) return 1;

        for(int i=0; i<dimen_; i++){
            xk[1][i] = x[i] + 0.5*dt_*k[1][i];
        }
        status = model_.rhs(t+0.5*dt_,xk[1],k[2]);
        if(status != 0) return 1;

        for(int i=0; i<dimen_; i++){
            xk[2][i] = x[i] +     dt_*k[2][i];
        }
        status = model_.rhs(t+    dt_,xk[2],k[3]);
        if(status != 0) return 1;

        for(int i = 0; i < dimen_; i++){
                x[i] += dt_/6*(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i]);
        }

        return 0;
}
