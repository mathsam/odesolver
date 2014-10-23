#include "linear-oscillator.h"
#include<math.h>

LinearOscillator::LinearOscillator(double beta,double gamma,double omega)
    :beta_(beta),
     gamma_(gamma),
     omega_(omega)
{}

LinearOscillator::~LinearOscillator()
{}

/* ode for linear oscillator:                           */
/* \dot \dot x + 2*beta*\dot x + x = gamma*cos(omega*t) */
/* x[0] -> x; x[1] -> \dot x                            */
int LinearOscillator::rhs(const double t,const double *x,
		double *fx) const{
    fx[0] = x[1];
    fx[1] = -2*beta_*x[1] - x[0] + gamma_*cos(omega_*t);
    return 0;
}
