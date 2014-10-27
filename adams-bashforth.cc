#include "adams-bashforth.h"
#include "model.h"
#include "runge-kutta.h"

AdamsBashforth::AdamsBashforth(double dt, const Model &model)
    : dimen_(model.dimen()),
      dt_(dt),
      model_(model),
      is_first_step_(true){
	fx_past_ = new double[dimen_];
}

AdamsBashforth::~AdamsBashforth(){
	delete [] fx_past_;
}

/* ODE solver using 2rd order Adams-Bashforth             */
/* for an ODE dx/dt = f(x,t), it is integrated as         */
/* x_j+1 = x_j + 1.5*dt*f(x_j,t_j)-0.5*dt*f(x_j-1,t_j-1)  */

int AdamsBashforth::Step(double t,double *x){
	if(is_first_step_){
        Rk4 integ_first_step = Rk4(dt_,model_);
        int status = integ_first_step.Step(t,x);
        if(status!=0) return status;

        model_.rhs(t,x,fx_past_);
        is_first_step_ = false;
        return status;
	}
	else{
        double fx[dimen_];
        model_.rhs(t,x,fx);
        for(int i=0;i<dimen_;i++){
        	x[i] += dt_*(1.5*fx[i] -0.5*fx_past_[i]);
        }
        for(int i=0;i<dimen_;i++){
        	fx_past_[i] = fx[i];
        }
        return 0;
	}
}
