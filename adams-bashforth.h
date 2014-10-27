#ifndef ADAMS_BASHFORTH_H_
#define ADAMS_BASHFORTH_H_
#include "integrator.h"
class Model;

class AdamsBashforth: public Integrator{
public:
	AdamsBashforth(double dt, const Model &model);
	~AdamsBashforth();
	int Step(const double t, double *x);
private:
	const int dimen_;
	const double dt_;
	const Model &model_;
	double *fx_past_;     // f(x_j-1,t_j-1)
	bool  is_first_step_; // whether this is first time step
};


# endif //ADAMS_BASHFORTH_H
