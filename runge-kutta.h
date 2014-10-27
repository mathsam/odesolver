#ifndef RK4_H_
#define RK4_H_
#include "integrator.h"
class Model;

class Rk4: public Integrator {
public:
	Rk4(double dt, const Model &model);
	~Rk4();
	int Step(const double t, double *x);
private:
	const int dimen_;
	const double dt_;
	const Model &model_;
};


#endif //RK4_H
