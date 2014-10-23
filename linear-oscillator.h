#ifndef LINEAR_OSCILLATOR_H
#define LINEAR_OSCILLATOR_H
#include "model.h"

/* Linear oscillator                                */
/* \dot x = y                                       */
/* \dot y = -2*beta_*y - x - gamma_*cos(omega_*t)   */

class LinearOscillator: public Model{
public:
	LinearOscillator(const double beta,const double gamma,const double omega);
	~LinearOscillator();
	int rhs(const double t,const double *x, double *fx) const;
	int dimen() const {return kDimen;}
private:
	const double beta_;
	const double gamma_;
	const double omega_;
	static const int kDimen = 2;
};

#endif //LINEAR_OSCILLATOR_H
