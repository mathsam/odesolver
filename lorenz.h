#ifndef LORENZ_H
#define LORENZ_H
#include "model.h"
/* Classic Lorenz system       */
/* \dot x = sigma_ * (y - x)   */
/* \dot y = rho_ * x - y - x*z */
/* \dot z = -beta_ * z + x*y   */

class Lorenz: public Model{
public:
	Lorenz(const double sigma,
			const double rho,const double beta);
	~Lorenz();
	int rhs(const double t,const double *x, double *fx) const;
	int dimen() const {return kDimen;}
private:
	const double sigma_;
	const double rho_;
	const double beta_;
	static const int kDimen = 3;
};


#endif  // LORENZ_H
