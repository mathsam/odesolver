#ifndef MODEL_H_
#define MODEL_H_

class Model {
 public:
  virtual ~Model() {}

  // fx = f(x,t)
  virtual int rhs(double t, const double *x, double *fx) const = 0;

  // number of states (size of x)
  virtual int dimen() const = 0;
};

#endif  // MODEL_H_
