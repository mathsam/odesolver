#ifndef DUFFING_H_
#define DUFFING_H_

#include "model.h"

// Forced Duffing oscillator:
//   \dot x = y
//   \dot y = x - x^3 - delta * y + gamma * cos(omega t)
class Duffing : public Model {
 public:
  Duffing(double delta, double gamma, double omega);
  ~Duffing();
  int rhs(double t, const double *x, double *fx) const;
  int dimen() const { return kDimen; }
 private:
  const double delta_;                        // damping constant
  const double gamma_;                        // amplitude of forcing
  const double omega_;                        // frequency of forcing
  static const int kDimen = 2;
};

#endif  // DUFFING_H_
