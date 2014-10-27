#include "euler.h"
#include "model.h"

Euler::Euler(double dt, const Model &model)
    : dimen_(model.dimen()),
      dt_(dt),
      model_(model) {
  fx_ = new double[dimen_];
}

Euler::~Euler() {
  delete [] fx_;
}

int Euler::Step(double t, double *x) {
  model_.rhs(t, x, fx_);
  for (int i = 0; i < dimen_; ++i) {
    x[i] += dt_ * fx_[i];
  }
  return 0;
}
