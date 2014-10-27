integrators = euler.o runge-kutta.o adams-bashforth.o
equations = duffing.o lorenz.o linear-oscillator.o
objects = ode_solve.cc $(integrators) $(equations)
test_obj= test_convergence.cc $(integrators) linear-oscillator.o

CXXFLAGS = -g -Wall -std=c++11
CXX      = g++

all: ode_solve test_convergence

ode_solve : $(objects)
	$(CXX) -o $@ $^

test_convergence: $(test_obj)
	$(CXX) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
