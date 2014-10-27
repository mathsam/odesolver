#!/usr/bin/env python
import os

integrators = ['euler','rk4','ab2'];
time_steps  = [0.1, 0.01, 0.001, 0.0001, 0.00001];
total_time  = 25;

for whic_integrator in integrators:
    for whic_time_step in time_steps:
        num_time = total_time/whic_time_step
        os.system('printf ' + str(whic_time_step)+'>>' + whic_integrator + '_conv.out')
        os.system('printf "  " >> ' +  whic_integrator + '_conv.out')
        os.system('./test_convergence ' + whic_integrator + ' ' \
                   + str(whic_time_step) + ' ' + str(num_time)  \
                   + '>>' + whic_integrator + '_conv.out')
#        os.system('printf "\n"' + '>>' + whic_integrator + '_conv.out')
