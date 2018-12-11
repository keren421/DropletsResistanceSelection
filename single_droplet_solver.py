from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

#dn0/dt = -sigma(di*ni)
#dni/dt = lambda(n0)/c*ni*(1-sigma(nj))-k(no)*ni


def bacteria_growth(t, n):
    d_rate = [0, 0, 10, 100]
    mics = [0, 1, 10, 100]
    nominal_growth_rate = 1
    nominal_killing_rate = 1

    growth_rate = [nominal_growth_rate*np.heaviside(mic-n[0],0) for mic in mics]
    killing_rate = [nominal_killing_rate*np.heaviside(n[0] - mic,0) for mic in mics]
    C = 50
    total_bacteria = n[1] + n[2] + n[3]
    return [-n[1]*d_rate[1] - n[2]*d_rate[2] - n[3]*d_rate[3],
            growth_rate[1]/C*n[1]*(1-total_bacteria) - killing_rate[1]*n(1),
            growth_rate[2]/C*n[2]*(1-total_bacteria) - killing_rate[2]*n(2),
            growth_rate[3]/C*n[3]*(1-total_bacteria) - killing_rate[3]*n(3)]


sensitive_bacteria = {'amount': 5,
                      'mic': 1,
                      'd_rate': 0}
resistant_bacteria = {'amount': 1,
                      'mic': 10,
                      'd_rate': 10}
super_resistant_bacteria = {'amount': 1,
                            'mic': 100,
                            'd_rate': 100}

t_end = 20
initial_pop = [1, 5,0,1]
sol = solve_ivp(bacteria_growth, (0, t_end), initial_pop)
plt.plot(sol.t, sol.y)
plt.ylabel('Population Size')
plt.xlabel('Time')
plt.show()