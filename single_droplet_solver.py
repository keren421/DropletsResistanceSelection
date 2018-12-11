from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

#dn0/dt = -sigma(di*ni)
#dni/dt = lambda(n0)*ni*(1-sigma(nj)/c)-k(no)*ni


def bacteria_growth(t, n,d_rate, mics):
    nominal_growth_rate = 1
    nominal_killing_rate = 10

    growth_rate = [nominal_growth_rate*np.heaviside(mic-n[0],1) for mic in mics]
    killing_rate = [nominal_killing_rate*np.heaviside(n[0] - mic,1) for mic in mics]
    C = 50
    total_bacteria = n[1] + n[2] + n[3]
    return [-n[0] * (n[1]*d_rate[1] + n[2]*d_rate[2] + n[3]*d_rate[3]),
            growth_rate[1]*n[1]*(1-total_bacteria/C) - killing_rate[1]*n[1],
            growth_rate[2]*n[2]*(1-total_bacteria/C) - killing_rate[2]*n[2],
            growth_rate[3]*n[3]*(1-total_bacteria/C) - killing_rate[3]*n[3]]


def single_droplet_solver(initial_pop, d_rate, mics):
    t_end = 20
    growth_in_droplet = lambda t,n: bacteria_growth(t, n, d_rate, mics)
    sol = solve_ivp(growth_in_droplet, (0, t_end), initial_pop)
    return sol


def plot_solution(sol):
    plt.figure(1)
    plt.plot(sol.t, sol.y[0], label='Antibiotic')
    plt.plot(sol.t, sol.y[1], label='Sensitive')
    plt.plot(sol.t, sol.y[2], label='Resistant')
    plt.plot(sol.t, sol.y[3], label='Super Resistance')
    plt.legend()
    plt.ylabel('Population Size')
    plt.xlabel('Time')
    plt.show()

sensitive_bacteria = {'amount': 5,
                      'mic': 1,
                      'd_rate': 0}
resistant_bacteria = {'amount': 1,
                      'mic': 10,
                      'd_rate': 10}
super_resistant_bacteria = {'amount': 1,
                            'mic': 100,
                            'd_rate': 100}

initial_pop = [8, 5, 1,0]
d_rate = [0, 0, 10, 100]
mics = [0, 1, 10, 100]
sol = single_droplet_solver(initial_pop, d_rate, mics)
plot_solution(sol)
