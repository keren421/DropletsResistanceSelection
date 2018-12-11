from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
#plt.ion()

#dn0/dt = -n0 * sigma(di*ni)
#dni/dt = lambda(n0)*ni*(1-sigma(nj)/c)-k(no)*ni


def bacteria_growth(t, n, d_rate, mics):
    nominal_growth_rate = 1
    nominal_killing_rate = 10 #
    C = 50
    total_bacteria = sum(n[1:])

    derivatives = []
    for j in range(len(n)):
        if j == 0:
            derivatives.append(-n[0] * sum([n[bacteria_num]*d_rate[bacteria_num] for bacteria_num in range(1,len(n))]))
        else:
            growth_rate = nominal_growth_rate*np.heaviside(mics[j]-n[0], 1)
            killing_rate = nominal_killing_rate*np.heaviside(n[0] - mics[j], 1)
            derivatives.append(growth_rate*n[j]*(1-total_bacteria/C) - killing_rate*n[j])
    return derivatives


def single_droplet_solver(initial_pop, d_rate, mics):
    t_end = 20
    growth_in_droplet = lambda t,n: bacteria_growth(t, n, d_rate, mics)
    sol = solve_ivp(growth_in_droplet, (0, t_end), initial_pop)
    plot_solution(sol)
    return sol


def plot_solution(sol):
    f = plt.figure()
    plt.plot(sol.t, sol.y[0], label='Antibiotic')
    plt.plot(sol.t, sol.y[1], label='Sensitive')
    plt.plot(sol.t, sol.y[2], label='Resistant')
    plt.plot(sol.t, sol.y[3], label='Super Resistance')
    plt.legend()
    plt.ylabel('Population Size')
    plt.xlabel('Time')


initial_pop = [8, 5,  1,   0]
d_rate = [0, 0,  1,  10]
mics = [0, 1, 10, 100]
sol = single_droplet_solver([11, 5, 1, 0], d_rate, mics)
#sol = single_droplet_solver([8, 5, 1, 0], d_rate, mics)
print(sol.y[:,-1])
plt.show()