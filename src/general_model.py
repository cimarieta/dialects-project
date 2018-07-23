#coding: utf-8
"""
NOT STOCHASTIC!
Code to the private dimension of the model
- integrates ode
vec[0] = allB
vec[n] = allA
"""
import numpy as np
from scipy.integrate import odeint


class DialectDynamicsSimulator(object):
    def __init__(self, params, conformity_mode='peer'):
        self.params = params
        self.conformity_mode = conformity_mode
        self.private_dz = self.peer_private_dz if conformity_mode == 'peer' else self.unanimity_private_dz

    def calc_zA(self, vec, n):
        return (np.sum([i*vec[i] for i in range(n)]) + n*(1.-np.sum(vec)))/n

    def common_dz(self, vec, term_A, term_B):
        n = self.params['n']
        if len(vec) == n:
            curr_vec = np.empty(len(vec)+1)
            curr_vec[:-1] = np.copy(vec)
            curr_vec[-1] = 1-np.sum(vec)
        else:
            curr_vec = np.copy(vec)
        diff_vec = np.empty(n)  # + 1)
        diff_vec[0] = -term_A(0)*curr_vec[0] + term_B(1)*curr_vec[1]
        for i in range(1, n):
            diff_vec[i] = term_A(i-1)*curr_vec[i-1] - (term_A(i) + term_B(i))*curr_vec[i] + \
                term_B(i+1)*curr_vec[i+1]
        #diff_vec[n] = -term_B(n)*curr_vec[n] + term_A(n-1)*curr_vec[n-1]
        return diff_vec

    def peer_private_dz(self, vec):
        n = self.params['n']
        xi = self.params['xi']
        pa = self.params['pa']
        pb = self.params['pb']
        term_pa = lambda i: (n-i)*(pa+i*xi/(n-1))
        term_pb = lambda i: i*(pb+(n-i)*xi/(n-1))
        diff_vec = self.common_dz(vec, term_pa, term_pb)
        return diff_vec

    def unanimity_private_dz(self, vec):
        n = self.params['n']
        xi = self.params['xi']
        pa = self.params['pa']
        pb = self.params['pb']
        term_pa = lambda i: (n-i)*pa if i < n - 1 else (n-i)*pa + xi
        term_pb = lambda i: i*pb if i > 1 else i*pb + xi
        diff_vec = self.common_dz(vec, term_pa, term_pb)
        return diff_vec

    def public_dz(self, vec):
        n = self.params['n']
        xi = self.params['xi']
        pi_A = self.params['pi_A']
        pi_B = self.params['pi_B']
        zA = self.calc_zA(vec, n)
        term_pi_A = lambda i: (n-i)*(pi_A+xi*zA)
        term_pi_B = lambda i: i*(pi_B+xi*(1-zA))
        diff_vec = self.common_dz(vec, term_pi_A, term_pi_B)
        return diff_vec

    def public_private_dz(self, vec):
        lpub = self.params['lpub']
        lpriv = self.params['lpriv']
        return lpub*self.public_dz(vec) + lpriv*self.private_dz(vec)


def _calc_peer_freqA(vec, pi_A, pi_B, xi, lpub, lpriv, pa, pb, n):
    params = {'pi_A': pi_A, 'pi_B': pi_B, 'pa': pa, 'pb': pb,
              'xi': xi, 'lpub': lpub, 'lpriv': lpriv, 'n': n}
    dialect_simulator = DialectDynamicsSimulator(params, conformity_mode='peer')
    return dialect_simulator.public_private_dz(vec)


def calc_peer_freqA(vec, t, pi_A, pi_B, xi, lpub, lpriv, pa, pb, n):
    return _calc_peer_freqA(vec, pi_A, pi_B, xi, lpub, lpriv, pa, pb, n)


def _calc_unanimity_freqA(vec, pi_A, pi_B, xi, lpub, lpriv, pa, pb, n):
    params = {'pi_A': pi_A, 'pi_B': pi_B, 'pa': pa, 'pb': pb,
              'xi': xi, 'lpub': lpub, 'lpriv': lpriv, 'n': n}
    dialect_simulator = DialectDynamicsSimulator(params, conformity_mode='unanimity')
    return dialect_simulator.public_private_dz(vec)


def calc_unanimity_freqA(vec, t, pi_A, pi_B, xi, lpub, lpriv, pa, pb, n):
    return _calc_unanimity_freqA(vec, pi_A, pi_B, xi, lpub, lpriv, pa, pb, n)


def integrate_peer_ode(initial_freqs, params_nums):
    return integrate_ode_select_fun(initial_freqs, params_nums, calc_peer_freqA)


def integrate_unanimity_ode(initial_freqs, params_nums):
    return integrate_ode_select_fun(initial_freqs, params_nums, calc_unanimity_freqA)


def integrate_ode_select_fun(initial_freqs, params_nums, function):
    t0=0
    tEnd=200
    dt=0.01
    t_values = np.arange(t0, tEnd, dt)
    freq_values = np.array(odeint(function, initial_freqs, t_values, args=params_nums))
    return freq_values, t_values


def gen_init_conditions(n):
    initial_freqs1 = np.zeros(n)
    initial_freqs2 = np.zeros(n)
    #initial_freqs3 = np.array([1./(n+1) for i in range(n)])
    initial_freqs4 = np.zeros(n)

    initial_freqs1[0] = 0.15
    initial_freqs2[0] = 0.35
    initial_freqs4[0] = 0.5

    initial_freqs1[-1] = 0.35
    initial_freqs2[-1] = 0.15
    initial_freqs4[-1] = 0.5

    for i in range(1, n-1):
        initial_freqs1[i] = (1.-initial_freqs1[0]-initial_freqs1[-1])/(n-1)
        initial_freqs2[i] = (1.-initial_freqs2[0]-initial_freqs2[-1])/(n-1)

    #init_conditions = [initial_freqs1, initial_freqs2, initial_freqs3, initial_freqs4]
    #init_conditions = [initial_freqs3, initial_freqs2]
    init_conditions = [[0.99, 0., 0.], [0.5, 0., 0.], [0.01, 0., 0.]]
    return init_conditions
