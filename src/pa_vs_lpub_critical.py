#coding: utf-8
from utils import read_params, gen_str_params
import numpy as np
from itertools import product
from scipy.optimize import brentq

import pylab_config

plt = pylab_config.setup()
params2latex = {'pi_A': r'$\pi_A$', 'pi_B': r'$\pi_B$',
                'xi': r'$\xi$', 'lpub': r'$\lambda_{pub}$',
                'lpriv': r'$\lambda_{priv}$', 'pa': r'$\pi_a$',
                'pb': r'$\pi_b$'}


def analytical_eq_point(lpub, pa, pb, pi_A, pi_B):
    lpriv = 1. - lpub
    pub_ctb = lpub*(pi_A + pi_B)
    priv_ctb = lpriv*(pa + pb)
    zA = (lpub*pi_A + lpriv*pa)/(pub_ctb + priv_ctb)
    return zA - 0.5


def find_critical(pa, pb, pi_A, pi_B):
    a = 0.
    b = 1.
    point_a = analytical_eq_point(a, pa, pb, pi_A, pi_B)
    point_b = analytical_eq_point(b, pa, pb, pi_A, pi_B)
    if point_a*point_b > 0:
        if point_a > 0:
            return 0.
        else:
            return 1.
    return brentq(analytical_eq_point, 0., 1., args=(pa, pb, pi_A, pi_B))


def plot():
    params_file = 'params_txtfiles/params_peer_critical_lpub.txt'
    all_params = read_params(params_file, int_params=['n'])
    pi_A, pi_B = all_params['pi_A'][0], all_params['pi_B'][0]
    params_order = ['pi_A', 'pi_B']
    params_to_use = [all_params[lab] for lab in params_order]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for params_set in product(*params_to_use):
        params = {lab: params_set[i] for i, lab in enumerate(params_order)}
        all_lpub_values = []
        colors = ['#ffb347', '#b347ff', '#47ffb3']
        linestyles = ['--', '-', '-.']
        for i, pb in enumerate([0.3, 0.5, 0.7]):

            pa_values = np.arange(0.01, 1.01, 0.01)
            # find critical lpub for each pa_value
            lpub_values = [find_critical(pa, pb, pi_A, pi_B) for pa in pa_values]

            ax.plot(pa_values, lpub_values, color=colors[i], label=r'$\pi_b=%.1f$' % pb,
                linestyle=linestyles[i])

            ax.grid(True, linestyle=':')
            all_lpub_values += lpub_values

        #ax.legend(bbox_to_anchor=(0.5,0.17), loc='center')
        ax.legend()
        title = gen_str_params(params, params_order=params_order, symbols=params2latex).replace('\n', '')

        ax.set_xlabel(r'$\pi_a$')
        ax.set_ylabel(r'critical $\lambda_{pub}$')

        ax.set_xlim([-0.05, 1.05])
        ax.set_ylim([-0.05, np.max(all_lpub_values)+0.05])
        ax.set_xticks([0., 0.2, 0.4, 0.6, 0.8, 1.])

        plt.tight_layout()

        str_params = gen_str_params(params, params_order=params_order)
        uniq_param = '_'.join(str_params.replace('\n', '').replace('$', '').split('; '))

        plt.savefig('figs/peer/critical_lpub_%s.pdf' % uniq_param)
        ax.set_title(title)
        plt.savefig('figs/peer/critical_lpub_%s.png' % uniq_param)


if __name__ == '__main__':
    plot()
