#coding: utf-8
from utils import read_params, gen_str_params
import click
import numpy as np
import pylab_config
from itertools import product
from general_plot import number_equilibria


plt = pylab_config.setup()

params2latex = {'pi_A': r'$\pi_A$', 'pi_B': r'$\pi_B$',
                'xi': r'$\xi$', 'lpub': r'$\lambda_{pub}$',
                'lpriv': r'$\lambda_{priv}$', 'pa': r'$\pi_a$',
                'pb': r'$\pi_b$'}


def inside_interval(x, interval, tol=0.01):
    """
        >>> inside_interval(1.5, [4, 8])
        False
        >>> inside_interval(3.2, [1.8, 5.3])
        True
    """
    start, end = sorted(interval)
    return (x >= start-tol) and (x <= end+tol)


def get_unstable_point(list_eq, init_cond):
    """
        >>> list_eq = [1, 1, 2, 2]
        >>> init_cond = [0.2, 0.4, 0.6, 0.8]
        >>> get_unstable_point(list_eq, init_cond)
        0.5
    """
    previous_eq = list_eq[0]
    for i, eq in enumerate(list_eq[1:]):
        if round(eq, 3) != round(previous_eq, 3):
            break
        previous_eq = eq
    return (init_cond[i] + init_cond[i+1])/2


@click.command()
@click.option('--peer/--unanimity', default=False)
@click.option('--lpub', default=None)
@click.option('--debug/--normal', default=False)
def run(peer, lpub, debug):
    params_file = 'params_txtfiles/params_arrow_equilibria'
    if debug:
        params_file += '_debug'
    elif lpub:
        params_file += '_lpub_{}'.format(lpub)
    params_file += '.txt'
    all_params = read_params(params_file, int_params=['n'])
    n = all_params['n'][0]
    params_lab1 = ['pi_B', 'xi', 'pa', 'pb', 'lpub', 'n']
    params_lab2 = ['pi_A']
    list_set1 = [all_params[lab] for lab in params_lab1]
    list_set2 = [all_params[lab] for lab in params_lab2]
    for params_set1 in product(*list_set1):
        params = {lab: params_set1[i] for i, lab in enumerate(params_lab1)}
        fig = plt.figure(figsize=(7, 4.9), dpi=600)
        ax = fig.add_subplot(111)
        for params_set2 in product(*list_set2):
            for i, lab in enumerate(params_lab2):
                plotted_points = []
                params[lab] = params_set2[i]
                print(lab, params_set2[i])
                #print params
                neq, list_eq, init_cond = number_equilibria(params, peer=peer, n=n)
                list_eq = np.array(list_eq, dtype=float)

                # stable
                stable_points = list(set(list_eq[1:])) if neq == 1 else [min(list_eq), max(list_eq)]
                for eq_point in stable_points:
                    ax.scatter(params[lab], eq_point, s=50, c='#1C86EE', edgecolor='',
                            marker='o')

                plotted_points += list(stable_points)

                if neq > 1:
                    if neq == 3:
                        unstable_point = [eq_pt for eq_pt in list_eq if eq_pt not in stable_points][0]
                    else:
                        unstable_point = get_unstable_point(list_eq, init_cond)

                    ax.scatter(params[lab], unstable_point, s=50, c='#ff4040', edgecolor='',
                            marker='x', linewidth=3)

                    plotted_points.append(unstable_point)
                all_points = sorted([0.] + plotted_points + [1.])
                for i, plotted_pt in enumerate(all_points[:-1]):
                    curr_point = plotted_pt
                    next_point = all_points[i+1]
                    upwards = (i % 2 == 0)
                    diff = 0.2 if (curr_point == 0 or next_point == 1.) else next_point - curr_point
                    if diff < 0.15:
                        if diff < 0.06:
                            little_diff = 0.02 if diff > 0.035 else diff/3
                            head_width = 0.0025 if diff > 0.035 else 0.0015
                            head_length = 0.003 if diff > 0.035 else 0.0015
                        else:
                            little_diff = 0.03
                            head_width = 0.005
                            head_length = 0.01
                        if upwards:
                            start = curr_point + little_diff
                            end = next_point - little_diff
                        else:
                            start = next_point - little_diff
                            end = curr_point + little_diff
                        ax.arrow(params[lab], start, 0, end-start,
                            head_width=head_width, head_length=head_length,
                            color='black')
                    else:
                        if upwards:
                            get_start = lambda rank, pt: pt + 0.03 if rank == 0 else pt - 0.07
                            get_end = lambda rank, pt: pt + 0.07 if rank == 0 else pt - 0.03
                        else:
                            get_start = lambda rank, pt: pt + 0.07 if rank == 0 else pt - 0.03
                            get_end = lambda rank, pt: pt + 0.03 if rank == 0 else pt - 0.07
                        for j, pt in enumerate([curr_point, next_point]):
                            if pt in [0., 1]:
                                continue
                            start = get_start(j, pt)
                            start = max(start, 0) if upwards else min(start, 1)
                            end = get_end(j, pt)
                            end = min(end, 1) if upwards else max(end, 0)
                            ax.arrow(params[lab], start, 0, end-start,
                                head_width=0.005, head_length=0.01,
                                color='black')

        excluded_keys = ['lpriv', 'n'] + params_lab2
        str_params = gen_str_params(params, params_order=params_lab1, exclude_keys=excluded_keys)

        uniq_param = '_'.join(str_params.replace('\n', '')
                            .replace('$', '').split('; '))

        ax.grid(True)
        ax.set_xlabel(r'$\pi_A$')
        ax.set_ylabel(r'$z_A$')
        ax.set_xlim([-0.01, 1.01])
        ax.set_ylim([-0.01, 1.01])

        plt.tight_layout()

        path_to_fig_dir = 'figs/peer/number_equilibria' if peer else 'figs'
        plt.savefig('%s/equilibria_%s.pdf' % (path_to_fig_dir, uniq_param), bbox_inches='tight')

        titulo = gen_str_params(params, params_order=params_lab1, symbols=params2latex,
            exclude_keys=excluded_keys).replace('\n', '')
        ax.set_title(titulo, fontsize=30)
        plt.savefig('%s/equilibria_%s.png' % (path_to_fig_dir, uniq_param), bbox_inches='tight')


if __name__ == '__main__':
    run()
