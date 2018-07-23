import numpy as np
import pandas as pd

from utils import read_params, gen_str_params, dynamics_to_dataframe
from general_model import gen_init_conditions, integrate_unanimity_ode, params2latex
import pylab_config

plt = pylab_config.setup()


def main():
    integrate_ode = integrate_unanimity_ode
    params_file = 'params_txtfiles/params_unanimity_zA_vs_xi.txt'
    all_params = read_params(params_file, int_params=['n'])
    n = all_params['n'][0]
    #print all_params
    params_order = ['pi_A', 'pi_B', 'xi', 'lpub', 'lpriv', 'pa', 'pb']
    types_names = [''.join(['z']+['a']*(n-i)+['b']*i) for i in range(n+1)]
    tipo2latex = {elem: r'$z_{%s}$' % (elem[1:].upper()) for elem in types_names}
    my_colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
    more_than_one = ['xi']
    init_conditions = gen_init_conditions(n)[::2]
    str_vec_in = [str(vec_in) for vec_in in init_conditions]
    df_zA = pd.DataFrame()
    linestyles = ['-', '--', '-.', ':']
    handles = []
    labels = []
    fig = plt.figure(figsize=(16, 14))
    ax = fig.add_subplot(111)
    out_count = 0
    for xi in all_params['xi']:
        params = {lab: all_params[lab][0] for i, lab in enumerate(params_order) if lab != 'xi'}
        params['xi'] = xi
        print(params)
        lpriv = 1-params['lpub']
        params['lpriv'] = lpriv
        params_nums = tuple([params[lab] for lab in params_order]+[n])
        for i, vec_in in enumerate(init_conditions):
            df_temp = pd.DataFrame()
            curve_points, time_points = integrate_ode(vec_in, params_nums)
            df_temp = dynamics_to_dataframe(n, curve_points, time_points, initial_vec_str=str_vec_in[i])
            df_temp['xi'] = xi
            df_zA = pd.concat([df_zA, df_temp])

            df_aux = df_temp.sort_values(by='t')
            tipo = 'zA'
            freq_tipo = df_aux[tipo].tolist()
            time_points = df_aux['t'].tolist()
            linha, = ax.plot(time_points, freq_tipo, color=my_colors[i], lw=2, linestyle=linestyles[out_count])
            out_count += 1
            handles.append(linha)

            condicao = []
            for i, elem in enumerate(vec_in[:1]):
                pedaco_str = tipo2latex[types_names[i]] + \
                    (r'$^{init}=%.0f' % round(elem, 0))+'$'
                condicao.append(pedaco_str)

            pedaco_str = r'$z_{BBB}$' + (r'$^{init}=%.0f' % round(1-sum(vec_in), 0))+'$'
            condicao.append(pedaco_str)
            pedaco_str = (r'$ \xi=%.7f' % xi).rstrip('0').rstrip('.')+'$'
            condicao.append(pedaco_str)
            labels.append(';\t'.join(condicao))

    fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.17), loc='center')
    default_order = ['lpub', 'lpriv', 'pi_A', 'pi_B', 'pa', 'pb', 'xi']
    default_order = [col for col in default_order if col not in more_than_one]
    titulo = gen_str_params(params, params_order=default_order, symbols=params2latex,
        exclude_keys=more_than_one).replace('\n', '')

    max_za = np.max(df_zA['zA'].tolist())
    min_za = np.min(df_zA['zA'].tolist())
    dist = 0.05*(max_za-min_za)

    ax.set_xlim([0., 160.1])
    #ax.set_ylim([0., 0.81])
    ax.set_ylim([min_za-dist, max_za+dist])
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$z_A$')
    ax.grid(True)
    plt.subplots_adjust(left=0., right=0.95, top=0.75, bottom=0.4)
    str_params = gen_str_params(params, params_order=params_order, exclude_keys=more_than_one)

    uniq_param = '_'.join(str_params.replace('\n', '').replace('$', '').split('; '))
    uniq_param += '_n=%d' % n

    plt.savefig('figs/same_plot_zA_dynamics_%s.pdf' % (uniq_param), bbox_inches='tight')
    ax.set_title(titulo)

    #plt.tight_layout()
    plt.savefig('figs/same_plot_zA_dynamics_%s.png' % (uniq_param), bbox_inches='tight')

    plt.close()


if __name__ == '__main__':
    main()
