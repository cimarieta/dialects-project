#coding: utf-8
from utils import read_params, gen_str_params, dynamics_to_dataframe
import click
import numpy as np
import pandas as pd
from itertools import product

from general_model import gen_init_conditions, integrate_peer_ode, integrate_unanimity_ode
import pylab_config

plt = pylab_config.setup()

params2latex = {'pi_A': r'$\pi_A$', 'pi_B': r'$\pi_B$',
                'xi': r'$\xi$', 'lpub': r'$\lambda_{pub}$',
                'lpriv': r'$\lambda_{priv}$', 'pa': r'$\pi_a$',
                'pb': r'$\pi_b$'}

lpub_values = [0., 0.3, 0.7, 1.]


@click.group()
def cli():
    pass


@cli.command()
@click.option('--peer/--unanimity', default=False)
def plot_several_lpubs(peer=False):
    integrate_ode = integrate_peer_ode if peer else integrate_unanimity_ode
    if peer:
        params_file = 'params_txtfiles/params_peer_zA_vs_lpub.txt'
    else:
        params_file = 'params_txtfiles/params_unanimity_zA_vs_lpub.txt'
    all_params = read_params(params_file, int_params=['n'])
    n = all_params['n'][0]
    init_conditions = gen_init_conditions(n)
    #print all_params
    params_order = ['pi_A', 'pi_B', 'xi', 'pa', 'pb']
    l = [all_params[lab] for lab in params_order]
    types_names = [''.join(['z']+['a']*(n-i)+['b']*i) for i in range(n+1)]
    my_colors = ['#e41a1c', '#45cea2', '#377eb8', '#984ea3']
    linestyles = ['--', '-', ':', '-.']
    for params_set in product(*l):
        params = {lab: params_set[i] for i, lab in \
                    enumerate(params_order)}
        full_ordered_params = ['pi_A', 'pi_B', 'xi', 'lpub', 'lpriv', 'pa', 'pb']
        for vec_in in init_conditions:
            fig = plt.figure(figsize=(10,6), dpi=600)
            ax = fig.add_subplot(111)
            maximos = []
            for i, lpub in enumerate(lpub_values):
                params['lpub'] = lpub
                lpriv = 1-params['lpub']
                params['lpriv'] = lpriv
                params_nums = tuple([params[lab] for lab in full_ordered_params]+[n])

                curve_points, time_points = integrate_ode(vec_in, params_nums)
               
                df_temp = dynamics_to_dataframe(n, curve_points, time_points, initial_vec_str=None)
                lpub_label = (r'$\lambda_{pub}=%.7f'%lpub).rstrip('0').rstrip('.')+'$'
                ax.plot(time_points, df_temp['zA'].tolist(), label=lpub_label, color=my_colors[i%4],
                        linestyle=linestyles[i%4], lw=1.8)
                max_y = np.max([np.max(df_temp[group_type].tolist()) for group_type in types_names])
                maximos.append(max_y)

            max_y = np.max(maximos)
            if params['xi']==10 or peer:
                ax.legend(ncol=2, bbox_to_anchor=(0.5,-0.57), loc='center')
            default_order = ['pi_A', 'pi_B', 'pa', 'pb', 'xi']
            excluded_keys = ['lpub', 'lpriv', 'n']
            titulo = gen_str_params(params, \
                        params_order=default_order, symbols=params2latex, exclude_keys=excluded_keys).replace('\n', '')
            #print titulo

            #ax.set_xlim([0.,max(time_points)])
            ax.set_xlim([0.,50])
            ax.set_ylim([-0.01, 1.01])
            ax.set_xlabel(r'$t$')
            ax.set_ylabel(r'$z_A$')
            ax.grid(True)

            #plt.tight_layout()
            plt.subplots_adjust(left=0., right=0.95, top=0.95, bottom=0.55)
            str_params = gen_str_params(params, params_order=params_order, exclude_keys=excluded_keys)

            init_za = (3*vec_in[0]+2*vec_in[1]+vec_in[2])/3.

            uniq_param = '_'.join(str_params.replace('\n', '') \
                                .replace('$','').split('; '))
            uniq_param += '_n=%d_za=%.2f'%(n,init_za)

            path_to_fig_dir = 'figs/peer' if peer else 'figs'
            plt.savefig('%s/different_lpubs_zA_dynamics_%s.pdf'%(path_to_fig_dir, uniq_param),\
                 bbox_inches='tight')
            ax.set_title(titulo, fontsize=30)
            plt.savefig('%s/different_lpubs_zA_dynamics_%s.png'%(path_to_fig_dir, uniq_param),\
                 bbox_inches='tight')

            plt.close()


@cli.command()
@click.option('--peer/--unanimity', default=False)
def plot_several_lpubs_zis(peer=False):
    integrate_ode = integrate_peer_ode if peer else integrate_unanimity_ode   
    if peer:
        params_file = 'params_txtfiles/params_peer_zA_vs_lpub.txt'
    else:
        params_file = 'params_txtfiles/params_unanimity_zA_vs_lpub.txt'
    all_params = read_params(params_file, int_params=['n'])
    n = all_params['n'][0]
    init_conditions = gen_init_conditions(n)
    #print all_params
    params_order = ['pi_A', 'pi_B', 'xi', 'pa', 'pb']
    l = [all_params[lab] for lab in params_order]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    types_names = [''.join(['z']+['a']*(n-i)+['b']*i) for i in range(n+1)]
    tipo2latex = {elem: r'$z_%d$'%(elem.count('a')) for elem in types_names}
    my_colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
    linestyles = ['--', '-', ':', '-.']
    for params_set in product(*l):
        params = {lab: params_set[i] for i, lab in \
                    enumerate(params_order)}

        params_order_fun = ['pi_A', 'pi_B', 'xi', 'lpub', 'lpriv', 'pa', 'pb']
        for vec_in in init_conditions:
            fig = plt.figure(figsize=(10,12), dpi=600)
            df = pd.DataFrame()
            for i, lpub in enumerate(lpub_values):
                params['lpub'] = lpub
                lpriv = 1-params['lpub']
                params['lpriv'] = lpriv
                params_nums = tuple([params[lab] for lab in params_order_fun]+[n])

                curve_points, time_points = integrate_ode(vec_in, params_nums)
               
                df_temp = dynamics_to_dataframe(n, curve_points, time_points, initial_vec_str=None)
                df_temp['lpub'] = lpub
                df = pd.concat([df, df_temp])

            for i, col in enumerate(types_names[::-1]):
                ax = fig.add_subplot(2,2,i+1)
                lines = []
                labels = []
                for cont, lpub in enumerate(lpub_values):
                    df_aux = df[df['lpub']==lpub]
                    df_aux = df_aux.sort_values(by='t')
                    y_values = df_aux[col].tolist()
                    lpub_label = (r'$\lambda_{pub}=%.7f'%lpub).rstrip('0').rstrip('.')+'$'
                    line, = ax.plot(time_points, y_values, color=my_colors[cont%4], linestyle=linestyles[cont%4])
                    lines.append(line)
                    labels.append(lpub_label)

                #ax.set_xlim([0.,max(time_points)])
                ax.set_xlim([0.,50.1])
                #ax.set_ylim([-0.01, 1.01])
                ax.set_xlabel(r'$t$')
                ax.set_ylabel(r'%s'%tipo2latex[col])
                ax.grid(True)

            default_order = ['pi_A', 'pi_B', 'pa', 'pb', 'xi']
            excluded_keys = ['lpub', 'lpriv', 'n']
            titulo = gen_str_params(params, \
                        params_order=default_order, symbols=params2latex, exclude_keys=excluded_keys).replace('\n', '')
            #print titulo

            plt.tight_layout()

            if params['xi'] > 1:
                plt.legend(lines, labels, bbox_to_anchor=(0.,-0.7), loc='center', ncol=2)
            plt.subplots_adjust(left=0., right=0.95, top=0.92, bottom=0.55, wspace=0.3)
            str_params = gen_str_params(params, params_order=params_order, exclude_keys=excluded_keys)

            init_za = (vec_in[1]+2*vec_in[2]+3*(1-sum(vec_in)))/3.

            uniq_param = '_'.join(str_params.replace('\n', '') \
                                .replace('$','').split('; '))
            uniq_param += '_n=%d_za=%.2f'%(n,init_za)

            path_to_fig_dir = 'figs/peer' if peer else 'figs'
            plt.savefig('%s/zis_lpubs_zA_dynamics_%s.pdf'%(path_to_fig_dir, uniq_param),\
                 bbox_inches='tight')
            plt.suptitle(titulo, fontsize=30)
            plt.savefig('%s/zis_lpubs_zA_dynamics_%s.png'%(path_to_fig_dir, uniq_param),\
                 bbox_inches='tight')

            plt.close()


def find_za(init, params, peer=False, n=3):
    integrate_ode = integrate_peer_ode if peer else integrate_unanimity_ode   
    params['lpriv'] = 1 - params['lpub']

    params_order_fun = ['pi_A', 'pi_B', 'xi', 'lpub', 'lpriv', 'pa', 'pb']
    params_nums = tuple([params[lab] for lab in params_order_fun]+[n])

    vec_in = [1-init, 0., 0.]
    curve_points, time_points = integrate_ode(vec_in, params_nums)
    last_point = curve_points[-1]
    zA = (np.sum([i*last_point[i] for i in range(n)]) + n*(1-np.sum(last_point)))/n
    return zA


def number_equilibria(params, peer=False, n=3):
    # zA level
    init_conditions = np.linspace(0.01, 0.99, 100)

    list_zA = []
    for init_zA in init_conditions:
        zA = find_za(init_zA, params, peer=peer, n=n)
        list_zA.append(zA)

    #print list_zA
    unique_zA = set([round(zA, 3) for zA in list_zA])
    list_zA = [str(round(zA, 5)) for zA in list_zA]
    return len(unique_zA), list_zA, init_conditions


@cli.command()
@click.option('--peer/--unanimity', default=False)
def bar_plot(peer=False, single=True):
    integrate_ode = integrate_peer_ode if peer else integrate_unanimity_ode
    if peer:
        params_file = 'params_txtfiles/params_peer_barplot.txt'
    else:
        params_file = 'params_txtfiles/params_unanimity_barplot.txt'
    all_params = read_params(params_file, int_params=['n'])
    n = all_params['n'][0]
    init_conditions = gen_init_conditions(n)
    #print all_params
    params_order = ['lpub', 'pi_A', 'pi_B', 'xi', 'pa', 'pb']
   
    l = [all_params[lab] for lab in params_order]
    types_names = [''.join(['z']+['a']*i+['b']*(n-i)) for i in range(n+1)]
    #my_colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
    params_order_fun = ['pi_A', 'pi_B', 'xi', 'lpub', 'lpriv', 'pa', 'pb']
    df_zA = pd.DataFrame()
    for params_set in product(*l):
        params = {lab: params_set[i] for i, lab in \
                    enumerate(params_order)}

        lpriv = 1-params['lpub']
        params['lpriv'] = lpriv
        params_nums = tuple([params[lab] for lab in params_order_fun]+[n])
        df_temp_zA = pd.DataFrame()
        df = pd.DataFrame(columns=['init', 't', 'zaaa', 'zaab', 'zabb', 'zbbb'])
        str_vec_in = [str(vec_in) for vec_in in init_conditions]
        for i, vec_in in enumerate(init_conditions[:1]):
            df_temp = pd.DataFrame()
            curve_points, time_points = integrate_ode(vec_in, params_nums)
            df_temp = dynamics_to_dataframe(n, curve_points, time_points, initial_vec_str=str_vec_in[i])
            df = pd.concat([df, df_temp])

        df_temp_zA['zA'] = df.apply(lambda row: row['zaaa']+2*row['zaab']/3 +\
                                    row['zabb']/3, axis=1)
        df_temp_zA['t'] = df['t'].tolist()
        for lab in params_order:
            df_temp_zA[lab] = [params[lab]]*len(df_temp_zA)
        df_temp_zA['init'] = df['init']
        df_zA = pd.concat([df_zA, df_temp_zA])

        df_aux = df[df['init']==str_vec_in[0]]

        fig_bar = plt.figure(figsize=(6,8), dpi=600)
        ax_bar = fig_bar.add_subplot(111)

        freq_types = [df_aux[group_type].tolist()[-1] for group_type in types_names]
        #print freq_types
       
        width = 1.
        ax_bar.bar(np.arange(n+1), freq_types, width, color='#499DF5', edgecolor='black')

        zA_value = np.dot(np.arange(n+1), freq_types)/3.
        ax_bar.axhline(y=zA_value, \
                        c='#EE30A7', lw=3., ls='dashed')
        ax_bar.text(3.8, zA_value-0.007, '%.2f'%zA_value, fontsize=24,
                        color='#EE30A7')
        ax_bar.set_ylim([0.,1.])
        ax_bar.set_ylabel('freq.')
        ax_bar.set_xticks(np.arange(n+1)+width/2)
       
        xtickNames = ax_bar.set_xticklabels([r'$z_%d$'%num for num \
                                        in range(n+1)])

        default_order = ['lpub', 'lpriv', 'pi_A', 'pi_B', 'pa', 'pb', 'xi']

        excluded_keys = ['lpub', 'lpriv', 'pi_A', 'pi_B', 'n'] \
                        if params['lpub']<0.01 else ['n']

        titulo = gen_str_params(params, \
                    params_order=default_order, symbols=params2latex,
                    exclude_keys=excluded_keys).replace('\n', '')
        #print titulo
        str_params = gen_str_params(params, params_order=params_order,
                        exclude_keys=excluded_keys)

        uniq_param = '_'.join(str_params.replace('\n', '') \
                            .replace('$','').split('; '))

        plt.setp(xtickNames, fontsize=24)
        plt.subplots_adjust(left=0.16, right=0.88)
        path_to_fig_dir = 'figs/peer' if peer else 'figs'
        plt.savefig('%s/bar_plots_%s.pdf'%(path_to_fig_dir, uniq_param))
        ax_bar.set_title(titulo)
        plt.savefig('%s/bar_plots_%s.png'%(path_to_fig_dir, uniq_param))
        fig_bar.clear()

@cli.command()
@click.option('--peer/--unanimity', default=False)
@click.option('--force_private/--use_lpub_from_file', default=False)
@click.option('--var', default='pa')
@click.option('--other_var', default='lpub')
def plot_za_vs_var(peer, force_private, var, other_var):
    integrate_ode = integrate_peer_ode if peer else integrate_unanimity_ode   
    if peer:
        params_file = 'params_txtfiles/params_peer_barplot.txt'
    else:
        params_file = 'params_txtfiles/params_unanimity_barplot.txt'
    all_params = read_params(params_file, int_params=['n'])
    if force_private:
        all_params['lpub'] = [0.]
        all_params['lpriv'] = [1.]

    n = all_params['n'][0]
    init_conditions = gen_init_conditions(n)
    params_order = ['pi_A', 'pi_B', 'xi', 'pa', 'pb', 'lpub', 'lpriv']
    params_lab1 = [lab for lab in params_order if lab.lower()\
                     not in [other_var, var]]
    params_lab2 = [var]
    list_set1 = [all_params[lab] for lab in params_lab1]
    list_set2 = [np.linspace(0.01,1.,20)]
    other_var_values = all_params[other_var]
    my_colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
    markers_lst = ['o', 'D', 'v', '*']
    if force_private:
        fig = plt.figure(figsize=(8,7), dpi=600)
        ax = fig.add_subplot(111)
    for params_set1 in product(*list_set1):
        if not force_private:
            fig = plt.figure(figsize=(8,6), dpi=600)
            ax = fig.add_subplot(111)
        lines = []
        labels=[]

        params = {lab: params_set1[i] for i, lab in \
                    enumerate(params_lab1)}
        for vec_in in init_conditions[:1]:
            for cnt, other_var_val in enumerate(other_var_values):
                values_y = []
                if other_var == 'lpub':
                    params['lpub'] = other_var_val
                    params['lpriv'] = 1 - params['lpub']
                else:
                    params[other_var] = other_var_val
                for params_set2 in product(*list_set2):
                    for i, lab in enumerate(params_lab2):
                        params[lab] = params_set2[i]
               
                    params_order_fun = ['pi_A', 'pi_B', 'xi', 'lpub', 'lpriv', 'pa', 'pb']
                    params_nums = tuple([params[lab] for lab in params_order_fun]+[n])

                    curve_points, time_points = integrate_ode(vec_in, params_nums)
                    df_temp = dynamics_to_dataframe(n, curve_points, time_points, initial_vec_str=None)
                    values_y.append(df_temp['zA'].tolist()[-1])

                # plot a curve for each other_var
                if force_private:
                    line, = ax.plot(list_set2[0][cnt::len(other_var_values)], values_y[cnt::len(other_var_values)],
                        color=my_colors[cnt], alpha=0.8)
                else:
                    values_x = list(list_set2[0])
                    if cnt == 0:
                        x_values = values_x[cnt::len(other_var_values)] + values_x[-1:]
                        y_values = values_y[cnt::len(other_var_values)] + values_y[-1:]
                    elif cnt == len(other_var_values) - 1:
                        x_values = values_x[:1] + values_x[cnt::len(other_var_values)]
                        y_values = values_y[:1] + values_y[cnt::len(other_var_values)]
                    else:
                        x_values = values_x[:1] + values_x[cnt::len(other_var_values)] + values_x[-1:]
                        y_values = values_y[:1] + values_y[cnt::len(other_var_values)] + values_y[-1:]
                line, = ax.plot(x_values, y_values, color=my_colors[cnt],
                        marker=markers_lst[cnt], alpha=0.8, mec='black')
                   
                label=(r'${}={:.2f}'.format(
                    params2latex[other_var].strip('$'), other_var_val)).rstrip('0').rstrip('.') + r'$'
                lines.append(line)
                labels.append(label)
            if not force_private:
                if other_var == 'lpub':
                    excluded_keys = ['lpub', 'lpriv']
                else:
                    excluded_keys = [other_var]
                excluded_keys += ['n'] + params_lab2
                xlabel = r'%s'%params2latex[params_lab2[0]]
                #left, right, top, bottom, wspace, hspace
                subplots_params = [0.2, 0.8, 0.9, 0.1, 0.2, 0.2]
                subplots_params = [0.125, 0.9, 0.9, 0.1, 0.2, 0.2]
                save_fig(ax, prefix_name='za_vs_%s_lpubs' % var, xlabel=xlabel,
                    ylabel=r'$z_A$', xlims=[-0.01, 1.01], ylims=[-0.01, 1.01],
                    excluded_keys=excluded_keys, peer=peer,
                    lines_labels=[lines, labels],
                    subplots_params=subplots_params, params=params)

    excluded_keys = ['xi', 'pi_A', 'pi_B', 'lpub', 'lpriv', 'n'] + params_lab2
    xlabel = r'%s' % params2latex[params_lab2[0]]
    if force_private:
        save_fig(ax, prefix_name='private_za_vs_%s_lpubs' % var, xlabel=xlabel,
             ylabel=r'$z_A$', xlims=[-0.01, 1.01], ylims=[-0.01, 1.01],
             excluded_keys=excluded_keys, peer=peer, params=params)

   
def save_fig(ax, prefix_name='', xlabel=r'$t$', ylabel=r'$z_A$',
             xlims=[], ylims=[], com_legenda=True, lines_labels=[],
             subplots_params=[],
             params_order=['pi_A', 'pi_B', 'xi', 'pa', 'pb', 'lpub'],
             excluded_keys=[], peer=False, subplots=False,
             params=[]):

    str_params = gen_str_params(params, params_order=params_order, exclude_keys=excluded_keys)

    uniq_param = '_'.join(str_params.replace('\n', '').replace('$', '').split('; '))

    ax.grid(True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlims != []:
        ax.set_xlim(xlims)
    if ylims != []:
        ax.set_ylim(ylims)

    plt.tight_layout()

    if com_legenda:
        if subplots_params == []:
            plt.legend()
        else:
            lines, labels = lines_labels
            plt.legend(lines, labels, bbox_to_anchor=(0.5, -0.3), loc='center', ncol=2)
            # default params: 0.125, 0.9, 0.9, 0.1, 0.2, 0.2
            left, right, top, bottom, wspace, hspace = subplots_params
            plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom, wspace=wspace)

    path_to_fig_dir = 'figs/peer' if peer else 'figs'
    plt.savefig('%s/%s_%s.pdf' % (path_to_fig_dir, prefix_name, uniq_param), bbox_inches='tight')

    titulo = gen_str_params(params, params_order=params_order,
        symbols=params2latex, exclude_keys=excluded_keys).replace('\n', '')
    if subplots:
        ax.set_suptitle(titulo, fontsize=30)
    else:
        ax.set_title(titulo, fontsize=30)
    plt.savefig('%s/%s_%s.png' % (path_to_fig_dir, prefix_name, uniq_param), bbox_inches='tight')


if __name__ == '__main__':
    cli()
