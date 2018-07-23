import numpy as np
import pandas as pd


def read_params(params_file='params.txt', int_params=[]):
    # read parameters from file
    with open(params_file, "r") as f:
        filetext = f.readlines()

    # remove lines that start with #
    lines = [line for line in filetext if line[0] != '#']

    # remove newline char
    lines = [line.strip("\n").split() for line in lines if len(line.strip("\n"))]
    # convert str -> float for the numbers
    list_params = [[float(p) if "-" not in p else p[1:] for p in lista if "#" not in p]
        for lista in lines]
    # generates list of numbers for the lines started with '-'
    list_params = [internal_list if type(internal_list[0])==float
       else np.arange(internal_list[1], internal_list[2]+float(internal_list[3])/10,
             internal_list[3]) for internal_list in list_params]
    try:
        # make list of params names searching for the pattern '#<name>'
        symbols = [line[-1].split('#')[-1].strip() for line in lines]
    except Exception as e:
        print('no symbols found')
        return
    # create dictionary to store the values
    all_params = {}
    for i, sym in enumerate(symbols):
        if sym in int_params:
            list_params[i] = [int(elem) for elem in list_params[i]]
        all_params[sym] = list_params[i]

    if 'lpub' in all_params:
        all_params['lpriv'] = [1 - elem for elem in all_params['lpub']]
    return all_params


def gen_str_params(params, exclude_keys=[], params_order=[], symbols={}, add={}):
    """ Generates string with parameters
    """
    intern_copy_params = params.copy()
    for key in exclude_keys:
        intern_copy_params.pop(key, 0)

    params_order = [col for col in params_order if col not in exclude_keys]

    list_params = []
    list_keys = list(intern_copy_params.keys()) + list(add.keys())
    list_keys = [k for k in list_keys if k not in params_order]
    list_keys = params_order + list_keys

    count = 1

    num_breakl = 5 if len(list_keys) % 2 == 0 else 4

    for k in list_keys:
        name_param = symbols.get(k, r'$%s$' % k)
        # try to get the name from the dictionary
        # if there is no key k, just keep k
        v = intern_copy_params.get(k, add.get(k, -1))
        with_linebreak = '\n' if count % num_breakl == 0 else ''
        chunk_str = with_linebreak + (r'%s=%.5f' % (name_param, v))\
            .rstrip('0').rstrip('.')
        list_params.append(chunk_str)
        count += 1

    str_params = '; '.join(list_params)
    return str_params


def dynamics_to_dataframe(n, curve_points, time_points, initial_vec_str=None):
    types_names = [''.join(['z']+['a']*i+['b']*(n-i)) for i in range(n+1)]
    curve_points_transp = np.asarray(curve_points).T.tolist()
    df = pd.DataFrame()
    df['t'] = time_points
    for j, col in enumerate(types_names):
        if col == 'zaaa':
            continue
        df[col] = curve_points_transp[j]
    df['zaaa'] = 1 - df['zbbb'] - df['zaab'] - df['zabb']
    if initial_vec_str is not None:
        df['init'] = [initial_vec_str]*len(df)
    df['zA'] = df.apply(lambda row: row['zaaa']+2*row['zaab']/3 + row['zabb']/3, axis=1)
    return df
