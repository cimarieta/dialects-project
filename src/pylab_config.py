import matplotlib.pyplot as plt


def setup():
    params = {
        'text.usetex': True,
        'font.family': 'serif',
        'font.serif': 'cm',
        'font.size': 16,
        'text.latex.unicode': True,
        'lines.linewidth': 2.,
        'axes.axisbelow': True
    }

    plt.rcParams.update(params)
    return plt
