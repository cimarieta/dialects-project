#!/bin/bash

mkdir -p src/peer
cd src
python general_plot.py  bar_plot --peer
python general_plot.py  plot_za_vs_var --peer --var pa --other_var xi
python general_plot.py plot_several_lpubs --peer
python pa_vs_lpub_critical.py

python general_plot.py  bar_plot --unanimity
python general_plot.py  plot_za_vs_var --unanimity --var pa --other_var xi
python plot_graphs.py
python general_plot.py plot_several_lpubs --unanimity
python general_plot.py plot_several_lpubs_zis --unanimity
python plot_arrow_equilibria.py
python plot_arrow_equilibria.py --lpub 0.8
