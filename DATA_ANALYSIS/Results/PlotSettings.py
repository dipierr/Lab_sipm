# PlotSettings.py

import matplotlib.pyplot as plt

def PlotSettings():
    # Direct input
    plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    plt.rcParams['text.latex.preamble']=[r"\usepackage{siunitx}"]
    # Options, rcParams
    params = {  # LaTeX font
                'text.usetex' : True,
                'font.size'   : 12,
                'font.family' : 'serif',
                # Font Size in Plots
                'axes.titlesize' : 18,
                'axes.labelsize' : 18,
                'xtick.labelsize': 16,
                'ytick.labelsize': 16,
                'savefig.bbox'   : 'tight',
                # Plot Grid
                'axes.grid' : True,
                # errorbars
                'errorbar.capsize' : 3
              }
    # Update rcParams
    plt.rcParams.update(params)
