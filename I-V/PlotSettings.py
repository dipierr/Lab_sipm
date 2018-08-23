# PlotSettings.py

import matplotlib.pyplot as plt

def PlotSettings():
    # Direct input
    plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    # Options, rcParams
    params = {  # LaTeX font
                'text.usetex' : True,
                'text.latex.unicode': True,
                'font.size' : 12,
                # Font Size in Plots
                'axes.titlesize': 18,
                'axes.labelsize': 18,
                'xtick.labelsize': 16,
                'ytick.labelsize': 16,
                'savefig.bbox': 'tight'
              }
    # Update rcParams
    plt.rcParams.update(params)
