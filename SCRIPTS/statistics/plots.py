import numpy as np
import matplotlib.pyplot as plt


def format_axis(x, pos):
    # change ticks to radians 
    # this function is used by class matplotlib.ticker.FuncFormatter 

    n = int(np.round(x * 2 / np.pi)) # half of pi multiples

    if n == 0: return "$0$"
    if n == 1: return "$\pi/2$"
    if n == 2: return "$\pi$"
    if n == -1: return "$-\pi/2$"
    if n == -2: return "-$\pi$"
    else: return None


def theta_histogram(values, file):
    # return distribution of dihedral angles given in radians
    axes = plt.axes()
    axes.xaxis.set_major_formatter(plt.FuncFormatter(format_axis))
    axes.set_xlim(-np.pi, np.pi)
    axes.set_xlabel("$\Theta$")
    axes.set_yticks([])
    plt.hist(values, alpha=0.45)
    plt.savefig(file)
    plt.clf()
