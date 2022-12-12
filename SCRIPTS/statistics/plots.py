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
