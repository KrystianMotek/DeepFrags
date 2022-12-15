import matplotlib.pyplot as plt


def alpha_histogram(values, file):
    fig = plt.gcf()
    axes = plt.axes()

    fig.set_size_inches(17.5, 7.0) # figure size

    # axes configuration
    axes.set_xlim(0, 180)
    axes.xaxis.set_label_coords(0.5007, -0.0787)
    axes.xaxis.set_major_locator(plt.MaxNLocator(19))
    axes.set_xlabel(r"$\alpha \hspace{0.25} (^\circ)$", fontsize=15)

    axes.hist(values, bins=30, density=True, fill=False) # create histogram
    plt.savefig(file)
    plt.clf()


def theta_histogram(values, file):
    fig = plt.gcf()
    axes = plt.axes()

    fig.set_size_inches(17.5, 7.0) # figure size

    # axes configuration
    axes.set_xlim(-180, 180)
    axes.xaxis.set_label_coords(0.5007, -0.0787)
    axes.xaxis.set_major_locator(plt.MaxNLocator(37))
    axes.set_xlabel(r"$\Theta \hspace{0.25} (^\circ)$", fontsize=15)

    axes.hist(values, bins=70, density=True, fill=False) # create histogram
    plt.savefig(file)
    plt.clf()


def correlation_plot(alpha, theta, file):
    fig = plt.gcf()
    axes = plt.axes()

    fig.set_size_inches(16.5, 16.5) # both axes are the same size

    # horizontal axis
    axes.set_xlim(0, 180)
    axes.set_xlabel(r"$\alpha \hspace{0.25} (^\circ)$", fontsize=25)
    axes.xaxis.set_major_locator(plt.MaxNLocator(18))
    axes.xaxis.set_label_coords(0.50, -0.04)

    # vertical axis 
    axes.set_ylim(-180, 180)
    axes.set_ylabel(r"$\Theta \hspace{0.25} (^\circ)$", fontsize=25)
    axes.yaxis.set_major_locator(plt.MaxNLocator(19))

    # create and save
    axes.scatter(x=alpha, y=theta, s=0.015)
    plt.savefig(file)
    plt.clf()
