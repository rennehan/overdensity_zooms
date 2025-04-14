import numpy as np
import matplotlib.pyplot as plt
import pickle
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument("histogram_file")
parser.add_argument("number_of_bins")
parser.add_argument("redshift")
parser.add_argument("little_h")
parser.add_argument("box_size_cMpc_h")
parser.add_argument("save_file")
args = parser.parse_args()

histogram_file = args.histogram_file
number_of_bins = int(args.number_of_bins)
redshift = float(args.redshift)
scale_factor = 1.0 / (1.0 + redshift)
box_size = float(args.box_size_cMpc_h) / float(args.little_h)
save_file_fmt = "{:s}"

figsize = (6.4/1.5, 4.8/1.5)
dpi = 300
axes_fontsize = 14
lw = 3.5

def set_aspect(ax, log_x, log_y):
    if(log_x and log_y):
        ax.set_aspect(np.log10(ax.get_xlim()[1]/ax.get_xlim()[0])/np.log10(ax.get_ylim()[1]/ax.get_ylim()[0]))
    elif(log_x):
        ax.set_aspect(np.log10(ax.get_xlim()[1]/ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))
    elif(log_y):
        ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/np.log10(ax.get_ylim()[1]/ax.get_ylim()[0]))
    else:
        ax.set_aspect(np.diff(ax.get_xlim())/np.diff(ax.get_ylim()))

hist, edges = pickle.load(open(histogram_file, "rb"))

overdensity_hist, be = np.histogram(np.log10(1.0 + hist.flatten()), bins=len(edges)**3//number_of_bins)
scaled_be = be

plt.figure(figsize=figsize, dpi=dpi)
plt.title("Cell size L = {:.1f} cMpc".format((1.0 + redshift) * (edges[1] - edges[0])), fontsize=axes_fontsize)
plt.xlabel("log(1 + $\delta$)", fontsize=axes_fontsize)
plt.ylabel("log(# of cells)", fontsize=axes_fontsize)
plt.xlim([scaled_be.min(), scaled_be.max()])
plt.stairs(np.log10(overdensity_hist), scaled_be, color="blue", lw=lw, fill=True)
plt.axvline(0, ls=":", lw=lw, c="k", label="Mean")
plt.minorticks_on()
plt.grid(which="minor", alpha=0.3, ls=":")
plt.grid(which="major", alpha=0.7, ls=":")
plt.legend(fontsize=axes_fontsize - 8, loc="upper left",
           frameon=True, facecolor="w", framealpha=1.0)

set_aspect(plt.gca(), False, False)
plt.savefig(save_file_fmt.format(args.save_file), bbox_inches="tight")
plt.close()
