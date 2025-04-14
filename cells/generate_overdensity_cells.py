import h5py
import pickle
import numpy as np
import seaborn as sb
import argparse as ap
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM

def set_aspect(ax, log_x, log_y):
    if(log_x and log_y):
        ax.set_aspect(np.log10(ax.get_xlim()[1]/ax.get_xlim()[0])/np.log10(ax.get_ylim()[1]/ax.get_ylim()[0]))
    elif(log_x):
        ax.set_aspect(np.log10(ax.get_xlim()[1]/ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))
    elif(log_y):
        ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/np.log10(ax.get_ylim()[1]/ax.get_ylim()[0]))
    else:
        ax.set_aspect(np.diff(ax.get_xlim())/np.diff(ax.get_ylim()))

plot_style = {"axes.axisbelow": True,
              "axes.edgecolor": "0",
              "axes.facecolor": "white",
              "axes.grid": False,
              "axes.labelcolor": ".15",
              "axes.spines.bottom": True,
              "axes.spines.left": True,
              "axes.spines.right": True,
              "axes.spines.top": True,
              "figure.facecolor": "white",
              "font.family": ["serif"],
              "font.sans-serif": ["Times New Roman"],
              "grid.color": ".8",
              "grid.linestyle": "-",
              "image.cmap": "rocket",
              "lines.solid_capstyle": "round",
              "patch.edgecolor": "w",
              "patch.force_edgecolor": True,
              "text.color": ".15",
              "xtick.bottom": True,
              "xtick.color": ".15",
              "xtick.direction": "in",
              "xtick.top": True,
              "ytick.color": ".15",
              "ytick.direction": "in",
              "ytick.left": True,
              "ytick.right": True}

axes_fontsize = 16
figsize = (6.4/1.5, 4.8/1.5)
dpi = 300
marker_size = 30
lw = 4.5

sb.set_style(plot_style)
colors = sb.color_palette("crest", 5)

parser = ap.ArgumentParser()
parser.add_argument("data_dir")
parser.add_argument("num_files")
parser.add_argument("snapshot_idx")
parser.add_argument("file_type")
args = parser.parse_args()

data_dir = args.data_dir
save_file_plot = "{:s}/overdensity_plot_{:03d}.pdf"
save_file_pkl = "{:s}/overdensity_cells_{:03d}.pkl"
data_file_fmt_multi = "{:s}/snapdir_{:03d}/snapshot_{:03d}.{:d}.hdf5"
data_file_fmt_single = "{:s}/snapshot_{:03d}.hdf5"

file_type = args.file_type
num_files = int(args.num_files)
snapshot_idx = int(args.snapshot_idx)
#file_type = "gadget"
#num_files = 8
#snapshot_idx = 111

if num_files == 1:
    data_file_fmt = data_file_fmt_single
else:
    data_file_fmt = data_file_fmt_multi

# open the first file to get the information about how big
# our array actually needs to be
if num_files == 1:
    header_file = data_file_fmt.format(data_dir, snapshot_idx)
else:
    header_file = data_file_fmt.format(data_dir, snapshot_idx, snapshot_idx, 0)
    
if file_type == "swift":
    with h5py.File(header_file, "r") as f:
        num_parts = f["Header"].attrs["NumPart_Total"][1]
        box_size = f["Header"].attrs["BoxSize"][0]
        h = f["Cosmology"].attrs["h"][0]
        mass = f["PartType1/Masses"][0] * 1.0e10
        omega_matter = f["Cosmology"].attrs["Omega_m"]
        omega_lambda = f["Cosmology"].attrs["Omega_lambda"]
        omega_baryon = f["Cosmology"].attrs["Omega_b"]
        redshift = f["Header"].attrs["Redshift"][0]
        length_to_Mpc = 1.0
        scale_factor = 1.0 / (1.0 + redshift)
else:
    with h5py.File(header_file, "r") as f:
        num_parts = f["Header"].attrs["NumPart_Total"][1]
        h = f["Parameters"].attrs["HubbleParam"]
        box_size = f["Parameters"].attrs["BoxSize"] / (1000.0 * h) # in cMpc units
        mass = f["Header"].attrs["MassTable"][1] * 1.0e10 / h
        omega_matter = f["Parameters"].attrs["Omega0"]
        omega_lambda = f["Parameters"].attrs["OmegaLambda"]
        redshift = f["Header"].attrs["Redshift"]
        length_to_Mpc = 1.0 / (1000.0 * h)
        scale_factor = 1.0 / (1.0 + redshift)

cosmo = FlatLambdaCDM(Om0=omega_matter, H0=100.*h)

print("The snapshot is at z = {:.3f}.".format(redshift))
print("There are {:d} particles total.".format(num_parts))
print("The box size is L = {:.3f} cMpc.".format(box_size))
print("The Hubble parameter is h = {:.3f}.".format(h))
print("The mass of the particles is M = {:g} Msun.".format(mass))

bins = np.arange(0, box_size * scale_factor, 50.0 * scale_factor)
if num_files == 1:
    print("Opening the single file {:s}".format(data_dir, snapshot_idx))
    with h5py.File(data_file_fmt.format(data_dir, snapshot_idx), "r") as f:
        num_parts_this_file = f["Header"].attrs["NumPart_ThisFile"][1]
        coords = np.array(f["PartType1/Coordinates"], dtype=np.float32) * scale_factor * length_to_Mpc
        hist, edges = np.histogramdd(coords, bins=[bins, bins, bins])
        del coords
else:
    hist = None
    edges = None
    print("Opening {:d} files and building the histogram.".format(num_files))
    for file_num in range(num_files):
        file_name = data_file_fmt.format(data_dir, snapshot_idx, snapshot_idx, file_num)
        print("Opening {:s}".format(file_name))
        with h5py.File(file_name, "r") as f:
            num_parts_this_file = int(f["Header"].attrs["NumPart_ThisFile"][1])
            print("There are {:d} particles in file #{:d}.".format(num_parts_this_file, file_num))
            coords = np.array(f["PartType1/Coordinates"], dtype=np.float32) * scale_factor * length_to_Mpc

            if hist is None:
                print("Generating the first histogram.")
                hist, edges = np.histogramdd(coords, bins=[bins, bins, bins])
                del coords
            else:
                print("Generating histogram #{:d}.".format(int(file_num + 1)))
                hist_to_add, _ = np.histogramdd(coords, bins=[bins, bins, bins])
                del coords
                hist += hist_to_add

print("Scaling the histogram to CGS density units.")
# Histogram should be the entire mass in the bin to compute
# the overdensity. Divide by the volume of each cell.
# Also, convert to g/cm**3. Divide by the mean and subtract one
# for overdensity!
hist *= mass / (bins[1] - bins[0])**3 * (1.989e33 / (3.086e24)**3)
print("Scaling histogram to the mean density.")
hist /= np.mean(hist)
print("Subtracting unity to get overdensity.")
hist -= 1.0

print("Saving pickle file.")
with open(save_file_pkl.format(data_dir, snapshot_idx), "wb") as f:
    pickle.dump({"hist3D_counts": hist,  # hist3d_counts are actually delta
                 "hist3D_edges": edges}, f, protocol=pickle.HIGHEST_PROTOCOL)
print("Histogram data saved to {:s}".format(save_file_pkl.format(data_dir, snapshot_idx)))

print("Creating histogram of overdensities.")
overdensity_hist, be = np.histogram(np.log10(1.0 + hist.flatten()), bins=len(bins)**3//25)
scaled_be = be

print("Producing and saving figure in {:s}".format(save_file_plot.format(data_dir, snapshot_idx)))
plt.figure(figsize=figsize, dpi=dpi)
plt.title("Cell size L = {:.1f} cMpc".format((1.0 + redshift) * (bins[1] - bins[0])), fontsize=axes_fontsize)
plt.xlabel("$\delta(\mathbf{x})$", fontsize=axes_fontsize)
plt.ylabel("# of cells", fontsize=axes_fontsize)
plt.xlim([scaled_be.min(), scaled_be.max()])
plt.stairs(overdensity_hist, scaled_be, color=colors[0], lw=lw, fill=True)
plt.axvline(0, ls=":", lw=lw, c="k", label="Mean")
plt.minorticks_on()
plt.grid(which="minor", alpha=0.3, ls=":")
plt.grid(which="major", alpha=0.7, ls=":")
plt.legend(fontsize=axes_fontsize - 8, loc="upper left",
           frameon=True, facecolor="w", framealpha=1.0)

set_aspect(plt.gca(), False, False) 
plt.savefig(save_file_plot.format(data_dir, snapshot_idx), bbox_inches="tight")
plt.close()

