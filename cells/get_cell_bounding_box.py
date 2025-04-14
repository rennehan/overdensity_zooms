import numpy as np
import pickle
import argparse as ap
import h5py
from astropy import units as u
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM

parser = ap.ArgumentParser()
parser.add_argument("data_dir")
parser.add_argument("num_files")
parser.add_argument("snapshot_idx")
parser.add_argument("file_type")
args = parser.parse_args()

data_dir = args.data_dir
data_file_fmt_multi = "{:s}/snapdir_{:03d}/snapshot_{:03d}.{:d}.hdf5"
data_file_fmt_single = "{:s}/snapshot_{:03d}.hdf5"
hist_file_fmt = "{:s}/overdensity_cells_{:03d}.pkl"
save_file_fmt = "{:s}/lowest_cell_pids_{:03d}.pkl"

file_type = args.file_type
num_files = int(args.num_files)
snapshot_idx = int(args.snapshot_idx)

hist_file = hist_file_fmt.format(data_dir, snapshot_idx)

data = pickle.load(open(hist_file, "rb"))
hist = data["hist3D_counts"]
edges = data["hist3D_edges"]

lowest_idx = np.argmin(hist)
ix, iy, iz = np.unravel_index(lowest_idx, hist.shape)
x_low, x_high = edges[0][ix], edges[0][ix + 1]
y_low, y_high = edges[1][iy], edges[1][iy + 1]
z_low, z_high = edges[2][iz], edges[2][iz + 1]

print("Cell location low: (x, y, z) = ({:g}, {:g}, {:g})".format(x_low, y_low, z_low))
print("Cell location high: (x, y, z) = ({:g}, {:g}, {:g})".format(x_high, y_high, z_high))

