import numpy as np
import h5py
import argparse as ap
from gizmorw import read_gizmo_file
import pickle

parser = ap.ArgumentParser()
parser.add_argument("data_dir")
parser.add_argument("num_files")
parser.add_argument("snapshot")
parser.add_argument("ics_file")
parser.add_argument("index_list_file")
parser.add_argument("box_size_in_kpc")
parser.add_argument("ic_file_idx")
args = parser.parse_args()

num_files = int(args.num_files)
snapshot = int(args.snapshot)
ICs_file = args.ics_file
box_size_in_kpc = float(args.box_size_in_kpc)
ic_file_idx = int(args.ic_file_idx)

if box_size_in_kpc <= 0.:
    raise ValueError("box_size_in_kpc must be >0")

save_file_fmt = "{:s}/lowest_coordinates_{:03d}.{:d}.dat"

print("Reading the IC data {:s}".format(ICs_file))

with h5py.File(ICs_file, "r") as f:
    particle_ids = np.array(f["PartType1/ParticleIDs"], dtype=np.uint64)
    coordinates = np.array(f["PartType1/Coordinates"], dtype=np.float32)

print("Reading pickle file {:s}".format(args.index_list_file))
pids_region = pickle.load(open(args.index_list_file, "rb"))

print("Finding particle IDs in the ICs.")
particle_idx = np.isin(particle_ids, pids_region, assume_unique=True)

if not np.any(particle_idx):
    with open(save_file_fmt.format(args.data_dir, snapshot, ic_file_idx), "w") as f:
        f.write("\n")

    raise ValueError("No matches in this snapshot")

print("Grabbing and normalize coordinates to the box size.")
coords_list = coordinates[particle_idx] / box_size_in_kpc  # Normalize to box coords
del particle_ids, coordinates, particle_idx, pids_region

print("Bounding box (in L=1 units): ")
print("Min. x: %g" % np.amin(coords_list[:, 0]))
print("Max. x: %g" % np.amax(coords_list[:, 0]))

print("Min. y: %g" % np.amin(coords_list[:, 1]))
print("Max. y: %g" % np.amax(coords_list[:, 1]))

print("Min. z: %g" % np.amin(coords_list[:, 2]))
print("Max. z: %g" % np.amax(coords_list[:, 2]))

out_of_bounds_idx = np.where((coords_list < 0) | (coords_list > 1))[0]

print("There are %d particles out of bounds." % len(out_of_bounds_idx))

if len(out_of_bounds_idx) > 0:
    # Shift everything to [0, 1]
    for i in range(0, 3):
        neg_idx = np.where(coords_list[:, i] < 0)
        pos_idx = np.where(coords_list[:, i] > 1)

        coords_list[neg_idx, i] += 1.0
        coords_list[pos_idx, i] -= 1.0

with open(save_file_fmt.format(args.data_dir, snapshot, ic_file_idx), "w") as f:
    for i in range(len(coords_list)):
        f.write("%g %g %g\n" % (coords_list[i, 0], coords_list[i, 1], coords_list[i, 2]))


