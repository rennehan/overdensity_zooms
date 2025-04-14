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
#data_dir = "/mnt/home/drennehan/ceph/simulating_the_universe/overdensity_intersections/sims/4096_800/output"
data_file_fmt_multi = "{:s}/snapdir_{:03d}/snapshot_{:03d}.{:d}.hdf5"
data_file_fmt_single = "{:s}/snapshot_{:03d}.hdf5"
hist_file_fmt = "{:s}/overdensity_cells_{:03d}.pkl"
save_file_fmt = "{:s}/lowest_cell_pids_{:03d}.pkl"

file_type = args.file_type
num_files = int(args.num_files)
snapshot_idx = int(args.snapshot_idx)
#file_type = "gadget"
#num_files = 8
#snapshot_idx = 111

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

if num_files == 1:
    print("Opening the single file {:s}".format(data_dir, snapshot_idx))
    with h5py.File(data_file_fmt.format(data_dir, snapshot_idx), "r") as f:
        num_parts_this_file = f["Header"].attrs["NumPart_ThisFile"][1]
        particle_ids = np.array(f["PartType1/ParticleIDs"], dtype=np.uint64)
        coords = np.array(f["PartType1/Coordinates"], dtype=np.float32) * scale_factor * length_to_Mpc
        del coords
else:
    particle_ids = np.zeros(num_parts, dtype=np.uint64)
    coords = np.zeros((num_parts, 3), dtype=np.float32)
    print("Opening {:d} files and building the arrays.".format(num_files))
    start_idx = 0
    for file_num in range(num_files):
        file_name = data_file_fmt.format(data_dir, snapshot_idx, snapshot_idx, file_num)
        print("Opening {:s}".format(file_name))
        with h5py.File(file_name, "r") as f:
            num_parts_this_file = int(f["Header"].attrs["NumPart_ThisFile"][1])
            end_idx = start_idx + num_parts_this_file
            print("There are {:d} particles in file #{:d}.".format(num_parts_this_file, file_num))
            particle_ids[start_idx:end_idx] = np.array(f["PartType1/ParticleIDs"], dtype=np.uint64)
            coords[start_idx:end_idx] = np.array(f["PartType1/Coordinates"], dtype=np.float32) * scale_factor * length_to_Mpc
            start_idx = end_idx

x_cond = (coords[:, 0] > x_low) & (coords[:, 0] < x_high)
y_cond = (coords[:, 1] > y_low) & (coords[:, 1] < y_high)
z_cond = (coords[:, 2] > z_low) & (coords[:, 2] < z_high)
del coords

select_idx = x_cond & y_cond & z_cond

with open(save_file_fmt.format(data_dir, snapshot_idx), "wb") as f:
    pickle.dump(particle_ids[select_idx], f, protocol=pickle.HIGHEST_PROTOCOL)


