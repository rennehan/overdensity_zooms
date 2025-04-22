import numpy as np
import h5py
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument("file_prefix")
parser.add_argument("start_index")
parser.add_argument("end_index")
args = parser.parse_args()

file_prefix = args.file_prefix
start_index = int(args.start_index)
end_index = int(args.end_index)

numpart_total = np.zeros(6, dtype=np.uint64)

print("Rewriting NumPart_ThisFile in all files.")
for i in range(start_index, end_index + 1):
    print("File index={:d}".format(i))

    with h5py.File("{:s}.{:d}.hdf5".format(file_prefix, i), "a") as f:
        pids = np.array(f["PartType1/ParticleIDs"])
        count_1 = len(pids)

        pids = np.array(f["PartType2/ParticleIDs"])
        count_2 = len(pids)

        numpart_thisfile = np.array([0, count_1, count_2, 0, 0, 0], dtype=np.uint64)
        print("NumPart_ThisFile=", numpart_thisfile)
        numpart_total += numpart_thisfile
        f["Header"].attrs["NumPart_ThisFile"] = numpart_thisfile

print("Rewriting NumPart_Total in all files.")
for i in range(start_index, end_index + 1):
    print("File index={:d}".format(i))
    print("NumPart_Total=", numpart_total)

    with h5py.File("{:s}.{:d}.hdf5".format(file_prefix, i), "a") as f:
        f["Header"].attrs["NumPart_Total"] = numpart_total

