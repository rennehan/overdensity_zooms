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

print("Rewriting Masses for PartType1 in all files.")
for i in range(start_index, end_index + 1):
    print("File index={:d}".format(i))

    with h5py.File("{:s}.{:d}.hdf5".format(file_prefix, i), "a") as f:
        mass_1 = np.array(f["Header"].attrs["MassTable"])[1]
        count_1 = np.array(f["Header"].attrs["NumPart_ThisFile"], dtype=np.uint64)[1]
        print("Create dataset with masses all equal to {:.3f}".format(mass_1))

        try:
            p = f["PartType1"]
            p.create_dataset("Masses", data=np.ones(count_1)*mass_1, dtype=np.float32)

            print("Reset the mass table to zeros")        
            mass_table = np.zeros(6, dtype=np.float32)
            f["Header"].attrs["MassTable"] = mass_table
        except:
            print("Failed in index={:d}".format(i))
            pass


