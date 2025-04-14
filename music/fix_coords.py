import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument("file_name")
args = parser.parse_args()

# Open the coordinates file and remove any lines
# that only have newline characters. This can
# happen if you concatenate some empty files together
# since the parent volume will only have some
# IC files contribute to the zoom region.
with open(args.file_name, "r") as f:
    lines = f.readlines()

non_empty_lines = [line for line in lines if line.strip()]

with open(args.file_name, "w") as f:
    f.writelines(non_empty_lines)

