import argparse, subprocess, shlex

parser = argparse.ArgumentParser(description = 'Perform bwa index analysis')
parser.add_argument("-in", "--input", help="Input file", required=True)

#Storing argument on variables
args = parser.parse_args()
input_file = args.input

try:
    bwa_index = (f"bwa index {input_file}")
    bwa_index = shlex.split(bwa_index)
    cmd_bwa_index = subprocess.Popen(bwa_index)
    cmd_bwa_index.wait()
except:
    print("Error in bwa index step")