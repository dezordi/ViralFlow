import argparse, os, subprocess, shlex

parser = argparse.ArgumentParser(description = 'Perform bwa index analysis')
parser.add_argument("-f", "--fasta", help="fasta_file", required=True)
parser.add_argument("-pr", "--prefixout", help="prefixout", required=True)
parser.add_argument("-dp", "--depth", help="depth", required=True)
parser.add_argument("-p", "--threads", help="threads", required=True)

#Storing argument on variables
args = parser.parse_args()
fasta_file = args.fasta
prefixout = args.prefixout
depth = args.depth
threads = args.threads

try:
    with open(f"{prefixout}.depth{depth}.fa.bc","w") as bamread_output: #abm = annotation bed merged
        bwa_read = (f"bam-readcount -d 50000 -b 30 -q 30 -w 0 -f {fasta_file} {prefixout}.sorted.bam")
        bwa_read = shlex.split(bwa_read)
        cmd_bwa_read = subprocess.Popen(bwa_read, stdout=bamread_output)
        cmd_bwa_read.wait()
    with open(f"{prefixout}.depth{depth}.fa.algn","w") as mafft_output:
        mafft = (f"mafft --thread {threads} --keeplength --add {prefixout}.depth{depth}.fa {fasta_file}")
        mafft = shlex.split(mafft)
        cmd_mafft = subprocess.Popen(mafft,stdout=mafft_output)
        cmd_mafft.wait()
    os.system(f"python /home/python_scripts/intrahost.py -in {prefixout}.depth{depth}.fa.bc -al {prefixout}.depth{depth}.fa.algn")
except:
    print("Error in minor variants step")
