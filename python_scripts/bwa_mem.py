import argparse, subprocess, shlex, os

parser = argparse.ArgumentParser(description = 'Perform bwa index analysis')
parser.add_argument("-f", "--fasta", help="fasta_file", required=True)
parser.add_argument("-pr", "--prefixout", help="prefixout", required=True)
parser.add_argument("-p", "--threads", help="threads", required=True)

#Storing argument on variables
args = parser.parse_args()
fasta_file = args.fasta
prefixout = args.prefixout
threads = args.threads

try:
    bwa_mem = (f"bwa mem -t {threads} {fasta_file} {prefixout}.R1.fq.gz {prefixout}.R2.fq.gz -o {prefixout}.bam")
    bwa_mem = shlex.split(bwa_mem)
    cmd_bwa_mem = subprocess.Popen(bwa_mem)
    cmd_bwa_mem.wait()
    sam_sort = (f"samtools sort -o {prefixout}.sorted.bam {prefixout}.bam")
    sam_sort = shlex.split(sam_sort)
    cmd_sam_sort = subprocess.Popen(sam_sort)
    cmd_sam_sort.wait()
    sam_index = (f"samtools index {prefixout}.sorted.bam")
    sam_index = shlex.split(sam_index)
    cmd_sam_index = subprocess.Popen(sam_index)
    cmd_sam_index.wait()
except:
    print("Error in bwa_mem step")

os.remove(prefixout+'.bam')