import argparse, subprocess, shlex

parser = argparse.ArgumentParser(description = 'Perform bwa index analysis')
parser.add_argument("-f1", "--fastq1", help="fastq1", required=True)
parser.add_argument("-f2", "--fastq2", help="fastq2", required=True)
parser.add_argument("-pr", "--prefixout", help="prefixout", required=True)
parser.add_argument("-mi", "--min_len", help="min_len", required=True)
parser.add_argument("-p", "--threads", help="threads", required=True)
parser.add_argument("-a", "--adapters", help="adapters", required=True)

#Storing argument on variables
args = parser.parse_args()
fastq1 = args.fastq1
fastq2 = args.fastq2
prefixout = args.prefixout
min_len = args.min_len
threads = args.threads
adapters = args.adapters

try:
    fastp = (f"fastp -i {fastq1} -I {fastq2} -o {prefixout}.R1.fq.gz -O {prefixout}.R2.fq.gz --cut_front --cut_tail --qualified_quality_phred 20 -l {min_len} -h {prefixout}.quality.html --thread {threads} --adapter_fasta {adapters}")
    fastp = shlex.split(fastp)
    cmd_fastp = subprocess.Popen(fastp)
    cmd_fastp.wait()
except:
    print("Error in fastp step")