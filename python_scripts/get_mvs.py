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
    os.system(f"python /home/IAM_SARSCOV2/minor_finder.py -in {prefixout}.depth{depth}.fa.bc \
                && sed -i -e 's/__/\//g' -e 's/--/|/g' {prefixout}.depth{depth}.fa.bc.fmt.minors.tsv \
                && python /home/IAM_SARSCOV2/major_minor.py -in {prefixout}.depth{depth}.fa.bc.fmt.minors.tsv")
    with open(f"{prefixout}.depth{depth}.fa.algn","w") as mafft_output:
        mafft = (f"mafft --thread {threads} --keeplength --add {prefixout}.depth{depth}.fa {fasta_file}")
        mafft = shlex.split(mafft)
        cmd_mafft = subprocess.Popen(mafft,stdout=mafft_output)
        cmd_mafft.wait()
    os.system(f"python /home/IAM_SARSCOV2/put_minor.py -in {prefixout}.depth{depth}.fa.algn -mv {prefixout}.depth{depth}.fa.bc.fmt.minors.tsv.fmt \
            && mv {prefixout}.depth{depth}.fa.algn.minor.fa {prefixout}.depth{depth}.minor.fa")
except:
    print("Error in minor variants step")