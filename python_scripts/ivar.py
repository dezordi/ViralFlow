import argparse, os

parser = argparse.ArgumentParser(description = 'Perform bwa index analysis')
parser.add_argument("-f", "--fasta", help="fasta_file", required=True)
parser.add_argument("-pr", "--prefixout", help="prefixout", required=True)
parser.add_argument("-dp", "--depth", help="depth", required=True)

#Storing argument on variables
args = parser.parse_args()
fasta_file = args.fasta
prefixout = args.prefixout
depth = args.depth

try:
    os.system(f"samtools mpileup -d 50000 --reference {fasta_file} -a -B {prefixout}.sorted.bam | ivar variants -p {prefixout} -q 30 -t 0.05 \
                && samtools mpileup -d 50000 --reference {fasta_file} -a -B {prefixout}.sorted.bam | ivar consensus -p {prefixout} -q 30 -t 0 -m {depth} -n N \
                && samtools mpileup -d 50000 --reference {fasta_file} -a -B {prefixout}.sorted.bam | ivar consensus -p {prefixout}.ivar060 -q 30 -t 0.60 -m {depth} -n N \
                && mv {prefixout}.fa {prefixout}.depth{depth}.fa \
                && sed -i -e 's/>.*/>'{prefixout}'/g' {prefixout}.depth{depth}.fa \
                && sed -i -e 's/__/\//g' -e 's/--/|/g' {prefixout}.depth{depth}.fa \
                && mv {prefixout}.ivar060.fa {prefixout}.depth{depth}.amb.fa \
                && sed -i -e 's/>.*/>'{prefixout}'/g' {prefixout}.depth{depth}.amb.fa \
                && sed -i -e 's/__/\//g' -e 's/--/|/g' {prefixout}.depth{depth}.amb.fa")
except:
    print("Error in ivar step")