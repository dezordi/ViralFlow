import argparse, os

parser = argparse.ArgumentParser(description = 'Perform bwa index analysis')
parser.add_argument("-pr", "--prefixout", help="prefixout", required=True)

#Storing argument on variables
args = parser.parse_args()
prefixout = args.prefixout

try:
    os.system(f"bedtools bamtobed -i {prefixout}.sorted.bam > {prefixout}.sorted.bed \
                    && samtools view {prefixout}.sorted.bam -u | bamdst -p {prefixout}.sorted.bed -o . \
                    && gunzip ./region.tsv.gz \
                    && gunzip ./depth.tsv.gz \
                    && sed -i -e 's/NC_045512\.2/'{prefixout}'/g' chromosomes.report \
                    && sed -i -e 's/__/\//g' -e 's/--/|/g' chromosomes.report")
except:
    print("Error in assembly metrics step")