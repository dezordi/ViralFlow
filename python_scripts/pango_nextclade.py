import argparse, os

parser = argparse.ArgumentParser(description = 'Perform bwa index analysis')
parser.add_argument("-pr", "--prefixout", help="prefixout", required=True)
parser.add_argument("-dp", "--depth", help="depth", required=True)
parser.add_argument("-p", "--threads", help="threads", required=True)

#Storing argument on variables
args = parser.parse_args()
prefixout = args.prefixout
depth = args.depth
threads = args.threads

try:
    num_lines = sum(1 for line in open(f"{prefixout}.depth{depth}.fa.bc.intrahost.short.tsv"))
    if num_lines > 1:
        os.system(f"cat {prefixout}.depth{depth}.fa {prefixout}.depth{depth}.fa.algn.minor.fa > {prefixout}.depth{depth}.all.fa \
            && nextclade -i {prefixout}.depth{depth}.all.fa -c {prefixout}.depth{depth}.all.fa.nextclade.csv --jobs {threads} \
            && pangolin {prefixout}.depth{depth}.all.fa -t {threads} --outfile {prefixout}.depth{depth}.all.fa.pango.csv")
    else:
        os.system(f"nextclade -i {prefixout}.depth{depth}.fa -c {prefixout}.depth{depth}.fa.nextclade.csv --jobs {threads} \
            && pangolin {prefixout}.depth{depth}.fa -t {threads} --outfile {prefixout}.depth{depth}.fa.pango.csv")
except:
    print("Error in pangolin and nextclade step")