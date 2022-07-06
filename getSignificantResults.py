# -----------------
#  getSignificantResults.py
#  -----------------
#
# Author: Maira R. Rodrigues   (Lab. Bioest. e Biol. Comp.)
# Contributor(s):
#
# Get significant (skato-value <= 1x10-5) results 
# for SKAT output files in input folder ending with .skat
#
## Updated at/by: 03/12/2020
#
# Command:
#  python ./getSignificantResults.py -t 0.00001 -c [Rare|Ultrarare] -d <path>
# Sample cmd:
# python ./getSignificantResults.py -t 0.00001 -c Rare -d D:\Project-EpimiRNA\\skat_robust\\rare

#!/usr/bin/python
# coding: utf-8
from os import listdir
from os.path import isfile, join
from argparse import ArgumentParser

# INPUTS: 
parser = ArgumentParser()
parser.add_argument("-t", "--thr", help="P-value threshold", type=float,required=True)
parser.add_argument("-c", "--cla", help="Variant class [Rare,Ultrarare", required=True)
parser.add_argument("-d", "--dir", help="Directory with SKATO association results.", required=True)

args = parser.parse_args()

# Set threshold by total number of genes tested
#ngenes =  5446
#thr = 0.05/ngenes 
# Set threshold by total number of genes teste (Matt)
thr      = args.thr
varclass = args.cla
resultdir= args.dir 

# VARIABLES 
DEBUG=False  # [True,False]

outputfile = "significant_results_SKATO_"+varclass+".txt"

# script printed header
print( "# ----------------------------------------------")
print( "# Showing significant results for SKATO analyses")
print( "# ----------------------------------------------\n")

print("Significance threshold: "+str(thr)+"\n")

# OUTPUTS

### Functions
########################################################################
def debug(mystring):
    if (DEBUG):
        print(mystring)


def read_dir(label,mypath):
#    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
#    for file in onlyfiles:
#        debug(file)
    # print table header
#   reads only analyses with all genes in each phenotype and all phenotypes together. 
    header="Variant_Class\tCohort\tPhenotype\tGene\tSKATO_gene_pvalue\tN_SNPs"
    print(header+"\n")
    # open for writing output table  
    with open(join(mypath, outputfile), 'w') as o:
        o.write("#significance threshold: "+str(thr)+"\n")
        o.write(header+"\n")
        for file in listdir(mypath):
            if (file.endswith(".skat.adj")):
                debug(file+"\n")
                debug(join(mypath, file))
                # read file
                with open(join(mypath, file), 'r') as f:
                    next(f) # SKIPS FIRST LINE
                    for line in f:
                        # REMOVE NEWLINE 
                        if line.count('\n') > 0:
                            line = line.strip('\n')
                        if line.count('\r') > 0:
                            line = line.strip('\r')
                        # SPACE-SEPARATED COLUMNS
                        #gene	SKATO	linear.weighted	IBS	IBS.weighted	quadratic	2wayIX	nsnps	snps
                        line  = line.split('\t')
                        gene  = line[0]
                        skato = float(line[1])
                        nsnps = line[7]
                        if (skato<=thr):
                            parts = file.split('.') 
                            cohort = parts[0]
                            if (parts[-1] == "autosomes"):
                                phenotype="All"
                            else:
                                phenotype=parts[-1]
                            print(  label+"\t"+cohort+"\t"+phenotype+"\t"+gene+"\t"+str(skato)+"\t"+nsnps)
                            o.write(label+"\t"+cohort+"\t"+phenotype+"\t"+gene+"\t"+str(skato)+"\t"+nsnps+"\n")

### MAIN CODE 
##########################################################################

read_dir(varclass,resultdir)

print("\n>> Table also saved in: "+outputfile)

print("\n")
