# -----------------
#  ped2scikitlearn.py
#  -----------------
#
# Author: Maira R. Rodrigues
# Contributor(s):
# Updated at/by: 21/01/2020 (Maira)
#
# Description: Prepare input for scikit-learn library. 
# Input: 
#	RAW file with recodeA option (encoded for additive model 0,1,2)
#	FAM
#   BIM
#   FRQ.CC file with Case and Control frequencies (calculated by plink)
#   Missing code [0,3]
#   PCA evec file (optional)
# Requirements:
#
# Output: 
#	.target file  (dimension  [1,samples])
#	.data    file (dimensions [samples,features])
#	.weights file (dimensions [3,features]) (columns SNP_INDEX MAF WEIGHT) * Here weight is calculated as the difference in MAF of cases and controls.
#   .covar   file (dimensions [samples,covars]) * gender + pca
#   .samples file (dimensions [samples,1]) SAMPLE_IDS SAME AS PED_FILE
#
# Command:
#  python .\ped2scikitlearn.py -m 0 -r ..\ped\combGVFS.g.clean2.raw -f ..\ped\combGVFS.g.clean2.fam -b ..\ped\combGVFS.g.clean2.bim -q ..\ped\combGVFS.g.clean2.frq.cc -p ..\pca\combGVFS.g.clean1a.evec


#!/usr/bin/python
# coding: utf-8
from argparse import ArgumentParser
import sys
import subprocess
from math import *
import time

start = time.time();

# INPUTS: 
parser = ArgumentParser()
parser.add_argument("-r", "--raw", help="File .raw", required=True)
parser.add_argument("-m", "--missing", help="Numerical code for missing value [0 or 3]",type=int, required=True)
parser.add_argument("-f", "--fam", help="Full path and name of FAM file", required=True)
parser.add_argument("-b", "--bim", help="Full path and name of BIM file", required=True)
parser.add_argument("-q", "--frq", help="Full path and name of frequency .frq.cc file", required=True)
parser.add_argument("-p", "--pca", help="Full path and name of smartpca .evec file", required=False)
parser.add_argument("-k", "--plinkpca", help="Full path and name of plink .eigenvec file", required=False)
#MISSING = '3' # ENCODING FOR MISSING DATA
#MISSING = '0' # ENCODING FOR MISSING DATA EQUAL TO MOST COMMON GENOTYPE.

args = parser.parse_args()

rawfile= args.raw
famfile= args.fam
bimfile= args.bim
frqfile= args.frq
MISSING = str(args.missing)

#nsnps = 105229

# FUNCTION DEFINITIONS
def openlog(filename):
	f = open(filename, 'w')	
	return f;
  
def printme(filehandler,str):
	filehandler.write(str+'\n')
	return;

def closelog(filehandler):
	filehandler.close();
	return;

def countsnps(filename):
    count=0
    with open(filename, 'r') as f:
        # READ OUT HEADER
        #next(f)
        for line in f:
            count+=1
                
    return count;

####################	
##### MAIN CODE ####
####################
print("# -------------------")
print("# ped2scikitlearn.py")
print("# -------------------\n")

classindex = 5  # column number for class and gender in fam file
genderindex= 4

nsnps = countsnps(bimfile)
print("> BIM has "+str(nsnps)+" SNPs")

# FOR PCA EIGENVALUES AS COVARIATES
if(str(args.pca)!="None"):
    totalpcs = 10
    eigen = {}  #key [id][pc] - value 
    evecflag = True
else:
    evecflag = False 


if(evecflag):
	evecfile=args.pca
	with open(evecfile, 'r') as f:
		next(f) # SKIPS FIRST LINE
		for line in f:
			# REMOVE NEWLINE 
			if line.count('\n') > 0:
				line = line.strip('\n')
			if line.count('\r') > 0:
				line = line.strip('\r')
			# SPACE-SEPARATED COLUMNS	
			line = line.split()
			iid = line[0]
			for pc in range(1,totalpcs+1):
				eigen[(iid,pc)] = line[pc]


# FOR PLINK PCA EIGENVALUES AS COVARIATES
if(str(args.plinkpca)!="None"):
    totalpcs = 10
    eigen = {}  #key [id][pc] - value 
    pevecflag = True
else:
    pevecflag = False 


if(pevecflag):
	evecfile=args.plinkpca
	with open(evecfile, 'r') as f:
		#next(f) # SKIPS FIRST LINE
		for line in f:
			# REMOVE NEWLINE 
			if line.count('\n') > 0:
				line = line.strip('\n')
			if line.count('\r') > 0:
				line = line.strip('\r')
			# SPACE-SEPARATED COLUMNS	
			line = line.split()
			fid = line[0]
			iid = line[1]
			for pc in range(2,totalpcs+2):
				eigen[(iid,(pc-1))] = line[pc]

# OUTPUTS:
data_filename   = rawfile[:-4] + ".data"      # ALL FILES ARE SPACE-SEPARATED
class_filename  = rawfile[:-4] + ".target"
weights_filename= rawfile[:-4] + ".weights"     
covariates_filename= rawfile[:-4] + ".covar"     
samples_filename= rawfile[:-4] + ".samples"     

# VARIABLES
idclassmap = {}
gender = {}
classcount = [0]*2 # Binary

totallables = 0 
totalsamples = 0
totalcovarsamples = 0
# mylist = [[0 for x in range(columns)] for x in range(rows)]
#frequency_list = [[0 for x in range(2)] for y in range(nsnps)] #coluna,linha  

## Codes in fam file:
## 1 control ; 2 affected
## Sex code ('1' = male, '2' = female, '0' = unknown)
##

## Codes for Target file:
## 0=control/negative  1=case/positive
#READING CLASS FILE TO MAP ID->CLASS
with open(famfile, 'r') as e:
	# READ OUT HEADER
	#next(e)
	for line in e:
		# REMOVE NEWLINE 
		if line.count('\n') > 0:
			line = line.strip('\n')
		if line.count('\r') > 0:
			line = line.strip('\r')
		# FAM FILE HAS SPACE-SEPARATED COLUMNS	
		line = line.split(' ')
		iid = line[1]
		cls = line[classindex]
		sex = line[genderindex]
		if iid not in idclassmap:
			
			if cls == "2":
				idclassmap[iid] = 1
			if cls == "1":
				idclassmap[iid] = 0
			gender[iid] = sex 
			totallables+=1
			

# COVARS FILE HEADER (IN ORDER)
# READING RECODE raw FILE - SPACE SEPARATED
print("Reading RAW file... (it can take a while)")
with open(samples_filename, 'w') as samples:
	with open(covariates_filename, 'w') as covar:
		with open(data_filename, 'w') as data:
			with open(class_filename, 'w') as labels:
				with open(rawfile, 'r') as f:
					next(f) # SKIPS FIRST LINE
					for line in f:
						# REMOVE NEWLINE 
						if line.count('\n') > 0:
							line = line.strip('\n')
						if line.count('\r') > 0:
							line = line.strip('\r')
						# SPACE-SEPARATED COLUMNS	
						line = line.split(' ')
						iid = line[1]
						sex = line[4] # Sex code ('1' = male, '2' = female, '0' = unknown)
						if iid in idclassmap:
							# ITERATE TIL END OF LINE AND WRITE IN .DATA FILE
							for i in range(6,len(line)):
								# CHANGE MISSING DATA ENCODING
								if ((line[i] != '0') and (line[i] != '1') and (line[i] != '2')):	
									genotype = MISSING
								else:
									# CALCULATE SNP ALLELE FREQUENCIES FROM GENOTYPE
									# CONSIDERING COLUMN INDICES AS: MINOR=0 AND MAJOR=1 
									#if (line[i] == '0'): #0 MINOR AND 2 MAJOR ALLELES
									#	frequency_list[i-6][1]+=2
									#if (line[i] == '1'): #1 MINOR AND 1 MAJOR ALLELES
									#	frequency_list[i-6][0]+=1
									#	frequency_list[i-6][1]+=1
									#if (line[i] == '2'): #2 MINOR AND 0 MAJOR ALLELES
									#	frequency_list[i-6][0]+=2
									# ENCODING IS THE SAME
									genotype = line[i]
								data.write(genotype+' ');
							data.write('\n')
							totalsamples+=1
							# WRITE .LABELS FILE
							labels.write(str(idclassmap[iid])+' ')

							# WRITE .SAMPLES FILE
							samples.write(iid+'\n')
	
							# WRITE COVARS FILE 1=male 2=female
							covar.write(sex)
							if(evecflag or pevecflag):
								for pc in range(1,totalpcs+1):
									covar.write('\t'+eigen[(iid,pc)])
							covar.write('\n')
							totalcovarsamples+=1
							
							# COUNT SAMPLES PER CLASS
							classcount[idclassmap[iid]]+=1
						else:
							#sys.exit("ERROR: PED sample "+iid+" does NOT have a matching class in labels file.")
							print("ATTENTION: PED sample "+iid+" does NOT have a matching ID in FAM file.\n")

# CALCULATE WEIGHTS AND WRITE TO FILE
# need to open frequencys file, read and store. 
snp_index = 0
with open(weights_filename, 'w') as weights:
    with open(frqfile, 'r') as e:
        # READ OUT HEADER
        next(e)
        for line in e:
            # REMOVE NEWLINE 
            if line.count('\n') > 0:
                line = line.strip('\n')
            if line.count('\r') > 0:
                line = line.strip('\r')
            # Freq FILE HAS SPACE-SEPARATED COLUMNS	
            #  CHR             SNP   A1   A2        MAF_A        MAF_U  NCHROBS_A  NCHROBS_U
            # MAF_A	Allele 1 frequency in cases
            # MAF_U	Allele 1 frequency in controls
            line = line.split()
            chr = line[0]
            snp = line[1]
            
            maf_cases    = float(line[4])
            maf_controls = float(line[5])

            maf_dif = abs(maf_cases-maf_controls) 
            
            #print str(i)+" MINOR count:"+str(frequency_list[i][0])+" MAJOR count:"+str(frequency_list[i][1])+" sum:"+str(sum(frequency_list[i]))+" MAF:"+str(maf)+"\n"
            #if (maf==0):
            #    w = 0.000
            #else:
            #    w = (1.0)/sqrt(maf*(1-maf))
            output = str(snp_index) + "\t" + str(maf_cases)[0:5] + "\t" + str(maf_dif)[0:5] + "\n"
            weights.write(output)
            snp_index=snp_index+1


print( "Finished.\n")
print("#Class 0 (controls): "+str(classcount[0]))
print("#Class 1 (cases): "+str(classcount[1]))

print( "#Total samples with labels: "+str(totallables))	
print( "#Total samples: "+str(totalsamples))
print( "#Total covar samples: "+str(totalcovarsamples))	
print( "#Total number of snps: "+str(nsnps))

print( "Output:" )
print( ""+class_filename)
print( "	.target file  (dimension  [1,samples])")
print( ""+data_filename)
print( "	.data    file (dimensions [samples,features])")
print( ""+weights_filename)
print( "	.weights file (dimensions [3,features]) (columns SNP_INDEX MAF_CASES WEIGHT(MAF_DIF))")
print( ""+covariates_filename)
print( "	.covar   file (dimensions [samples,covars]) * here just 1 covar = gender")
print( ""+samples_filename)
print( "	.samples file (dimensions [samples,1]) SAMPLE_IDS SAME AS PED_FILE")

min= (time.time()-start)/60
hour= min/60
print( "Execution time: "+str(hour)[0:4]+"hs:"+str(min)[0:4]+"min\n")





