##################################################
## STEP 1: Matrix preparation
##################################################
#    Convert FINAL PED file to Sample X SNP matrix and Case-Control array.
#    ATTENTION: in this step (script ped2scikitlearn.py) the value of missing data is codified!
#    OBS: Missing data on scikit-learn: "Such datasets however are incompatible with scikit-learn 
#         estimators which assume that all values in an array are numerical, and that all have
#         and hold meaning."
#    MISSING GENOTYPES: (i)  can be coded as 3, indicating a different value from the {0,1,2} additive code; 
#                       (ii) can be coded as 0, to be considered as the most frequent genotype following the additive model coding.
#                       Here we will code missing as 3.
#
### WRITE FREQUENCY FILES

# OUTPUTS (frq.cc):
# CHR	Chromosome code
# SNP	Variant identifier
# A1	Allele 1 (usually minor)
# A2	Allele 2 (usually major)
# MAF_A	Allele 1 frequency in cases
# MAF_U	Allele 1 frequency in controls
# NCHROBS_A	Number of case allele observations
# NCHROBS_U	Number of control allele observations

OUTPUTDIR=/home/thais/TCC/genes # path to output files
COHORT=helena #name of the cohort
DATASET=${OUTPUTDIR}/${COHORT}
#EVECPLINK=/home/thais/HelenaBIPMED/Common/FlowFilter/pca/helena.step8-3.eigenvec
EVECPLINK=/home/thais/TCC/genes/helena.step8-3.genes.eigenvec
GENES="/home/thais/TCC/helena.genesEED.txt"

PEDTOSK=ped2scikitlearn.py

mkdir -p ${OUTPUTDIR}/data
mkdir -p ${OUTPUTDIR}/frq

REFALLELES=${DATASET}.step8-3.genes.bim
# --a2-allele [filename] {A2 allele col. number} {variant ID col.}
A2COL=6
VARCOL=2


echo "##### STEP 1: MATRIZ PREPARATION"


MONO=/home/thais/HelenaBIPMED/Common/FlowFilter/manhattan/helena.mono.txt 
POLI=/home/thais/HelenaBIPMED/Common/FlowFilter/manhattan/helena.poli.txt

### Divide in Mono and Poli Datasets
plink --bfile ${DATASET}.step8-3.genes --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${MONO} --make-bed --recodeA --out ${DATASET}.step10.genes.mono

plink --bfile ${DATASET}.step8-3.genes --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${MONO} --make-bed --recodeA --out ${DATASET}.step10.genes.poli

plink --bfile ${DATASET}.step8-3.genes --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${MONO} --make-bed --recodeA --out ${DATASET}.step10.genes

##Datasets
plink --bfile ${DATASET}.step10.genes.poli --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --freq case-control --out ${DATASET}.step10.genes.poli

plink --bfile ${DATASET}.step10.genes.mono --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --freq case-control --out ${DATASET}.step10.genes.mono

plink --bfile ${DATASET}.step10.genes --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --freq case-control --out ${DATASET}.step10.genes

mv *.frq* ${OUTPUTDIR}/frq/

#### CODE MISSING TO 3 TO SET MISSING LATER ON IN R SKAT
#
# Output:
#        raw.data
#        raw.target
#        raw.weights
#        raw.covar
#        raw.samples
#        raw.log
#

##Datasets
LABEL=${COHORT}.step10.genes.poli.frq.cc
python ${PEDTOSK} -m 3 -r ${DATASET}.step10.genes.poli.raw -f ${DATASET}.step10.genes.poli.fam -b ${DATASET}.step10.genes.poli.bim -k ${EVECPLINK} -q ${OUTPUTDIR}/frq/${LABEL}
#ultrarare
LABEL=${COHORT}.step10.genes.mono.frq.cc
python ${PEDTOSK} -m 3 -r ${DATASET}.step10.genes.mono.raw -f ${DATASET}.step10.genes.mono.fam -b ${DATASET}.step10.genes.mono.bim -k ${EVECPLINK} -q ${OUTPUTDIR}/frq/${LABEL}

LABEL=helena.step10.genes.frq.cc
python ${PEDTOSK} -m 3 -r helena.step10.genes.raw -f helena.step10.genes.fam -b helena.step10.genes.bim -k ${EVECPLINK} -q /home/thais/TCC/genes/frq/helena.step10.genes.frq.cc

mv *.target ${OUTPUTDIR}/data
mv *.data ${OUTPUTDIR}/data
mv *.weights ${OUTPUTDIR}/data
mv *.covar ${OUTPUTDIR}/data
mv *.samples ${OUTPUTDIR}/data


##################################################
## STEP 2: Generating Random Forest AUROC and AURP
##################################################

# Change inside the script the path of the dataset
python /home/thais/HelenaBIPMED/RandonForest/scikit-learn-example-RF.py/

echo "Process ended. Finally :) "