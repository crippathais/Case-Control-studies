############################################################################
#############################################################################
#
# FILTER VCF GENOMIC/EXOME SAMPLES - RARE VARIANTS
#
# AUTHOR: Thais de Oliveira (PostDoc at FMC - UNICAMP)
#
# PAPERS: Anderson, et al., 2010, Nature 
#         Clarke, et al., 2011, PMC
#         Setakis, e tal., 2006, Genomes Research
#         Reed, et al., 2015, Statistics in Medicine
#
# INPUTS: .vcf file name and full path
#         file with sample status (mapping sample_id to status and gender)
#         full path of output directory
#
# REQUIREMENTS:
#          PLINK
#          BCFTOOLS
#          R
#          Perl
#
# FUNCTION: Convert a vcf file to ped/map and apply sample and SNP filtering
#
# Last modified: 12.March.2022
#
#############################################################################
#############################################################################


#############################################################################

# Order of data cleaning steps:
#For each cohort:
# 1.  Convert VCF 2 PED/MAP. VCF quality filters: keep SNPs with GQ > 30 and PASS
# 2.  Differential missingness (detect platform/batch differences between case and control genotype data)
# 3.  Sample filter: Excludes samples with > 10% MISSING genotypes.
#     Site filter: Excludes SNPs with > 5% MISSING genotypes.
# 4.  Sample filter: Heterozygosity (Identification of individuals with elevated missing data rates or outlying heterozygosity rate)
# 5.  Site filter: HWE Filter - (i) filter controls with more stringed value than cases
# 6.  Site filter: LD pruning.
# 7.  Sample filter: Relatedness filter.
# 8.  PCA analysis PRUNING - Recommended to do manually
# 9.  Logistic Regression - For commom variants
# 10. Site Filtering on Controls (rare alleles in Cases)
# 11. Divide in groups (if it's necessary for your own data)
# 12. Matrix preparation
# 13. Creating genelist
# 14. Run SKAT-O

# REQUIREMENTS
# imiss-vs-het.Rscript
# lmiss-hist.Rscript
# run-IBD-QC.pl
# MahatPlot.Rscript
# PCA.Rscript
# ped2scikitlearn.py
# collapseGenes.Rscript
# getSignificantResults.py
# genomic-inflation-lambda-QQplot.r
# run_SKAT.r
# run_SKATX.r

#############################################################################

#### OUTPUT Directory
OUTPUTDIR= # path to output files

##### INPUTS - Edit the Full paths for each variable
COHORT= #name of the cohort
INPUTVCF= #path to input VCF file
CASES= #sample_id sample_id
GENDER=    # sample_id sample_id {0,1,2}
#REMOVE= #samples to remove
DATASET=${OUTPUTDIR}/${COHORT}

### Other necessary scripts
PEDTOSK=ped2scikitlearn.py

#### Create new dir
mkdir -p ${OUTPUTDIR}/pca
mkdir -p ${OUTPUTDIR}/miss
mkdir -p ${OUTPUTDIR}/ibd
mkdir -p ${OUTPUTDIR}/data
mkdir -p ${OUTPUTDIR}/frq
mkdir -p ${OUTPUTDIR}/gene

## Number of samples original VCF files
echo "PRE-STEP" > GeneralLog.txt
echo "Samples" > GeneralLog.txt
bcftools query -l ${INPUTVCF} | wc -l >> GeneralLog.txt
echo "SNPs" > GeneralLog.txt
bcftools view -H ${INPUTVCF} | awk '{print $3}' | wc -l >> GeneralLog.txt

##################################################
## STEP 1: CONVERT VCF TO PED/MAP
##         Filter genotypes with min GQ 30
#          Filter SNPs with PASS filter
#          Keep VCF REF allele 
##################################################

echo "##### STEP 1: VCF FILTERS"

# HOW TO INTERPRET THE GQ
# GQ = SCALED PHRED
# Phred Quality Score	Error	Accuracy (1 - Error)
# 10	1/10 = 10%	90%
# 20	1/100 = 1%	99%
# 30	1/1000 = 0.1%	99.9%
# 40	1/10000 = 0.01%	99.99%
# 50	1/100000 = 0.001%	99.999%
# 60	1/1000000 = 0.0001%	99.9999%

# BIM FILE FORMAT: A text file with no header line, and one line per variant with the following six fields:
#Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
#Variant identifier
#Position in morgans or centimorgans (safe to use dummy value of '0')
#Base-pair coordinate (normally 1-based, but 0 ok; limited to 231-2)
#Allele 1 (corresponding to clear bits in .bed; usually minor)
#Allele 2 (corresponding to set bits in .bed; usually major)

# Note: VCF reference alleles are set to A2 by the autoconverter even when they appear to be minor. However, to maintain backwards compatibility with PLINK 1.07, PLINK 1.9 normally forces major alleles to A2 during its loading sequence. One workaround is permanently keeping the .bim file generated during initial conversion, for use as --a2-allele input whenever the reference sequence needs to be recovered. (If you use this method, note that, when your initial conversion step invokes --make-bed instead of just --out, you also need --keep-allele-order to avoid losing track of reference alleles before the very first write, because --make-bed triggers the regular loading sequence.)

# Note: to set proper IDs for variants which have not been assigned standard IDs.
#       assign them chromosome-and-position-based IDs with --set-missing-var-ids
# Note: To skip variants which failed one or more filters tracked by the FILTER field, use --vcf-filter.

plink --vcf ${INPUTVCF} --keep-allele-order --double-id --biallelic-only strict --allow-extra-chr --allow-no-sex --make-bed --out ${DATASET}

### TO KEEP VCF REF ALLELES #########################
REFALLELES=${DATASET}.bim
# --a2-allele [filename] {A2 allele col. number} {variant ID col.}
A2COL=6
VARCOL=2

plink --vcf ${INPUTVCF} --double-id --biallelic-only strict --allow-extra-chr --allow-no-sex --vcf-min-gq 30 --vcf-filter --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --update-sex ${GENDER} --not-chr MT,Y --set-missing-var-ids @:# --make-pheno ${CASES} '*' --make-bed --out ${DATASET}

echo "STEP 1: SITE FILTER - SNPs with GQ > 30 and PASS" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.log >> GeneralLog.txt


##################################################
## STEP 2: DIFFERENTIAL SNP MISSINGNES
###################################################

echo "##### STEP 2: Differential missingness"

# Given a case/control phenotype, --test-missing tries to detect platform/batch differences between case and control genotype data by performing Fisher's exact test on case/control missing call counts at each variant. (Variants with 0% or 100% missing call frequency are now skipped.) Results are written to plink.missing. Multiple-testing corrected p-values can be obtained from --adjust or permutation tests.
#
#The 'midp' modifier causes Lancaster's mid-p adjustment to be applied to Fisher's exact test. We recommend its use, since the conservative bias (which the mid-p adjustment eliminates) exhibited by the regular Fisher's exact test has no value in this context.
#
#Any variant which comes up as highly significant under this test should be treated with great caution; spurious association results are likely.

#Any variant which comes up as highly significant under this test should be treated with great caution; spurious association results are likely.
plink --bfile ${DATASET} --allow-extra-chr --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --test-missing midp --adjust --allow-no-sex --make-bed --out ${DATASET}.tmp

awk 'BEGIN {FS=" ";OFS="\t"} {if($4 < 0.05) print $2}' ${DATASET}.tmp.missing.adjusted > ${DATASET}.tmp.missing.remove
wc -l ${DATASET}.tmp.missing.remove

# Remove SNPs with missing less than 0.05 in bed format
plink --bfile ${DATASET} --allow-extra-chr  --double-id --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --exclude ${DATASET}.tmp.missing.remove --recode vcf --make-bed --out ${DATASET}.step2

printf "\nSTEP 2: DIFFERENTIAL SNP MISSINGNES\n" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.step2.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.step2.log >> GeneralLog.txt

##################################################
## STEP 3: SAMPLE FILTER: Exclude > 10% MISSING
##         SITE FILTER:   Exclude >  5% MISSING
##################################################

echo "##### STEP 3: SITE AND SAMPLE FILTERS"

plink --bfile ${DATASET}.step2 --geno 0.05 --double-id --allow-extra-chr  --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --out ${DATASET}.step3a

printf "\nSTEP 3: SITE FILTERING\n" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.step3a.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.step3a.log >> GeneralLog.txt

plink --bfile ${DATASET}.step3a --double-id --mind 0.1 --allow-extra-chr  --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --out ${DATASET}.step3

printf "\nSTEP 3: SAMPLE FILTERING\n" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.step3.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.step3.log >> GeneralLog.txt


##################################################
## STEP 4: SAMPLE FILTER: heterozygosity
##################################################

echo "##### STEP 4: SAMPLE FILTER - Heterozygosity"

#Similarly, the distribution of mean heterozygosity (excluding sex chromosomes) across all individuals should be inspected to identify individuals with an excessive or reduced proportion of heterozygote genotypes, which may be indicative of DNA sample contamination or inbreeding, respectively. Mean heterozygosity (which is given by (N − O)/N, where N is the number of non-missing genotypes and O the observed number of homozygous genotypes for a given individual) will differ between populations and SNP genotyping panels.

plink --bfile ${DATASET}.step3 --missing --out ${DATASET}.step3.Miss
plink --bfile ${DATASET}.step3 --het --out ${DATASET}.step3.het

#####Create R graph: raw-GWA-data.imiss-vs-het.pdf
R CMD BATCH /home/thais/HelenaBIPMED/Rare/FlowFilter/scripts/imissX-vs-het.Rscript

#####Create R graph: raw-GWA-data.lmiss-vs-het.pdf --> verify threshold
R CMD BATCH /home/thais/HelenaBIPMED/Rare/FlowFilter/scripts/lmissX-hist.Rscript

# Output imiss-vs-het_z.pdf (graph) and Exclude_SamplesHet.txt (het+2SD & het-2SD)
# Output lmiss.pdf (graph) Exclude_Samples_Miss.txt (log10(0.2)-->3%)
#plink --bfile ${DATASET}.step3 --double-id --biallelic-only strict --allow-extra-chr --allow-no-sex --remove ${OUTPUTDIR}/miss/Exclude_SamplesHet.txt --remove ${OUTPUTDIR}/miss/Exclude_Samples_Miss.txt --make-bed --out ${DATASET}.step3.step4

plink --bfile ${DATASET}.step3 --double-id --biallelic-only strict --allow-extra-chr --allow-no-sex --remove ${OUTPUTDIR}/miss/Exclude_SamplesHet.txt --make-bed --out ${DATASET}.step4

printf "\nSTEP 4: SAMPLE FILTER - HETEROZYGOSITY\n" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.step4.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.step4.log >> GeneralLog.txt

###############################################################
## STEP 5: SITE FILTER - HWE
################################################################

echo "##### STEP 5: SITE FILTER - HWE FILTER"

# --hwe filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold. We recommend setting a low threshold—serious genotyping errors often yield extreme p-values like 1e-50
# 'midp' modifier applies the mid-p adjustment described in Graffelman J, Moreno V (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium. The mid-p adjustment tends to bring the null rejection rate in line with the nominal p-value, and also reduces the filter's tendency to favor retention of variants with missing data.
## as in paper IJMPR Quality control in GWAS 2018 
## HWE p < 1e-10 for cases 
## HWE p < 1e-6  for controls
## Less strict case thr avoids discarding disease-associated SNPs

# filter controls 
plink --bfile ${DATASET}.step4 --filter-controls --hwe 1e-6 midp --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --out ${DATASET}.step5.controls

# filter cases
plink --bfile ${DATASET}.step4 --filter-cases --hwe 1e-10 midp --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --out ${DATASET}.step5.cases


# remove filtered out on controls 
plink --bfile ${DATASET}.step4 --extract <(cut -f2 ${DATASET}.step5.controls.bim) --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --out ${DATASET}.step5.tmp

# removed filtered out on cases for final hwe file.
plink --bfile ${DATASET}.step5.tmp --extract <(cut -f2 ${DATASET}.step5.cases.bim) --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --out ${DATASET}.step5

rm ${DATASET}.step5.controls*
rm ${DATASET}.step5.cases*
rm ${DATASET}.step5.tmp*

printf "\nSTEP 5: SITE FILTER - HWE\n" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.step5.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.step5.log >> GeneralLog.txt


################################################################
## STEP 6: SITE FILTER - LD PRUNING
################################################################

echo "##### STEP 6: SITE FILTER - LD pruning"

# variants that are in approximate linkage equilibrium

plink --bfile ${DATASET}.step5 --indep-pairwise 50 5 0.2 --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --out ${DATASET}.tmp

plink --bfile ${DATASET}.step5 --extract ${DATASET}.tmp.prune.in --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --out ${DATASET}.step6

rm ${DATASET}.tmp.prune.in

printf "\nSTEP 6: SITE FILTER - LD PRUNING\n" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.step6.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.step6.log >> GeneralLog.txt


###############################################################
## STEP 7: SAMPLE FILTER - RELATEDNESS FILTER (has to happen after LD pruning)
###############################################################

echo "##### STEP 7: SAMPLE FILTER - RELATEDNESS PRUNING"

#Relatedness among samples is analyzed by estimating the proportion of identical-by-descent (IBD) alleles between pairs of individuals by $\hat{pi}$ and relatedness matrix (GRM) estimations. In this context, it is expected that third-degree relatives have values of $\hat{pi}$ or GRM = 0.125. Therefore, pairs of samples with $\hat{pi}$ or GRM > 0.125 are removed from the data.
# The expectation is that IBD = 1 for duplicates or monozygotic twins, IBD = 0.5 for first-degree relatives, IBD = 0.25 for second-degree relatives and IBD = 0.125 for third-degree relatives. Due to genotyping error, LD and population structure there is often some variation around these theoretical values and it is typical to remove one individual from each pair with an IBD > 0.1875, which is halfway between the expected IBD for third- and second-degree relatives. For these same reasons an IBD > 0.98 identifies duplicates.

######## estimates IBD from the sample by PLINK
plink --bfile ${DATASET}.step6 --allow-extra-chr --allow-no-sex --genome --out ${DATASET}.step7

cp ${DATASET}.step7.genome test.genome
cp ${DATASET}.step3.Miss.imiss test.imiss

# Run to identify all pairs of individuals with IBD > 0.125
# change directories, if necessary, in perl script
perl /home/thais/HelenaBIPMED/Rare/FlowFilter/scripts/runX-IBD-QC.pl test

plink --bfile ${DATASET}.step6 --double-id --biallelic-only strict --allow-extra-chr --allow-no-sex --remove ${OUTPUTDIR}/ibd/fail-IBD-QC.txt --make-bed --out ${DATASET}.step7-2

mv test* ibd/

printf "\nSTEP 7: SAMPLE FILTER - RELATEDNESS PRUNING\n" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.step7-2.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.step7-2.log >> GeneralLog.txt


##################################################
## STEP 8: PCA analysis PRUNING
## This step is recommended to perform manualy and not in batch.
##################################################

###############################
## PLINK with Visual inspection
###############################

echo "##### STEP 8: PCA analysis PRUNING"

plink --bfile ${DATASET}.step7-2 --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --pca --out ${DATASET}.step8

printf "\nSTEP 8: PCA analysis PRUNING\n" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.step8.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.step8.log >> GeneralLog.txt

mv ${DATASET}.step8.eigenval pca/
mv ${DATASET}.step8.eigenvec pca/

# Run PCA 3D and PCA 2D
R CMD BATCH /home/thais/HelenaBIPMED/Rare/FlowFilter/scripts/PCA.Rscript # USE with graph interface

TOREMOVE=/home/thais/HelenaBIPMED/Rare/FlowFilter/X/pca/ExcludeSamplesPCA.txt

# To remove outliers identified with visual inspection
plink --bfile ${DATASET}.step8 --allow-extra-chr --allow-no-sex --remove ${TOREMOVE} --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed -out ${DATASET}.step8-2

printf "\nSTEP 8: PCA analysis PRUNING - No Outliers PC1\n" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.step8-2.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.step8-2.log >> GeneralLog.txt

# Run pca again
plink --bfile ${DATASET}.step8-2 --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --pca --make-bed -out ${DATASET}.step8-3

mv ${DATASET}.step8-3.eigenval pca/
mv ${DATASET}.step8-3.eigenvec pca/

# Run script PCA 3D and PCA 2D again
# PCA3D.Rscript # USE with graph interface

### Save final EVEC file
EVECPLINK=/home/thais/HelenaBIPMED/Rare/FlowFilter/X/pca/helenaX.step8-3.eigenvec

##################################################
## STEP 9: LOGISTIC REGRESSION - PART 1
## This step is recommended to perform manualy and not in batch.
##################################################

echo "##### STEP 9:  LOGISTIC REGRESSION - PART 1"

######### Logistic Regression, Mahatman Plot and lambda
#Using eigenvec from PCA as covar
plink --bfile ${DATASET}.step8-3 --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --logistic --covar ${EVECPLINK} --ci 0.95 --adjust --out ${DATASET}.step9

printf "\nSTEP 9: LOGISTIC REGRESSION - PART 1\n" >> GeneralLog.txt
grep "Among remaining phenotypes" ${DATASET}.step9.log >> GeneralLog.txt
grep "people pass filters and QC" ${DATASET}.step9.log >> GeneralLog.txt
grep "Genomic inflation est. lambda" ${DATASET}.step9.log >> GeneralLog.txt

# R
R CMD BATCH /home/thais/HelenaBIPMED/Rare/FlowFilter/scripts/MahatPlot.Rscript


##################################################
## STEP 10: SITE FILTER ON CONTROLS
##         MAX-MAF 1% (keeps only rare variants)
##         MAX-MAF 0.05% (keeps only ultra-rare variants)
##################################################

echo "##### STEP 10: MINOR ALLELE FREQUENCY FILTER"

# Rare 
plink --bfile ${DATASET}.step8-3 --filter-controls --max-maf 0.01 --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --out ${DATASET}.step9.controls.rare

plink --bfile ${DATASET}.step8-3 --extract <(cut -f2 ${DATASET}.step9.controls.rare.bim) --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --recodeA --out ${DATASET}.step10.rare

# Ultra-rare
plink --bfile ${DATASET}.step8-3 --filter-controls --max-maf 0.0005 --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --out ${DATASET}.step9.controls.ultra

plink --bfile ${DATASET}.step8-3 --extract <(cut -f2 ${DATASET}.step9.controls.ultra.bim) --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --make-bed --recodeA --out ${DATASET}.step10.ultra


###############################################################
## STEP 11: DIVIDE MONO AND POLIGENIC 
###############################################################

echo "##### STEP 11: DIVIDE IN GROUPS"

MONO=/home/thais/HelenaBIPMED/Common/FlowFilter/manhattan/helena.mono.txt 
POLI=/home/thais/HelenaBIPMED/Common/FlowFilter/manhattan/helena.poli.txt

#POLI
plink --bfile ${DATASET}.step10.rare --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${MONO} --make-bed --recodeA --out ${DATASET}.step11.poli.rare

plink --bfile ${DATASET}.step10.ultra --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${MONO} --make-bed --recodeA --out ${DATASET}.step11.poli.ultra

#MONO
plink --bfile ${DATASET}.step10.rare --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${POLI} --make-bed --recodeA --out ${DATASET}.step11.mono.rare

plink --bfile ${DATASET}.step10.ultra --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${POLI} --make-bed --recodeA --out ${DATASET}.step11.mono.ultra

#### Recoding VCF files

#POLI
#plink --bfile ${DATASET}.step10.rare --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${MONO} --make-bed --recode vcf --out ${DATASET}.step11.poli.rare

#plink --bfile ${DATASET}.step10.ultra --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${MONO} --make-bed --recode vcf --out ${DATASET}.step11.poli.ultra

#MONO
#plink --bfile ${DATASET}.step10.rare --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${POLI} --make-bed --recode vcf --out ${DATASET}.step11.mono.rare

#plink --bfile ${DATASET}.step10.ultra --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --remove ${POLI} --make-bed --recode vcf --out ${DATASET}.step11.mono.ultra


##################################################
## STEP 12: Matrix preparation
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

echo "##### STEP 12: MATRIZ PREPARATION"

plink --bfile ${DATASET}.step10.rare --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --freq case-control --out ${DATASET}.step12.rare

plink --bfile ${DATASET}.step10.ultra --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --freq case-control --out ${DATASET}.step12.ultra

## subphenotypes
#RARE
#mono
plink --bfile ${DATASET}.step11.mono.rare --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --freq  case-control --out ${DATASET}.step12.mono.rare
#poli
plink --bfile ${DATASET}.step11.poli.rare --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --freq case-control --out ${DATASET}.step12.poli.rare

#ULTRARARE
#Poli
plink --bfile ${DATASET}.step11.poli.ultra --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --freq case-control --out ${DATASET}.step12.poli.ultra
#Mono
plink --bfile ${DATASET}.step11.mono.ultra --double-id --allow-extra-chr --allow-no-sex --a2-allele ${REFALLELES} ${A2COL} ${VARCOL} --freq case-control  --out ${DATASET}.step12.mono.ultra

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

#Normal
#Rare
LABEL=${COHORT}.step12.rare.frq.cc
python ${PEDTOSK} -m 3 -r ${DATASET}.step10.rare.raw -f ${DATASET}.step10.rare.fam -b ${DATASET}.step10.rare.bim -k ${EVECPLINK} -q ${OUTPUTDIR}/frq/${LABEL}
#ultrarare
LABEL=${COHORT}.step12.ultra.frq.cc
python ${PEDTOSK} -m 3 -r ${DATASET}.step10.ultra.raw -f ${DATASET}.step10.ultra.fam -b ${DATASET}.step10.ultra.bim -k ${EVECPLINK} -q ${OUTPUTDIR}/frq/${LABEL}

#RARE - MONO
LABEL=${COHORT}.step12.mono.rare.frq.cc
python ${PEDTOSK} -m 3 -r ${DATASET}.step11.mono.rare.raw -f ${DATASET}.step11.mono.rare.fam -b ${DATASET}.step11.mono.rare.bim -k ${EVECPLINK} -q ${OUTPUTDIR}/frq/${LABEL}

#RARE - POLI
LABEL=${COHORT}.step12.poli.rare.frq.cc
python ${PEDTOSK} -m 3 -r ${DATASET}.step11.poli.rare.raw -f ${DATASET}.step11.poli.rare.fam -b ${DATASET}.step11.poli.rare.bim -k ${EVECPLINK} -q ${OUTPUTDIR}/frq/${LABEL}

#ULTRARARE - MONO
LABEL=${COHORT}.step12.mono.ultra.frq.cc
python ${PEDTOSK} -m 3 -r ${DATASET}.step11.mono.ultra.raw -f ${DATASET}.step11.mono.ultra.fam -b ${DATASET}.step11.mono.ultra.bim -k ${EVECPLINK} -q ${OUTPUTDIR}/frq/${LABEL}

#ULTRARARE - POLI
LABEL=${COHORT}.step12.poli.ultra.frq.cc
python ${PEDTOSK} -m 3 -r ${DATASET}.step11.poli.ultra.raw -f ${DATASET}.step11.poli.ultra.fam -b ${DATASET}.step11.poli.ultra.bim -k ${EVECPLINK} -q ${OUTPUTDIR}/frq/${LABEL}

mv *.target ${OUTPUTDIR}/data
mv *.data ${OUTPUTDIR}/data
mv *.weights ${OUTPUTDIR}/data
mv *.covar ${OUTPUTDIR}/data
mv *.samples ${OUTPUTDIR}/data

echo "Process ended. Finally :) "

##################################################
## STEP 13: Creating a Genelist file
##################################################

echo "##### STEP 13: GENELIST CREATION"

##Create genelist files from BIM files

awk '{print $2}' ${DATASET}.step11.mono.ultra.bim > ${DATASET}.step11.mono.ultra.genelist.txt
awk '{print $2}' ${DATASET}.step11.poli.ultra.bim > ${DATASET}.step11.poli.ultra.genelist.txt
awk '{print $2}' ${DATASET}.step11.mono.rare.bim > ${DATASET}.step11.mono.rare.genelist.txt
awk '{print $2}' ${DATASET}.step11.mono.rare.bim > ${DATASET}.step11.poli.rare.genelist.txt

HEAD=/home/thais/HelenaBIPMED/Gotenks/merge_helenaBIPMED_annot.vcf.txt
OUT=${OUTPUTDIR}/gene/genelist

### Creating Gene List (CHR, SNP, POS, GENE) for SKAT-O
## Here use your own gene set
for m in AARS1 ACER3 ACY1 ADAM22 ADGRG1 ADGRL2 ADGRV1 ADNP ADORA2A ADRA2B ADSL AHDC1 AIMP1 AIMP2 ALDH7A1 ALG13 AMPD2 AMT ANKLE2 ANKRD11 AP2M1 AP3B2 ARFGEF2 ARHGEF15 ARHGEF9 ARID1A ARID1B ARV1 ARX ASAH1 ASPM ASXL1 ASXL3 ATN1 ATP13A2 ATP1A2 ATP1A3 ATP6AP2 ATP6V1A ATP7A ATP8A2 ATRX AUTS2 BCL11A BRAF BRAT1 BSCL2 CACNA1A CACNA1D CACNA1E CACNA1H CACNA2D2 CACNB4 CAD CASK CASR CBL CCDC88A CCDC88C CDK13 CDK19 CDK5 CDKL5 CENPE CENPJ CERS1 CHAMP1 CHD2 CHD4 CHD8 CHRNA2 CHRNA4 CHRNA7 CHRNB2 CILK1 CLCN2 CLCN4 CLDN5 CLN3 CLN5 CLN6 CLN8 CLP1 CLTC CNKSR2 CNNM2 CNPY3 CNTN2 CNTNAP2 COL4A2 CPA6 CPLX1 CPT2 CREBBP CSNK1E CSNK1G1 CSTB CTSD CYFIP2 DALRD3 DCX DEAF1 DENND5A DEPDC5 DHDDS DIAPH1 DIP2A DLAT DNAJC5 DNM1 DNM1L DOCK7 EEF1A2 EFHC1 EHMT1 EIF2S3 ELP4 EMX2 ENG EP300 EPM2A EPRS1 ERBB4 ERMARD EXOSC3 EXT2 FAN1 FARS2 FASN FGF12 FIG4 FLNA FMR1 FOLR1 FOXG1 FRRS1L GABBR2 GABRA1 GABRA2 GABRA3 GABRA5 GABRB1 GABRB2 GABRB3 GABRD GABRG2 GAD1 GAL GAMT GATAD2B GATM GBA GCSH GFAP GLDC GLS GNAO1 GNB5 GOSR2 GPAA1 GRIA3 GRIA4 GRIN1 GRIN2A GRIN2B GRIN2D GRN GUF1 HACE1 HADH HADHB HCN1 HDAC4 HECW2 HERC1 HEXA HEXB HNRNPU HOXD@ IER3IP1 IQSEC2 IRF2BPL ITPA KANSL1 KATNB1 KCNA1 KCNA2 KCNB1 KCNC1 KCNC2 KCND3 KCNH1 KCNH5 KCNJ10 KCNMA1 KCNQ2 KCNQ3 KCNQ5 KCNT1 KCNT2 KCTD17 KCTD3 KCTD7 KIF11 KIF2A KIF5C KIFBP KLF13 LAMB1 LAMC3 LGI1 LIAS LMNB2 MAGI2 MAPK10 MBD5 MDH2 ME2 MECP2 MED17 MEF2C MFSD2A MFSD8 MOCS1 MOCS2 MPDZ MTHFR MTMR1 MTOR NACC1 NDE1 NECAP1 NEDD4L NEUROD2 NEXMIF NHLRC1 NPC1 NPC2 NPRL2 NPRL3 NR4A2 NRG2 NRXN1 NRXN2 NSDHL NSF NTRK2 NUS1 OCLN OPHN1 OR10H2 OTUD6B OTUD7A PACS2 PAFAH1B1 PAX6 PCDH12 PCDH19 PCDHG@ PCLO PDHA1 PDHX PDP1 PHGDH PIGA PIGC PIGN PIGP PIGQ PIGT PIK3R2 PLAA PLCB1 PLEKHG2 PLPBP PNKP PNPO POLG POLG2 PPP1R15B PPP3CA PPT1 PRDM8 PRICKLE1 PRICKLE2 PRRT2 PTCH1 PTPN23 QARS1 RAB11A RARS2 RB1 RBFOX1 RELN RHOBTB2 RNASEH2A RNASEH2B RNASEH2C ROGDI RORB RPH3A RTN4IP1 RTTN RYR3 SAMHD1 SARS1 SASS6 SCARB2 SCN10A SCN1A SCN1B SCN2A SCN3A SCN8A SCN9A SEPSECS SERPINI1 SESN3 SGCE SHH SHOX SIK1 SIX3 SLC12A5 SLC12A6 SLC13A5 SLC16A2 SLC19A3 SLC1A2 SLC20A2 SLC25A12 SLC25A22 SLC2A1 SLC35A2 SLC35A3 SLC45A1 SLC6A1 SLC6A8 SLC6A9 SLC9A6 SMARCA2 SMARCA4 SMARCB1 SMARCE1 SMC1A SMS SNAP25 SNIP1 SNX27 SPAST SPATA5 SPG7 SPTAN1 SRGAP2 SRPX2 ST3GAL3 ST3GAL5 STAMBP STRADA STX1B STXBP1 SUOX SYN1 SYNGAP1 SYNJ1 SZT2 TBC1D24 TBCD TBCE TBL1XR1 TCF4 TMTC3 TNK2 TOR1A TPP1 TREX1 TRIM8 TRIO TRMT10A TRMT9B TRPM1 TSC1 TSC2 TSEN15 TSEN2 TSEN54 TUBA1A TUBA8 TUBB2A TUBB2B TUBB3 TUBG1 TWNK UBA5 UBE2A UBE3A UFC1 UFM1 UGDH UGP2 VARS1 VPS53 VRK2 WASF1 WDR45 WDR45B WDR62 WDR73 WWOX XPR1 YWHAG ZEB2 ZMYND8 ZNF182; 
do
    grep -w "GENEINFO=${m}" ${HEAD} | awk '{print $1 "\t" $3 "\t" $2}' | awk '{print $0, "'${m}'"}' >> ${OUT}.txt
done


### Changing names of SNPs
grep -v "\." ${OUT}.txt > ${OUT}_WoPoint.txt
grep "\." ${OUT}.txt | awk '{$3=$1":"$2} 1' > ${OUT}_WithPoint.txt
cat ${OUT}_WoPoint.txt ${OUT}_WithPoint.txt > ${OUT}_Fin.txt
LIST=sort -k1,1 -k3,3 ${OUT}_Fin.txt > ${OUT}_Fin_sort.txt

#${DATASET}.step11.mono.ultra.genelist.txt
#${DATASET}.step11.poli.ultra.genelist.txt
#${DATASET}.step11.mono.rare.genelist.txt
#${DATASET}.step11.poli.rare.genelist.txt

# Separate genes in autossomes and genes in chrX
for inp in *genelist.txt
do
 for typ in poli mono
 do
    for var in ultra rare
    do
        grep -wf ${inp} ${LIST} > ${DATASET}.step13.${typ}.${var}.genelist.txt
        grep "" ${DATASET}.step13.${typ}.${var}.genelist.txt > ${DATASET}.step13.${typ}.${var}.genelist_chrX.txt
        grep -v "" ${DATASET}.step13.${typ}.${var}.genelist.txt > ${DATASET}.step13.${typ}.${var}.genelist_auto.txt
    done
 done
done

printf "\nSTEP 13: GENELIST \n" >> GeneralLog.txt
printf "\nGenelist - mono.ultra - step 11 \n" >> GeneralLog.txt
wc -l ${DATASET}.step11.mono.ultra.genelist.txt >> GeneralLog.txt
printf "\nGenelist - mono.ultra - step 13 \n" >> GeneralLog.txt
wc -l ${DATASET}.step13.mono.ultra.genelist.txt >> GeneralLog.txt
printf "\nGenelist - poli.ultra - step 11 \n" >> GeneralLog.txt
wc -l ${DATASET}.step11.poli.ultra.genelist.txt >> GeneralLog.txt 
printf "\nGenelist - poli.ultra - step 13 \n" >> GeneralLog.txt
wc -l ${DATASET}.step13.poli.ultra.genelist.txt >> GeneralLog.txt
printf "\nGenelist - mono.rare - step 11 \n" >> GeneralLog.txt
wc -l ${DATASET}.step11.mono.rare.genelist.txt >> GeneralLog.txt 
printf "\nGenelist - mono.rare - step 13 \n" >> GeneralLog.txt
wc -l ${DATASET}.step13.mono.rare.genelist.txt >> GeneralLog.txt
printf "\nGenelist - poli.rare - step 11 \n" >> GeneralLog.txt
wc -l ${DATASET}.step11.poli.rare.genelist.txt >> GeneralLog.txt
printf "\nGenelist - poli.rare - step 13 \n" >> GeneralLog.txt
wc -l ${DATASET}.step13.poli.rare.genelist.txt >> GeneralLog.txt

###if necessary: if one or more genes are in the same SNP run this part.
# R CMD BATCH /home/thais/HelenaBIPMED/Rare/FlowFilter/scripts/collapseGenes.Rscript

##################################################
## STEP 14: SKAT-O
##################################################

echo "##### STEP 14: SKAT-O ANALYSES"

# RUN SKAT-O ANALYSES (scripts/) for autossomes:
R CMD run_SKAT.r

# RUN SKAT-O ANALYSES (scripts/) for chrX:
R CMD run_SKAT_X.r

# List significant results:
DIR= # output directory
python /home/thais/HelenaBIPMED/Rare/FlowFilter/scripts/getSignificantResults.py -t 0.00001 -c Rare -d ${DIR}

# For significant results run:
genomic-inflation-lambda-QQplot.r

# CLEAN THE HOUSE
rm *.tmp*
rm *.nosex*
rm *.log*
rm *.gif*
rm *.hh*

