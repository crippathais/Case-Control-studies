# Case-Control-studies

Basic flow with the possible softwares for rare variants and for common variants . To deepen and test I also have a script wrote in Bash and also in R Markdown.

```bash
CaseControl_RareVar.sh
```

```bash
CaseControl_CommonVar.sh
```

Auxiliaries scripts (QC pre-processing):
- imiss-vs-het.Rscript
- lmiss-hist.Rscript
- run-IBD-QC.pl
- PCA.Rscript
- MahatPlot.Rscript -- Plot of linear regressionand QQplot
  
Scripts (Analysis)
- CollapseGenes.R -- if one or more genes are in the same SNP run this part.
- run_SKAT_final.R and run_SKAT_finalX.R -- SKAT software ()
- ped2scikitlearn.py 
- getSignificantResults.py -- retrieve significant results from SKAT run.


