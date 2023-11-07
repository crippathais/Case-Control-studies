#################################################
# SKAT Analysis
# Maira R. Rodrigues, Benilton Carvalho, Barbara Henning
# LAB BIOEST. AND BIOCOMP. - IMECC - FCM -UNICAMP
# Updated by Thais Crippa
# Last update: 10.08.2020
#
#################################################
#install.packages("SKAT",type = "source")
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("RCurl")  
# if error, sudo apt-get install libcurl4-gnutls-dev
#BiocManager::install("GenomicRanges")  
##
library(SKAT)


readAndRun <- function(label,ann.tb){
  
  data.file   = paste0(label,".data")
  target.file = paste0(label,".target")
  covars.file = paste0(label,".covar") #PCA
  samples.file= paste0(label,".samples")
  bim.file    = paste0(label,".bim")
  
  # Output files:
  filename = unlist(strsplit(data.file, "/"))
  
  # Reading files
  print("Reading data files...")
  data.tb   <- read.table(file=data.file, na.strings = "3")       # genotype matrix
  bim.tb    <- read.table(file=bim.file)        # snp matrix in map format: chr<TAB>markerID<TAB>genetic_distance<TAB>position
  target.tb <- read.table(file=target.file)     # binary phenotype vector
  covars.tb <- read.table(file=covars.file,
                          colClasses =  c("factor",rep("numeric",10)))# covariates matrix
  covars.tb<-covars.tb[,-1]
  samples.tb<- read.table(file=samples.file)    # samples matrix
  
  # EPIMIRNA DATASET MTLE
  print(paste("-> n_samples =",nrow(data.tb)))  
  print(paste("-> n_snps    =",ncol(data.tb))) 
  
  # ADD marker id as column names and sample id as row names in data.tb
  colnames(bim.tb)  = c("CHR","ID","DISTANCE","POS","A1","A2")
  colnames(data.tb) = bim.tb$ID
  rownames(data.tb) = samples.tb$V1
  rownames(covars.tb) = samples.tb$V1
  
  target.tb = as.numeric(target.tb)
  status = cbind(samples.tb,target.tb)
  colnames(status) = c("sample","status")
  
  #data.sub = as.matrix(data.tb[,!(colnames(data.tb) %in% atgc.tb$V1)])
  #data.tb   = as.matrix(data.sub)
  #ncol(data.tb)
  #bim.tb = bim.tb[!(bim.tb$ID %in% atgc.tb$V1),]
  #nrow(bim.tb)
  data.tb   = as.matrix(data.tb)
  
  # batch.tb$SAMPLE = trimws(batch.tb$SAMPLE)
  # lib = merge(samples.tb, batch.tb, by.x = "V1", by.y = "SAMPLE", sort = FALSE) 
  # covars.tb$lib = as.factor(lib$library)
  
  epimirna = cbind(covars.tb,target.tb)
  colnames(epimirna) = c(sprintf("PC%d", 1:10),"status") 
  epimirna$status = as.factor(epimirna$status)
  
  # rm(list = c('lib'))
  
  
  ########### GENE SETS - SNP SUBSETING
  
  # filter ann.tb by snps in bim file. 
  nrow(ann.tb)
  ann.sub = ann.tb[ann.tb$snp %in% bim.tb$ID,]
  nrow(ann.sub)
  genes = unique(ann.sub$gene)
  length(genes)
  # create list with unique genes
  genes = na.omit(unique(genes))
  length(genes)
  genelist = genes
  
  #######Pegar primeira coluna da Tabela
  
  ############
  # SKAT RUN
  ############
  #
  # Obs on missing values:
  # If there are missing genotypes, SKAT automatically imputes them based on 
  # Hardy-Weinberg equilibrium. You can choose from "bestguess", "fixed" or "random".
  # The "bestguess" imputes missing genotypes as most likely values (0,1,2), 
  # the "fixed"imputes missing genotypes by assigning the mean genotype value (2p, p is the MAF) and
  # the "random" imputes missing genotypes by generating binomial(2,p) random variables. 
  # The default imputation method for the SKAT function is "fixed" and for the SKATBinary function is "bestguess".
  runSKAT <- function(map,data,epimirna,genelist,ann.sub,assoc.out){
    
    #map = map.tb
    #data = nlfe
    #epimirna = epimirna.nlfe
    #genelist = g
    print("Runing SKATO...")
    
    
    start.time <- Sys.time()
    ###################################
    # Separate Autosomes and X chr SNPs
    ###################################
    #xsnps = map[map$CHR == 23,]$ID
    #length(xsnps)
    #data.autosomes = as.matrix(data[,!colnames(data) %in% xsnps])
    #data.x         = as.matrix(data[,colnames(data) %in% xsnps])
    #ncol(data)
    print(paste("Dataset has:",ncol(data),"SNPs"))
    #ncol(data.x)
    
    # null model
    ### Alterar o sexo (tirar da amostra)
    formula = status ~ PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10
    #formula = status ~ sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9 +PC10 + batch
    #formula = status ~ sex +PC1 +PC2
    
    
    ##### MODELO NULO
    epim <-SKAT_Null_Model(formula, data=epimirna, out_type="D")
    #epimx =SKAT_Null_Model_ChrX(formula, SexVar="gender", out_type="D")
    
    # default method is davies
    #methods = c("davies","SKATO")
    # for Binary_Robust
    methods = c("SKATO")
    kernels = c("linear.weighted", "IBS", "IBS.weighted","quadratic", "2wayIX")
    kernels.SKATO = c("linear.weighted")
    
    # Select method by index (m)
    #m = 1 #davies (default)
    #m = 2 #SKATO
    
    # Outputs for Autosomes
    df     = data.frame(gene=character(),SKATO=double(),linear.weighted=double(),IBS =double(),IBS.weighted =double(),quadratic =double(),TwoWayIX=double(),nsnps= integer(),snps=character(),stringsAsFactors = FALSE)
    df.adj = data.frame(gene=character(),SKATO=double(),linear.weighted=double(),IBS =double(),IBS.weighted =double(),quadratic =double(),TwoWayIX =double(),nsnps= integer(),snps=character(),stringsAsFactors = FALSE)
    colnames(df)[7] <- "2wayIX"
    colnames(df.adj)[7] <- "2wayIX"
    # Outputs for ChrX
    #df.x     = data.frame(gene=character(),linear.weighted=double(),IBS =double(),IBS.weighted =double(),quadratic =double(),TwoWayIX =double(),SKATO =double(),nsnps= integer(),snps=character(),stringsAsFactors = FALSE)
    #df.x.adj = data.frame(gene=character(),linear.weighted=double(),IBS =double(),IBS.weighted =double(),quadratic =double(),TwoWayIX =double(),SKATO =double(),nsnps= integer(),snps=character(),stringsAsFactors = FALSE)
    #colnames(df.x)[6] <- "2wayIX"
    #colnames(df.x.adj)[6] <- "2wayIX"
    #assoc.x.out   = paste("epimirna/skat/",tail(filename,1),".chrX.c",c,sep = "")
    
    
    ##### VERIFICAR SE PRECISA DESSA PARTE!!!! SPLIT POR ;
    idx=1
    idx2=1
    total_loop = length(genelist)
    total_snps = 0
    
    for (i in 1:length(genelist)) {
      #i=3
      regex = paste("\\b",genelist[i],"\\b", sep = "")
      match = grepl(regex,ann.sub[, "gene"])
      snps = ann.sub[match,]$snp
      snps = snps[!is.na(snps)]
      snps =  unique(snps)
      length(snps)
      
      
      # RUN SKAT FOR AUTOSOME SNPS
      ############################
      # IMPUTATION: When variates are very rare and missing rates between cases and controls are
      # highly unbalanced, impute.method="fixed" can yield inflated type I error rate. In this case, we
      # recommend to use impute.method="bestguess", which does not suffer the same problem.
      data.sub = as.matrix(data[,colnames(data) %in% snps])
      print(paste("Gene:",genelist[i],"SNPS:",length(snps),"SNPs in dataset:",ncol(data.sub)))
      
      # if has at least 1 (>0) or 2 (>1) snps
      if(ncol(data.sub)>0){
        for(m in 1:length(methods)){
          
          #PARA VARIAR OS M´ETODOS RETIRAR O IF
          if(methods[m]=="SKATO"){    
            for(k in 1:length(kernels.SKATO)){
              print(paste("i=",i,"/",total_loop,"Autosomes","Kernel=",kernels.SKATO[k], "Method=",methods[m]))
              # RUN SKAT ON GENE SUBSET
              #epiSKAT = SKATBinary(data.sub, epim, method=methods[m], kernel = kernels.SKATO[k])
              # type of kernel (default= "linear.weighted"). The possible choices are "linear" and "linear.weighted".
              epiSKAT = SKATBinary_Robust(data.sub, epim, method=methods[m], kernel = kernels.SKATO[k])
              
              df[idx,methods[m]] = epiSKAT$p.value
            }  
          }# end if SKATO 
          
        }# end for methods
        
        df[idx,]$gene = genelist[i]
        df[idx,]$nsnps= ncol(data.sub)
        if(ncol(data.sub)==1){
          df[idx,]$snps = paste(intersect(colnames(data),snps), collapse = " ")
        }else{
          df[idx,]$snps = paste(colnames(data.sub), collapse = " ")
        }
        idx = idx+1
        total_snps=total_snps+ncol(data.sub)
        
        
      }# end if data>1
      else {
        print(paste("Gene:",genelist[i]," has less than 1 SNP for this dataset."))
      }
      
    }# end loop genes
    
    # Results for Autosomes
    ########################
    # Total genelist:
    print(paste("Total autosome genes:",nrow(df)))
    print(paste("Total snps:",total_snps))
    
    # how many genes with significant p-val
    nsigp    = nrow(df[df[,"SKATO"]<0.05,])
    print(paste("pval<0.05:",nsigp))
    # Adjust p-value
    for(k in 1:length(kernels)){
      df.adj[1:nrow(df),kernels[k]] = p.adjust(df[,kernels[k]], method = "bonferroni")
    }  
    df.adj[1:nrow(df),"SKATO"] = p.adjust(df[,"SKATO"], method = "bonferroni")
    df.adj$gene = df$gene
    df.adj$snps = df$snps
    df.adj$nsnps= df$nsnps
    
    # how many genes with significant adjusted p-val
    nsigpadj = nrow(df.adj[df.adj[,"SKATO"]<0.05,])
    print(paste("pval<0.05 bonferroni:",nsigpadj))
    # Order results by first p-value
    df.ordered = df[order(df[,"SKATO"]),]
    df.adj.ordered = df.adj[order(df.adj[,"SKATO"]),]
    
    print(paste("Results for autosome genes writen to:",assoc.out))
    write.table(df.ordered, file = assoc.out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(df.adj.ordered, file = paste(assoc.out,".adj",sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    # Results for Chr X
    ########################
    #print(paste("Total chrX genelist:",nrow(df.x)))
    # p-value correction
    # how many genes with significant p-val
    #nsigp    = nrow(df.x[df.x$linear.weighted<0.05,])
    #print(paste("pval<0.05:",nsigp))
    # Adjust p-value
    #for(k in 1:length(kernels)){
    #  df.x.adj[1:nrow(df.x),kernels[k]] = p.adjust(df.x[,kernels[k]], method = "bonferroni")
    #}  
    #df.x.adj[,"SKATO"] = p.adjust(df.x[,"SKATO"], method = "bonferroni")
    #df.x.adj$gene = df.x$gene
    #df.x.adj$snps = df.x$snps
    #df.x.adj$nsnps= df.x$nsnps
    # how many genes with significant adjusted p-val
    #nsigpadj = nrow(df.x.adj[df.x.adj$linear.weighted<0.05,])
    #print(paste("pval<0.05 bonferroni:",nsigpadj))
    #[1] 
    #[1] 
    
    #df.x.ordered = df.x[order(df.x$linear.weighted),]
    #df.x.adj.ordered = df.x.adj[order(df.x.adj$linear.weighted),]
    
    #print(paste("Results for chrX writen to:",assoc.x.out))
    #write.table(df.x.ordered, file = assoc.x.out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    #write.table(df.x.adj.ordered, file = paste(assoc.x.out,".adj",sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    end.time <- Sys.time()
    time.taken <- difftime(end.time, start.time, units='mins')
    print(paste("---> Finished in :",time.taken,"minutes"))
    
    print(paste("---> Covariates :",formula))
    
    
  } # end skat function
  
  # Call SKAT function 
  filename = unlist(strsplit(data.file, "/"))
  assoc.out   = paste0(outputdir,tail(filename,1),".skat")
  print(assoc.out)
  #assoc.x.out   = paste(outputdir,tail(filename,1),".chrX",sep = "")
  runSKAT(bim.tb,data.tb,epimirna,genelist,ann.sub,assoc.out)
  rm(assoc.out)
  
}

###### RUN General 
outputdir="~/Desktop/SKAT/NovoTest/Geral/"

# This file was created with script map_tiled.annot.py based on tiled.padded100.annot.hg38.bed
ann.file     = "Genelist_Fin_sort_m_wX.txt"
ann.tb  = read.table(file=ann.file,stringsAsFactors = FALSE, header = TRUE)   

# All phenotypes
label = "helenaX.step10.ultra"
readAndRun(label,ann.tb)








