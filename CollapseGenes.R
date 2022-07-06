#####
# 
######
pacotes <- c("readr", "dplyr", "tidyr", "stringr")

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}

tb1<-read_tsv("/home/thais/HelenaBIPMED/Rare/FlowFilter/X/helenaX.step13.mono.ultra.genelist.txt")
tb2<-read_tsv("/home/thais/HelenaBIPMED/Rare/FlowFilter/X/helenaX.step13.mono.rare.genelist.txt")
tb3<-read_tsv("/home/thais/HelenaBIPMED/Rare/FlowFilter/X/helenaX.step13.poli.ultra.genelist.txt")
tb4<-read_tsv("/home/thais/HelenaBIPMED/Rare/FlowFilter/X/helenaX.step13.poli.rare.genelist.txt")

#tb %>% 
  #separate_rows(gene,sep=";") %>% 
  #write_tsv("genesbyRows.txt")

#tb<-read_tsv("genesbyRows.txt")


tb1 %>% 
  group_by(chr, snp, pos) %>%
  summarise(gene = paste(gene, collapse = ";")) %>%
  write_tsv("helenaX.step13.mono.ultra.collapsed_genes.txt")

tb2 %>% 
  group_by(chr, snp, pos) %>%
  summarise(gene = paste(gene, collapse = ";")) %>%
  write_tsv("helenaX.step13.mono.rare.collapsed_genes.txt")

tb3 %>% 
  group_by(chr, snp, pos) %>%
  summarise(gene = paste(gene, collapse = ";")) %>%
  write_tsv("helenaX.step13.poli.ultra.collapsed_genes.txt")

tb4 %>% 
  group_by(chr, snp, pos) %>%
  summarise(gene = paste(gene, collapse = ";")) %>%
  write_tsv("helenaX.step13.poli.rare.collapsed_genes.txt")

