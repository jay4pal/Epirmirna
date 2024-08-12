
##### new gene scoring based on PADOG geneset scroing 
# remove dupicate rows by index ---> rat_db = rat_db[-c(3966),]
# add quant score to data 
##########################################################################
#rat_db$quant_score = apply(df, 1, function(x) {                         #
#  writeLines(sprintf("mirna: %s", x[1]));                               #
#  return((as.numeric(x[2])) * (df2$fc[match(x[1], df2$mir)]))           #
#})                                                                      #
##########################################################################
# Sql inner join - reduces the number of mTIs 
#test1=sqldf("select rat_db.*, (rat_db.Final_score * dg_sele_hm.log2_FC_16day) as quant_score
#       from rat_db inner join dg_sele_hm
#       ON rat_db.miRNA = dg_sele_hm.miRNA_H_M") 
#x[1] = rat_db["miRNA"]
#x[14] = rat_db["Final_score"]
rat_db$quant_score = apply(rat_db,1,function(x){
  return((as.numeric(x[14])) * (dg_sele_hm$log2_FC_16day[match(x[1], dg_sele_hm$rno_mir)]))
})

###############count no of mirnas target s gene in whole data without cut off ############
rat_db_all = sqldf("select S.*, C.n_miRNA
                   from rat_db S
                   inner join (Select target_gene, count(target_gene) as n_miRNA
                   from rat_db
                   group by target_gene) C on S.target_gene = C.target_gene") 

#data with cut off >0.5
rat_db_cut = sqldf("select * from rat_db where 
                   rat_db.Final_score > 0.5")
##########################################################################################
# find frequncy of gene in data set RAT_db
# decide whether u want to add new scoring to just Cut off > 0.5 data or whole data sample 
rat_db1 = sqldf("select S.*, C.n_miRNA
                from rat_db_cut S
                inner join (Select target_gene, count(target_gene) as n_miRNA
                from rat_db_cut
                group by target_gene) C on S.target_gene = C.target_gene") 
# add weight for frequency of genes
# max(rat_db1$n_miRNA) 
# min(rat_db1$n_miRNA) 
#can also try this apply function
for(i in 1:nrow(rat_db1)){
  rat_db1$gene_weight[i] = (1 - sqrt((max(rat_db1$n_miRNA)-rat_db1$n_miRNA[i])/
                                       (max(rat_db1$n_miRNA)-min(rat_db1$n_miRNA))))
}

## scale gene_weight from 1-2 as gene with 1 freq = 0 (product will be 0)
# old rabge 0-1  --------> new range 1-2
for(i in 1:nrow(rat_db1)){
  rat_db1$gene_weight_scale[i] = 1+(((rat_db1$gene_weight[i] - 0)/(1-0)*(2-1)))
  # simply 1 + gene_weight
}


# gene scoring 
# calc mean score if same MTI occurs more than once due to different binding sites
data = rat_db1
rat_db2 = data %>%
  group_by(target_gene,n_miRNA, gene_weight_scale) %>%
  summarise(Genescore = (sum((quant_score*gene_weight_scale), na.rm = TRUE))/
              length(unique(target_gene)))
# merge the genes scores with with gene names
#merge the genescore with gene data 
rat_db1 = merge(rat_db1, rat_db2, by.x = 2, by.y = 1, all.x = TRUE)
#gene_expression like term
rat_db1$gene_exp = rat_db1$Genescore*(-1)
#remove duplicate  columns
rat_db1 = rat_db1[ , !names(rat_db1) %in% c("n_miRNA.y","gene_weight_scale.y")]

rm(rat_db2, rat_db_cut, rat_db)

rat_geneScore = rat_db1[c("target_gene","gene_exp")]
rat_geneScore = rat_geneScore[!duplicated(rat_geneScore),]


#########################################################################################
# similar scoring scheme for mouse data
# dg_sele_hm$Mature_ID_21 => mmu miRNA ids ; chnage colnames if desired by user
mouse_db$quant_score = apply(mouse_db,1,function(x){
  return((as.numeric(x[14])) * (dg_sele_hm$log2_FC_16day[match(x[1], dg_sele_hm$Mature_ID_21)]))
})
mmu_db_all = sqldf("select S.*, C.n_miRNA
                   from mouse_db S
                   inner join (Select target_gene, count(target_gene) as n_miRNA
                   from mouse_db
                   group by target_gene) C on S.target_gene = C.target_gene") 

mouse_db_cut = sqldf("select * from mouse_db where 
                     mouse_db.Final_score > 0.5")

mouse_db1 = sqldf("select S.*, C.n_miRNA
                  from mouse_db_cut S
                  inner join (Select target_gene, count(target_gene) as n_miRNA
                  from mouse_db_cut
                  group by target_gene) C on S.target_gene = C.target_gene") 
# add weight for frequency of genes
# max(rat_db1$n_miRNA) 
# min(rat_db1$n_miRNA) 
#can also try this apply function
for(i in 1:nrow(mouse_db1)){
  mouse_db1$gene_weight[i] = (1 - sqrt((max(mouse_db1$n_miRNA)-mouse_db1$n_miRNA[i])/
                                         (max(mouse_db1$n_miRNA)-min(mouse_db1$n_miRNA))))
}

## scale gene_weight from 1-2 as gene with 1 freq = 0 (product will be 0)
# old rabge 0-1  --------> new range 1-2
for(i in 1:nrow(mouse_db1)){
  mouse_db1$gene_weight_scale[i] = 1+(((mouse_db1$gene_weight[i] - 0)/(1-0)*(2-1)))
  # simply 1 + gene_weight
}


# gene scoring 
# calc mean score if same MTI occurs more than once due to different binding sites
data = mouse_db1
mouse_db2 = data %>%
  group_by(target_gene,n_miRNA, gene_weight_scale) %>%
  summarise(Genescore = (sum((quant_score*gene_weight_scale), na.rm = TRUE))/
              length(unique(target_gene)))
# merge the genes scores with with gene names
#merge the genescore with gene data 
mouse_db1 = merge(mouse_db1, mouse_db2, by.x = 2, by.y = 1, all.x = TRUE)
#gene_expression like term
mouse_db1$gene_exp = mouse_db1$Genescore*(-1)
#remove duplicate  columns
mouse_db1 = mouse_db1[ , !names(mouse_db1) %in% c("n_miRNA.y","gene_weight_scale.y")]
rm(data, mouse_db_cut, mouse_db, mouse_db2)

mmu_geneScore = mouse_db1[c("target_gene","gene_exp")]
mmu_geneScore = mmu_geneScore[!duplicated(mmu_geneScore),]
##############################################################################################
###############################SPIA preparartion ############################################
# script to download Kegg pathways
library(RCurl)
library(XML)
require(stringr)
kegg_rest = "http://www.kegg.jp/kegg-bin/download?entry="
kegg_end = "&format=kgml"
file = "C:/Users/ajaypal/Music/rno_kegg/" 
# Kegg_PathIDs containing list of kegg pathways ids :D:\mirnas94\Kegg_PathIDs.txt
#load this file
ext_mmu = "mmu"
ext_rno = "rno"
ext_hsa = "hsa"
for(i in 1:nrow(Kegg_PathIDs)){
  Sys.setenv(http_proxy="proxy.rcsi-internal.ie:80")
  kegg_url = paste(kegg_rest,ext_hsa,str_sub(Kegg_PathIDs$V1[i],1,5),kegg_end,sep = "")
  file1= paste(file,ext_hsa,Kegg_PathIDs$V1[i],sep = "")
  download.file(url = kegg_url, destfile = file1)}
################################################################################

#use above downloaded pathway kgml file to build global pathway in SPIA
#example : ### for rat - kegg Global pathway generated by for SPIA
library(SPIA)
makeSPIAdata(kgml.path = "folder with kegg kgml files", organism = "rno",
             out.path = "D:/python_comp/pathway_21miR/R_model/")
#output file = D:/python_comp/pathway_21miR/R_model/rnoSPIA.RData

#set the working directory where you need the Csv files to be saved per pathway
#load the RData file created by SPIA for user input pathways from SPIA
#####Gene information from each pathway: # complete optional for gene to gene interactions 
path_gene_info=function(file,org){
  e <- new.env(parent = emptyenv())
  load(file, envir = e)
  obj <- get('path.info', envir =e)
  lapply( names(obj), function(nam) {
    write.csv( obj[[nam]], file=paste(org,nam, ".csv", sep="") )
    #cat(sprintf('%s.csv', nam) )
  })
}
#load the output file #optional for 2nd level of gene-gene (molecular process) interactions
path_gene_info("D:/mirnas94/spia_kegg/mmuSPIA.RData","mmu")
#################################################################################
################### SPIA Pathway analysis ##################################################
library(SPIA)
library(clusterProfiler)
### for rat - kegg pathway data for SPIA
### user defined kegg pathways need to run it once only / when pathways are updated
makeSPIAdata(kgml.path = "D:/mirnas94/spia_kegg/rno_kegg/", organism = "rno",
             out.path = "D:/mirnas94/spia_kegg/")

#########################################################################################
###### convert gene symbols to entrez ids for all predicted genes without cutoff
spia_entrez_all = bitr(unique(rat_db_all$target_gene), fromType = "SYMBOL", 
                       toType = "ENTREZID", annoDb = "org.Rn.eg.db")
#### convert gene symbols to entrez ids  genes with cutoff 0.5
spia_geneset = sqldf("select distinct rat_db1.target_gene, rat_db1.Final_score, rat_db1.gene_exp,
                     rat_db1.Genescore from rat_db1")
spia_entrez = bitr(spia_geneset$target_gene, fromType = "SYMBOL", 
                   toType = "ENTREZID", annoDb = "org.Rn.eg.db")
spia_geneset = merge(spia_entrez,spia_geneset, by.x = 1, by.y = 1, all.x = TRUE)
#temporary data table for entrwz ids and gene scores only
tg1 = sqldf("select distinct spia_geneset.ENTREZID, spia_geneset.gene_exp from spia_geneset")
DE_epi = tg1$gene_exp # genescores in list variable with cut off 0.5
names(DE_epi)<- as.vector(tg1$ENTREZID) # vetorize the scores by entrezids
all_epi = unique(spia_entrez_all$ENTREZID) # entrez ids of all genes in list variable 
res_spia_rat = spia(de = DE_epi, all= all_epi, organism = "rno", 
                    data.dir = "D:/mirnas94/spia_kegg/",nB= 2000, plots = FALSE,
                    beta = NULL, verbose = TRUE, combine = "fisher")

rm(tg1, DE_epi, all_epi,spia_geneset,spia_entrez, spia_entrez_all)

# Script to find out the total number of genes in the pathway
library(stringr)
kegg_rest = "http://rest.kegg.jp/link/genes/"
res_spia_rat$Pathway_genes = "" 
#Sys.setenv(http_proxy="proxy.rcsi-internal.ie:80")

for(i in 1:nrow(res_spia_rat)){
  kegg_url = paste(kegg_rest,"rno",res_spia_rat$ID[i],sep = "")
  Sys.setenv(http_proxy="proxy.rcsi-internal.ie:80")
  res_spia_rat$Pathway_genes[i] = str_count(getURL(kegg_url), "path:") 
}
#convert data from character to numeric
res_spia_rat$Pathway_genes = as.numeric(as.character(res_spia_rat$Pathway_genes))

#Percentage = generatio 
res_spia_rat$Path_percentage = (res_spia_rat$NDE/res_spia_rat$Pathway_genes)*100
#################################SPIA for MOUSE####repeated as above##########################################
### for rat - kegg pathway data for SPIA
### user defined kegg pathways need to run it once only / when pathways are updated
makeSPIAdata(kgml.path = "D:/mirnas94/spia_kegg/mmu_kegg/", organism = "mmu",
             out.path = "D:/mirnas94/spia_kegg/")
#############################################################################################
spia_entrez_all = bitr(unique(mmu_db_all$target_gene), fromType = "SYMBOL", 
                       toType = "ENTREZID", annoDb = "org.Mm.eg.db")
spia_geneset = sqldf("select distinct mouse_db1.target_gene, mouse_db1.Final_score, 
                      mouse_db1.gene_exp,
                     mouse_db1.Genescore from mouse_db1")
spia_entrez = bitr(spia_geneset$target_gene, fromType = "SYMBOL", 
                   toType = "ENTREZID", annoDb = "org.Mm.eg.db")
spia_geneset = merge(spia_entrez,spia_geneset, by.x = 1, by.y = 1, all.x = TRUE)
tg1 = sqldf("select distinct spia_geneset.ENTREZID, spia_geneset.gene_exp from spia_geneset")
DE_epi = tg1$gene_exp
names(DE_epi)<- as.vector(tg1$ENTREZID)
all_epi = unique(spia_entrez_all$ENTREZID)
res_spia_mmu = spia(de = DE_epi, all= all_epi, organism = "mmu", 
                    data.dir = "D:/mirnas94/spia_kegg/",nB= 2000, plots = FALSE,
                    beta = NULL, verbose = TRUE, combine = "fisher")

rm(spia_entrez,spia_entrez_all,spia_geneset,tg1,DE_epi,all_epi,i)

# Script to find out the total number of genes in the pathway
library(stringr)
kegg_rest = "http://rest.kegg.jp/link/genes/"
res_spia_mmu$Pathway_genes = "" 
#Sys.setenv(http_proxy="proxy.rcsi-internal.ie:80")

for(i in 1:nrow(res_spia_mmu)){
  kegg_url = paste(kegg_rest,"mmu",res_spia_mmu$ID[i],sep = "")
  Sys.setenv(http_proxy="proxy.rcsi-internal.ie:80")
  res_spia_mmu$Pathway_genes[i] = str_count(getURL(kegg_url), "path:") 
}
#convert data from character to numeric
res_spia_mmu$Pathway_genes = as.numeric(as.character(res_spia_mmu$Pathway_genes))

#Percentage = generatio 
res_spia_mmu$Path_percentage = (res_spia_mmu$NDE/res_spia_mmu$Pathway_genes)*100
##################################################################################
#use entrezids for David analysis ---
#gene symbols have a limit of 3000 enteries; Entrez ids no limit 
mmu_entrez_ids = bitr(unique(mouse_db1$target_gene), fromType = "SYMBOL", 
                      toType = "ENTREZID", annoDb = "org.Mm.eg.db")
rno_entrez_ids = bitr(unique(rat_db1$target_gene), fromType = "SYMBOL", 
                      toType = "ENTREZID", annoDb = "org.Rn.eg.db")
rm(mmu_entrez_ids, rno_entrez_ids)
#save files to local file and use entrez ids for DAVID functional analysis
#David download Functional annotation file for kegg pathway (file by clicking on bar)
#remove Gene name and Species column and
#replace ,\n [last comma at kegg_pathway column] with \n in EMeditor or elsewhere
##########################################################################################
#DAVID pathway file from DAVID analysis
#rat david file = Davidrno11
#mouse david file = Davidmmu11
# keep check when separating: pathnames don't have , in the name otherwise repalce it with ; etc
#separate column into rows
require(tidyr)
data = Davidrno11
rno_David = data %>%
  mutate(KEGG_PATHWAY = strsplit(as.character(KEGG_PATHWAY),","))%>% 
  unnest(KEGG_PATHWAY)
#mmu_David= mmu_David[!(mmu_David$KEGG_PATHWAY==""),]
#count number of genes in 1 pathway
rno_David = sqldf("select S.*, C.cnt from
                  rno_David S
                  inner join (Select KEGG_PATHWAY, count(KEGG_PATHWAY) as cnt
                  from rno_David
                  group by  KEGG_PATHWAY) C on S. KEGG_PATHWAY = C. KEGG_PATHWAY")

#group genes by pathway
data= rno_David
rno_David = data %>% 
  group_by(KEGG_PATHWAY, cnt) %>%
  summarize(Pathways = paste0(unique(ID), collapse = "+"))
rm(data)
#########################################################################
# Script to find out the total number of genes in the pathway
library(stringr)
kegg_rest = "http://rest.kegg.jp/link/genes/"
rno_David$Pathway_genes = "" 
#Sys.setenv(http_proxy="proxy.rcsi-internal.ie:80")

for(i in 1:nrow(rno_David)){
  kegg_url = paste(kegg_rest,str_sub(rno_David$KEGG_PATHWAY[i],1,8),sep = "")
  Sys.setenv(http_proxy="proxy.rcsi-internal.ie:80")
  rno_David$Pathway_genes[i] = str_count(getURL(kegg_url), "path:") 
}
#convert data from character to numeric
rno_David$Pathway_genes = as.numeric(as.character(rno_David$Pathway_genes))
#############gene ratio
rno_David$Path_percentage = (rno_David$cnt/rno_David$Pathway_genes)*100
rno_David= rno_David[!rno_David$Pathway_genes==0,]
#EASE score / p-value by fisher exact test method can be calculated
#need to know the total number of genes in genome of organism
#fisher.test(rbind(c(297,29960),c(3,40)), alternative="less")
#user genes in pathway = 3 (rno_David$cnt)
#total number of genes in pathway = 40 # rno_David$Pathway_genes
#user input genes not in pathway = 297 (total input - no of genes in pathway)
#genes from total genome not in pathway = 29960 (total genome - total pathway genes)

# to get total genome value: using org.Rn.eg.db
# keep check on this; EASE score will be different if there is ANY UPDATATION of this db
genome = length(unique(keys(org.Rn.eg.db, keytype="ENTREZID")))
rno_David$EASE_score=""

for (i in 1:nrow(rno_David)) {
  first = length(unique(rat_db1$target_gene))-rno_David$cnt[i]
  second = genome - rno_David$Pathway_genes[i]
  rno_David$EASE_score[i] = fisher.test(rbind(c(first,second),c(rno_David$cnt[i],rno_David$Pathway_genes[i])), 
                                     alternative="less")
  } 
rno_David$EASE_score = as.numeric(rno_David$EASE_score)

# benjamini FDR 
rno_David$EASE_score_bh = p.adjust(rno_David$EASE_score, method = "BH", 
                                   n = length(rno_David$EASE_score))
################################################################################################
#Kegg link for pathways with proteins impacted (highlighted)
kegg_path = "http://www.genome.jp/kegg-bin/show_pathway?"
rno_David$Kegg_link = "" 

for(i in 1:nrow(rno_David)){
  kegg_l = paste(kegg_path,str_sub(rno_David$KEGG_PATHWAY[i],1,8),sep = "")
  rno_David$Kegg_link[i] = paste(kegg_l,rno_David$Pathways[i],sep = "+")
}
##########################David mouse data #####################################################
#separate column into rows
require(tidyr)
data = Davidmmu11
mmu_David = data %>%
  mutate(KEGG_PATHWAY = strsplit(as.character(KEGG_PATHWAY),","))%>% 
  unnest(KEGG_PATHWAY)
#mmu_David= mmu_David[!(mmu_David$KEGG_PATHWAY==""),]
#count number of genes in 1 pathway
mmu_David = sqldf("select S.*, C.cnt from
                  mmu_David S
                  inner join (Select KEGG_PATHWAY, count(KEGG_PATHWAY) as cnt
                  from mmu_David
                  group by  KEGG_PATHWAY) C on S. KEGG_PATHWAY = C. KEGG_PATHWAY")

#group genes by pathway
data= mmu_David
mmu_David = data %>% 
  group_by(KEGG_PATHWAY, cnt) %>%
  summarize(Pathways = paste0(unique(ID), collapse = "+"))
rm(data)


#########################################################################
# Script to find out the total number of genes in the pathway
library(stringr)
kegg_rest = "http://rest.kegg.jp/link/genes/"
mmu_David$Pathway_genes = "" 
#Sys.setenv(http_proxy="proxy.rcsi-internal.ie:80")

for(i in 1:nrow(mmu_David)){
  kegg_url = paste(kegg_rest,str_sub(mmu_David$KEGG_PATHWAY[i],1,8),sep = "")
  Sys.setenv(http_proxy="proxy.rcsi-internal.ie:80")
  mmu_David$Pathway_genes[i] = str_count(getURL(kegg_url), "path:") 
}
#convert data from character to numeric
mmu_David$Pathway_genes = as.numeric(as.character(mmu_David$Pathway_genes))
############### gene ratio
mmu_David$Path_percentage = (mmu_David$cnt/mmu_David$Pathway_genes)*100

#################################################################################################
# to get total genome value: using org.Mm.eg.db
# keep check on this; EASE score will be different if there is ANY UPDATATION of this db
genome = length(unique(keys(org.Mm.eg.db, keytype="ENTREZID")))
mmu_David$EASE_score=""

for (i in 1:nrow(mmu_David)) {
  first = length(unique(mouse_db1$target_gene))-mmu_David$cnt[i]
  second = genome - mmu_David$Pathway_genes[i]
  mmu_David$EASE_score[i] = fisher.test(rbind(c(first,second),
                                      c(mmu_David$cnt[i],mmu_David$Pathway_genes[i])), 
                                      alternative="less")
} 
mmu_David$EASE_score = as.numeric(mmu_David$EASE_score)

# benjamini FDR 
mmu_David$EASE_score_bh = p.adjust(mmu_David$EASE_score, method = "BH", 
                                   n = length(mmu_David$EASE_score))
#Kegg link for pthways with highlighted genes
kegg_path = "http://www.genome.jp/kegg-bin/show_pathway?"
mmu_David$Kegg_link = "" 
#Sys.setenv(http_proxy="proxy.rcsi-internal.ie:80")

for(i in 1:nrow(mmu_David)){
  kegg_l = paste(kegg_path,str_sub(mmu_David$KEGG_PATHWAY[i],1,8),sep = "")
  mmu_David$Kegg_link[i] = paste(kegg_l,mmu_David$Pathways[i],sep = "+")
}

rm(first,genome,i,kegg_path,kegg_l, kegg_url, second)
#########################################################################################
#try pathnet on this sample # works with humangenome only
# convert symbol to entrezid
#use david if want to change rat/ symbol to human/ mouse symbols

### pathway analysis using PAthnet
#taking gene data for rat - antilog2 values of expresiion data generated(scores)

testsamp1 = sqldf("select distinct test.targetGene, test.GeneExp from test")
testsamp1$antilogGeneE = 2^(testsamp1$GeneExp)
#convert gene symbols to EntrezIDs
test1 = bitr(testsamp1$targetGene, fromType = "SYMBOL", toType = "ENTREZID",
             annoDb = "org.Rn.eg.db")
#preform for human if reqd = no fine results from rat ids as pathnet use hsa pathways
#put values in david
write.table(testsamp1, file = "C:/Users/ajaypal/Music/test.txt",append = FALSE,
            sep = "\t", eol = "\n", na = "NA", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"))
#david converts 420 out 453 rat genes to human genes
#david ----> Main accessions -------> entrezgene id - download file
#read downloaded file
test1 = read.table("C:/Users/ajaypal/Music/test1.txt",header = T,
                   sep = "\t", stringsAsFactors = FALSE)
#convert to entrezid and emseble id 
#test1 = bitr(test1$ID, fromType = "SYMBOL", 
#             toType = c("ENTREZID","ENSEMBL"), annoDb = "org.Rn.eg.db")
#merge entrez id and genescores

testsamp1 = merge(testsamp1, test1, by.x = "targetGene", by.y = "ID",
                  all.y = TRUE)
rm(test1)
#bring colmn with entrez id to first
col_id = grep("ENTREZ_GENE_ID", names(testsamp1))
testsamp1 = testsamp1[, c(col_id,(1:ncol(testsamp1))[-col_id])]

#PAthnet is just a test run; not a favourable PAthway analysis method 
#pathnet enrichment analysis # works only with human genome data
A = as.matrix(read.table(file= "D:/Downloads/PathNet/adjacency_data.txt", sep = "\t",
                         header = T))
pathway = read.table(file = "D:/Downloads/PathNet/pathway_data.txt", sep = "\t",
                     header = T)
#test1 = read.table(file = "C:/Users/ajaypal/Music/test1.txt", sep = "\t",
#                   header = T, stringsAsFactors = FALSE)
# load Pathnet
#Column_DirectEvidence = 4 column with expression values antilog
load("D:/Downloads/PathNet/PathNet.RData")
PathNet(Enrichment_Analysis = TRUE, DirectEvidence_info = testsamp1,
        Adjacency = A, pathway = pathway, Column_DirectEvidence = 4,
        n_perm = 2000, threshold = 0.05)

#The results are sorted in increasing order based on the p_PathNet_FWER.
#Results include pathway name (Name), number of genes of the pathway present
#in the microarray data (No_of_Genes), number of genes significant from direct
#evidence (Sig_Direct), number of genes significant from combined evidence
#(Sig_Combi), significance of enrichment from hypergeometric test (p_Hyper),
#family wise error rate correction of p_Hyper (p_Hyper_FWER), significance
#of enrichment from PathNet (p_PathNet), family wise error rate correction
#of p_PathNet (p_PathNet_FWER). We used p_PathNet_FWER in our manuscript and
#we recommend users to use p_PathNet_FWER.
Pathnet_result <- read.table("D:/python_comp/pathway_21miR/R_model/PathNet_enrichment_results.txt",
                             header = T,sep = "\t", stringsAsFactors = FALSE)
rm(A, Pathnet_result, pathway,col_id)
##########################################################################################

##########################################################################################
####### Searching transcription factors in MTIs 
# database = http://bioinfo.life.hust.edu.cn/AnimalTFDB/
rat_tf <- read.delim("D:/mirnas94/Rattus_norvegicus_transcription_factors_gene_list.txt", stringsAsFactors=FALSE)
mmu_tf <- read.delim("D:/mirnas94/Mus_musculus_transcription_factors_gene_list.txt", stringsAsFactors=FALSE)

# convert targetgenes to ensemble or Entrez ids
require(clusterProfiler)
rat_ensemble = bitr(unique(rat_db1$target_gene), fromType = "SYMBOL", toType = c("ENTREZID"),
                    annoDb = "org.Rn.eg.db")
rat_ensemble = rat_ensemble[!duplicated(rat_ensemble),]


mmu_ensemble = bitr(unique(mouse_db1$target_gene), fromType = "SYMBOL", toType = c("ENTREZID"),
                    annoDb = "org.Mm.eg.db")
mmu_ensemble = mmu_ensemble[!duplicated(mmu_ensemble),]
#Entrez ids were used # as ensemble ids were having duplicate enteries
rat_tf <- read.delim("D:/mirnas94/Rattus_norvegicus_transcription_factors_gene_list.txt", stringsAsFactors=FALSE)
transcription_factors_rat = Reduce(intersect, list(rat_ensemble$ENTREZID,
                                                   rat_tf$Entrez.ID))
transcription_factors_rat = bitr(transcription_factors_rat, fromType = "ENTREZID", toType = c("SYMBOL","ENSEMBL"),
                                 annoDb = "org.Rn.eg.db")


mmu_tf <- read.delim("D:/mirnas94/Mus_musculus_transcription_factors_gene_list.txt", stringsAsFactors=FALSE)
transcription_factors_mmu = Reduce(intersect, list(mmu_ensemble$ENTREZID,
                                                   mmu_tf$Entrez.ID))
transcription_factors_mmu = bitr(transcription_factors_mmu, fromType = "ENTREZID", toType = c("SYMBOL","ENSEMBL"),
                                 annoDb = "org.Mm.eg.db")
##########################################################################################
# for cytoscape attribute file to type of genes as t-factors
rat_cut_attribute$type[rat_cut_attribute$name %in% transcription_factors_rat$SYMBOL] = "T-factor"
mmu_attrubute$type[mmu_attrubute$name %in% transcription_factors_mmu$SYMBOL] = "T-factor"

###################### Brazil overlap searching #########################################
#########preprocessed file with rat and mmu dat and expression levels in DG region of brain
Brazil_merged <- read.delim("//physiol-dt-32/project-shared/Epilepsy/Brazil_merged.txt",
                            stringsAsFactors=FALSE)

brazil_overlap_rat = sqldf("select distinct Brazil_merged.SYMBOL_RNO, 
                      Brazil_merged.log2FoldChange, rat_db1.gene_exp 
                      from Brazil_merged, rat_db1 
                      where Brazil_merged.SYMBOL_RNO = rat_db1.target_gene COLLATE NOCASE")

b_o_r_correalte = cor(brazil_overlap_rat$log2FoldChange, brazil_overlap_rat$gene_exp,
                      method = "pearson")

####################### mouse data ##################################################
brazil_overlap_mmu = sqldf("select distinct Brazil_merged.SYMBOL_MMU,
                      Brazil_merged.log2FoldChange, mouse_db1.gene_exp from
                      Brazil_merged, mouse_db1 where 
                      Brazil_merged.SYMBOL_MMU = mouse_db1.target_gene COLLATE NOCASE")

b_o_m_correalte = cor(brazil_overlap_rat$log2FoldChange, brazil_overlap_rat$gene_exp,
                      method = "pearson")

########################################################################################
## using entrez ids to get overlap as gene symbols are ambiguous
rat_entrez = sqldf("select rat_db1.target_gene, rat_db1.gene_exp from rat_db1")

library(clusterProfiler) # convert symbols to entrezids
entrez = bitr(unique(rat_entrez$target_gene), fromType = "SYMBOL", toType = c("ENTREZID"),
              annoDb = "org.Rn.eg.db")
rat_entrez = merge(rat_entrez, entrez , by.x = 'target_gene', by.y = 'SYMBOL', all.y = TRUE)

brazil_overlap_rat = sqldf("select distinct Brazil_merged.SYMBOL_RNO, Brazil_merged.ENTREZ_RNO,
      Brazil_merged.log2FoldChange, rat_entrez.gene_exp from Brazil_merged, rat_entrez where
      Brazil_merged.ENTREZ_RNO = rat_entrez.ENTREZID COLLATE NOCASE")

b_o_r_correalte = cor(brazil_overlap_rat$log2FoldChange, brazil_overlap_rat$gene_exp,
                      method = "pearson")

########## mouse
mmu_entrez = sqldf("select mouse_db1.target_gene, mouse_db1.gene_exp from mouse_db1")

entrez = bitr(unique(mmu_entrez$target_gene), fromType = "SYMBOL", toType = c("ENTREZID"),
              annoDb = "org.Mm.eg.db")

mmu_entrez = merge(mmu_entrez, entrez , by.x = 'target_gene', by.y = 'SYMBOL', all.y = TRUE)

brazil_overlap_mmu = sqldf("select distinct Brazil_merged.SYMBOL_MMU, Brazil_merged.ENTREZ_MMU,
  Brazil_merged.log2FoldChange, mmu_entrez.gene_exp from Brazil_merged, mmu_entrez where 
  Brazil_merged.ENTREZ_MMU = mmu_entrez.ENTREZID COLLATE NOCASE")

b_o_m_correalte = cor(brazil_overlap_mmu$log2FoldChange, brazil_overlap_mmu$gene_exp,
                      method = "pearson")

#################################################################################################
# r script to make barplots
# order the bars accordingly
res_Spia_mmu11$PathName =factor(res_Spia_mmu11$PathName,
                                levels = res_Spia_mmu11$PathName)
## generate the plot
ggplot(res_Spia_mmu11, aes(x=PathName, y=Gene_ratio))+
  stat_summary(geom = 'bar', colour="#F7DC6F",fill = "#F7DC6F") + coord_flip() +  
  theme(text=element_text(size=20)) +
  geom_text(aes(label = pG))
