############################################################
# ajaypal@rcsi.ie
############################################################


library(RCurl)
library(XML)
library(dplyr)
require(tidyr)
library(sqldf)
require(gdata)
perlPAth = "C:/Strawberry/perl/bin/perl5.24.0.exe"


##to  write any table to file 

# get list of miRNAs, rat/mouse
mir_list <- read.csv("D:/python_comp/pathway_21miR/mir_list.txt", sep="",
                         stringsAsFactors=FALSE)
#read mirbase21 to find mouse Mirna with same rat mirna fasta seq
mirbase21 <- read.delim("D:/python_comp/mirTARBASE/mirbase21.txt", stringsAsFactors=FALSE)

# extracting for hsa+mmu+rno
mirbase21 = sqldf("select * from mirbase21 where 
                  mirbase21.Mature_ID_21 like '%hsa%' or
                  mirbase21.Mature_ID_21 like '%mmu%' or 
                  mirbase21.Mature_ID_21 like '%rno%'")

rat_mirna = sqldf("select mirbase21.Mature_ID_21, mirbase21.Mature_Seq_21
                  from mir_list, mirbase21
                  where mirbase21.Mature_ID_21 like '%rno%'
                  and mir_list.miRNA_ID = mirbase21.Mature_ID_21")

#sele mouse mirnas with same fasta seq
mmu_mirna = sqldf("select mirbase21.Mature_ID_21, rat_mirna.Mature_Seq_21
                  from rat_mirna, mirbase21
                  where (rat_mirna.Mature_Seq_21 = mirbase21.Mature_Seq_21) 
                  and mirbase21.Mature_ID_21 like '%mmu%'")
miRNAs_r_m = rbind(rat_mirna, mmu_mirna)
rm(rat_mirna, mmu_mirna, mir_list)
colnames(miRNAs_r_m) = c("miRNA_ID","miR_seq")
#check for matureIDs of miRNAs w.r.t-the name version used in database TargetScan/miRDB/Diana
#TargetScan 
# Family + Seed + merged using slq db names from stefan's MTIs db
# create table TargetScan_h_m_r as 
#select S."MiRBase ID", S."Seed+m8" ,ts_summary_mmu.* from ts_summary_mmu, 
#(select distinct ts_miRFamily."MiRBase ID", ts_miRFamily."Seed+m8" from ts_miRFamily 
# where ts_miRFamily."MiRBase ID" like '%mmu%' or ts_miRFamily."MiRBase ID" like '%rno%' or 
#       ts_miRFamily."MiRBase ID" like '%hsa%') S
#where S."Seed+m8" = ts_summary_mmu."miRNA family" AND 
#((ts_summary_mmu."Representative miRNA" like '%hsa%' and S."MiRBase ID" like '%hsa%') or 
#(ts_summary_mmu."Representative miRNA" like '%mmu%' and S."MiRBase ID" like '%mmu%') or 
# (ts_summary_mmu."Representative miRNA" like '%rno%' and S."MiRBase ID" like '%rno%'))
############################################################################################
### Same in r 
#targetscan complete file preprocessed no need for following steps
targetScan = read.delim("D:/mirnas94/targetscan_mmu_jan16/targetscan_complete.txt",
                        stringsAsFactors=FALSE)
#######################################################################

ts_miRFamily= read.delim("D:/mirnas94/targetscan_mmu_jan16/miR_Family_Info.txt",
                         stringsAsFactors=FALSE)
# summary count targetScan human 7.1 (edition June 2016)
ts_summary= read.delim("D:/mirnas94/targetscan_mmu_jan16/Summary_Counts.all_predictions.txt",
                       stringsAsFactors=FALSE)
#chnge colnames to remove spaces
colnames(ts_miRFamily) = c("miR_family","Seed_m8","Species_ID","MiRBase_ID",        
                           "Mature_sequence","FamilyConservation","MiRBaseAccession")
colnames(ts_summary) = c("Transcript_ID","Gene_Symbol",                       
                         "miRNA_family","Species_ID",                          
                         "Total_conserved_sites","conserved_8mer" ,     
                         "conserved_7mer_m8","conserved_7mer_1a" ,  
                         "Total_nonconserved_sites","nonconserved_8mer",   
                         "nonconserved_7mer_m8","nonconserved_7mer_1a",
                         "Number_6mer","Representative_miRNA" ,               
                         "Total_context_score","Cum_weighted_score" ,
                         "Aggregate_PCT" )
## remove rows/colmns not reqd --> decrease computation load
ts_summary = ts_summary[ , !names(ts_summary) %in% 
                           c("Total_conserved_sites","conserved_8mer" ,     
                             "conserved_7mer_m8","conserved_7mer_1a" ,  
                             "Total_nonconserved_sites","nonconserved_8mer",   
                             "nonconserved_7mer_m8","nonconserved_7mer_1a",
                             "Number_6mer","Aggregate_PCT")]
ts_summary = sqldf("select ts_summary.* from ts_summary where
                   (ts_summary.Representative_miRNA like '%hsa%' and ts_summary.Species_ID = 9606) or
                   (ts_summary.Representative_miRNA like '%mmu%' and ts_summary.Species_ID = 10090) or
                   (ts_summary.Representative_miRNA like '%rno%' and ts_summary.Species_ID = 10116)")

targetScan = sqldf("select S.MiRBase_ID, S.Seed_m8 ,ts_summary.* from ts_summary, 
                   (select distinct ts_miRFamily.MiRBase_ID, ts_miRFamily.Seed_m8 
                   from ts_miRFamily 
                   where ts_miRFamily.MiRBase_ID like '%mmu%' 
                   or ts_miRFamily.MiRBase_ID like '%rno%' 
                   or ts_miRFamily.MiRBase_ID like '%hsa%') S
                   where S.Seed_m8 = ts_summary.miRNA_family AND 
                   ((ts_summary.Representative_miRNA like '%hsa%' and S.MiRBase_ID like '%hsa%') 
                   or (ts_summary.Representative_miRNA like '%mmu%' and S.MiRBase_ID like '%mmu%')
                   or (ts_summary.Representative_miRNA like '%rno%' and S.MiRBase_ID like '%rno%'))"
)

rm(ts_miRFamily, ts_summary)
## remove colmns not reqd --> decrease computation load
targetScan = targetScan[ , !names(targetScan) %in% c("Seed_m8","Transcript_ID",
                                                     "miRNA_family","Representative_miRNA")]
# merge miRNA_Gene to create keyword
targetScan$miR_target = paste(targetScan$MiRBase_ID,
                              targetScan$Gene_Symbol)
targetScan[, "miR_target"] = targetScan$miR_target

#### speciesID and mirna prefix check same
targetScan = sqldf("select TargetScan_mod.* from TargetScan_mod where
                   (TargetScan_mod.MiRBase_ID like '%hsa%' and TargetScan_mod.Species_ID = 9606) or
                   (TargetScan_mod.MiRBase_ID like '%mmu%' and TargetScan_mod.Species_ID = 10116) or
                   (TargetScan_mod.MiRBase_ID like '%rno%' and TargetScan_mod.Species_ID = 10090)")

############################################################################################
#### Mirdb target predictions
#Python script is also available to do the below conversion
####script to convert miRDB db - refseq ids to gense symbols and names
#mirdb complete file preprocessed no need for following steps
miRDB = read.delim("D:/mirnas94/targetscan_mmu_jan16/miRDB_complete.txt",
                        stringsAsFactors=FALSE)
# load database downloaded from mirDB file
miRDB_prediction_result <- read.delim("C:/Users/ajaypal/Downloads/miRDB_prediction_result.txt",
                                      header=FALSE, stringsAsFactors=FALSE)

# separate database according to species 
#rat
mirdb_rno = sqldf("select * from miRDB_prediction_result 
                  where miRDB_prediction_result.V1 like '%rno%'")
#mouse
mirdb_mmu = sqldf("select * from miRDB_prediction_result 
                  where miRDB_prediction_result.V1 like '%mmu%'")
#human
mirdb_hsa = sqldf("select * from miRDB_prediction_result 
                  where miRDB_prediction_result.V1 like '%hsa%'")

############ convert ref seq ids to genenames
require(clusterProfiler)
# rat, mouse, human 
mirdbRAT = bitr(mirdb_rno$V2, fromType = "REFSEQ", 
                toType = c("ENSEMBL", "GENENAME", "SYMBOL", "ENTREZID"), 
                annoDb = "org.Rn.eg.db")
mirdbMOUSE = bitr(mirdb_mmu$V2, fromType = "REFSEQ", 
                  toType = c("ENSEMBL", "GENENAME", "SYMBOL", "ENTREZID"), 
                  annoDb = "org.Mm.eg.db")
mirdbHUMAN = bitr(mirdb_hsa$V2, fromType = "REFSEQ", 
                  toType = c("ENSEMBL", "GENENAME", "SYMBOL", "ENTREZID"), 
                  annoDb = "org.Hs.eg.db")


# merge Gene info to scoring database
mirdb_rno = merge(mirdb_rno, mirdbRAT, by.x = 2, by.y = 1, all.x = TRUE)
mirdb_mmu = merge(mirdb_mmu, mirdbMOUSE, by.x = 2, by.y = 1, all.x = TRUE)
mirdb_hsa = merge(mirdb_hsa, mirdbHUMAN, by.x = 2, by.y = 1, all.x = TRUE)

# merge these three columns into one or treat one by one
#assuming they have similar co7ulmn names # bind_rows
library(dplyr)
miRDB_predicts = bind_rows(mirdb_hsa, mirdb_mmu, mirdb_rno)

rm(mirdbHUMAN, mirdbMOUSE, mirdbRAT, miRDB_prediction_result, 
   mirdb_rno, mirdb_hsa, mirdb_mmu)
# 1 MTI having different scorings # due to diff refseq ids 
#keep checking for duplicate gene enteries

# can change the column names if desire to more precise names
colnames(miRDB_predicts) = c("RefseqID", "miRNA", "miRDB_score", "ENSEMBL",
                             "GENENAME","GeneSymbol","ENTREZID")
#remove refseq ids - column name

miRDB_predicts = miRDB_predicts[ , !names(miRDB_predicts) %in% c("RefseqID","ENSEMBL")]

#make miR_target keyword
miRDB_predicts$miR_target = paste(miRDB_predicts$miRNA,miRDB_predicts$GeneSymbol)
miRDB_predicts[, "miR_target"] = miRDB_predicts$miR_target
miRDB_predicts$miRDB_score = as.numeric(miRDB_predicts$miRDB_score)
#remove enteries where GeneName or GeneSymbol is not available 
miRDB = miRDB_predicts[complete.cases(miRDB_predicts), ]

# calc mean score if same MTI occurs more than once due to different binding sites
data = miRDB
miRDB_predicts = data %>%
  group_by(miR_target,miRNA,GeneSymbol) %>%
  summarise(miRDBscores = mean(miRDB_score))

rm(data,miRDB)
miRDB_predicts = miRDB_predicts[!names(miRDB_predicts) %in% c("miRDB_score")]
miRDB_predicts = miRDB_predicts[!duplicated(miRDB_predicts),]
########################################################################
# query for miRNAs

temp_TS = sqldf("select * from targetscan
                where exists
                (select miRNAs_r_m.miRNA_ID from miRNAs_r_m where 
                targetscan.MiRBase_ID = miRNAs_r_m.miRNA_ID)") 
rm(targetscan)

temp_mirdb = sqldf("select miRDB.* from miRDB
                   where exists
                   (select miRNAs_r_m.miRNA_ID from miRNAs_r_m where 
                   miRDB.miRNA = miRNAs_r_m.miRNA_ID)") 

# merge the targetScan + mirdb query 
merged_ts_mdb = merge(temp_TS, temp_mirdb, by.x = 'miR_target',
                      by.y = 'miR_target', suffixes = c('_ts','_mirdb'),
                      all.x = TRUE, all.y = TRUE)
rm(miRDB,temp_TS,temp_mirdb)
gc()
# clean the data manualy in text editor to remove MTIS with NA/Null score in both dbs
# remove the column not required
colnames(merged_ts_mdb) # column kept in table
#[1] "miRNA"                          "target_gene"     # by spliting the keyword miR_target                   
#[3] "Total_context_score"            "Cum_weighted_score"
#[5] "miRDBscores"                    "miR_target"  
##########################################################################################
merged_ts_mdb = merged_ts_mdb[!names(merged_ts_mdb) %in% 
                                c("MiRBase_ID","Gene_Symbol","miRNA","GeneSymbol")]
#split miR_target to mirna and targetgene
require(tidyr)

merged_ts_mdb = separate(data = merged_ts_mdb, col = miR_target, 
                         into = c("miRNA","target_gene"), sep = " ")

merged_ts_mdb$miR_target = paste(merged_ts_mdb$miRNA,
                                 merged_ts_mdb$target_gene)
merged_ts_mdb[, "miR_target"] = merged_ts_mdb$miR_target

colnames(merged_ts_mdb) # match the sequence of columns with above

#scaling the scores from 0-1 #targetScan
#Check and convert colmn datatypes to numeric
merged_ts_mdb$Cum_weighted_score = as.numeric(merged_ts_mdb$Cum_weighted_score, na.rm = TRUE)
#not using total context score for TArgetsan predictions
merged_ts_mdb$miRDBscores = as.numeric(merged_ts_mdb$miRDBscores, na.rm = TRUE)

ts_max = max(abs(merged_ts_mdb$Cum_weighted_score), na.rm=TRUE)
ts_min = min(abs(merged_ts_mdb$Cum_weighted_score), na.rm=TRUE)
merged_ts_mdb$Cum_wt_score_scaled = (((abs(merged_ts_mdb$Cum_weighted_score))
                                      -ts_min)* (1/(ts_max-ts_min)))
rm(ts_max, ts_min)
########################################################################################
#scaling the scores from 0-1 # mirdb
mirdb_max = max(abs(merged_ts_mdb$miRDBscores), na.rm=TRUE)
mirdb_min = min(abs(merged_ts_mdb$miRDBscores), na.rm=TRUE)
merged_ts_mdb$mirdb_score_scaled = (((abs(merged_ts_mdb$miRDBscores))-mirdb_min)
                                    * (1/(mirdb_max-mirdb_min))
)
rm(mirdb_max, mirdb_min)
##########################################################################################
#get single prediction score from  targetScan and mirdb and scale it
predict_score = function(ts, mdb){
  
  if((is.na(ts)) || (is.na(mdb))) {
    score = sum(ts, mdb, na.rm = TRUE)
  }
  else{
    score = ((ts+mdb)*2)
  }
  return(score)
}

#########calling above function to get predict score

for(i in 1:nrow(merged_ts_mdb)){
  merged_ts_mdb$predict_score[i] = predict_score(abs(merged_ts_mdb$Cum_wt_score_scaled[i]),
                                                 abs(merged_ts_mdb$mirdb_score_scaled[i]))
}

#scaling the predictedscroe from 0-1
predict_max = max(abs(merged_ts_mdb$predict_score), na.rm=TRUE)
predict_min = min(abs(merged_ts_mdb$predict_score), na.rm=TRUE)
merged_ts_mdb$predict_score_scaled = (((abs(merged_ts_mdb$predict_score))-predict_min)
                                      * (1/(predict_max-predict_min))
)
rm(predict_max, predict_min, i)

######### Add validation data ## mtis_validated##########
mtis_validated <- read.delim("D:/mirnas94/validation_data/validationdata_db.txt",
                             stringsAsFactors=FALSE)
#Above file is preprocessed - no ne3ed for following steps => go to query part
# make keyword in validationdata mirna_target
mtis_validated$miR_target = paste(mtis_validated$mirna,
                                  mtis_validated$gene)
mtis_validated[, "miR_target"] = mtis_validated$miR_target

## count number of publication per mti
validation_db = sqldf("select S.*, C.cnt
                      from mtis_validated S
                      inner join (Select miR_target, count(miR_target) as cnt
                      from mtis_validated
                      group by miR_target) C on S.miR_target = C.miR_target")
#group pubmids for 1 mti in 1 colnm
data = validation_db
validation_db = data %>%
  group_by(miR_target,mirna,
           gene, cnt
  ) %>%
  summarize(References = paste0(unique(pubmid), collapse = ";"))
rm(data, mtis_validated)

# QUERY for user required mirna from validation_db
temp_valid = sqldf("select mtis_validated.* from mtis_validated
                   where exists
                   (select miRNAs_r_m.miRNA_ID from miRNAs_r_m where 
                   mtis_validated.mirna = miRNAs_r_m.miRNA_ID)")

# merge valid_temp to merged_ts_mdb db
merged_ts_mdb_valid = merge(merged_ts_mdb, temp_valid, by.x = 'miR_target',
                            by.y = 'miR_target', suffixes = c('.p','.v'), 
                            all.x = TRUE, all.y = TRUE)
rm(temp_valid,mtis_validated)

#drop columns not reqd
merged_ts_mdb_valid = merged_ts_mdb_valid[ , !names(merged_ts_mdb_valid) %in% 
                                             c("mirna","gene", "miRNA", "target_gene")]
#split the miR_target cloumn to get mirnas+genes
merged_ts_mdb_valid = separate(data = merged_ts_mdb_valid, col = miR_target, 
                               into = c("miRNA","target_gene"), sep = " ")
#make keywords
merged_ts_mdb_valid$miR_target = paste(merged_ts_mdb_valid$miRNA,
                                       merged_ts_mdb_valid$target_gene)
merged_ts_mdb_valid[, "miR_target"] = merged_ts_mdb_valid$miR_target


#score the mtis no of validations
validity_score = function(cnt){
  if(is.na(cnt)){score = 0}
  else if(as.numeric(cnt)==0){score = 0}
  else if(as.numeric(cnt)==1){score = 0.5}
  else if(as.numeric(cnt)==2){score = 1.0}
  else if(as.numeric(cnt)>=3){score = 1.5}
  return(score)
}
###########################################################################################
#for(i in 1:nrow(merged_ts_mdb_valid)){                                                   #
#  merged_ts_mdb_valid$Validity_score[i] = validity_score(abs(merged_ts_mdb_valid$cnt[i]))#
#}                                                                                        #
##########################################################################################
# calculating validation score for number of publications for each MTI
#x[10] = merged_ts_mdb_valid$cnt
merged_ts_mdb_valid$Validity_score = apply(merged_ts_mdb_valid,1, function(x){
  return(validity_score(x[10]))
})
##########################################################################################
#for(i in 1:nrow(merged_ts_mdb_valid)){                                                  #
#  merged_ts_mdb_valid$Final_score[i] = sum(merged_ts_mdb_valid$Validity_score[i],       #
#                                           merged_ts_mdb_valid$predict_score_scaled[i], #
#                                           na.rm = TRUE)}                               #
##########################################################################################
# Sum up Validity_score[colmn# 13] + predict_score_scaled[9] to get Final_score[14]

merged_ts_mdb_valid$Final_score = apply(merged_ts_mdb_valid,1, function(x){
  return(sum(as.numeric(x[9])+as.numeric(x[13]),na.rm = TRUE))
})
rm(merged_ts_mdb)

## sepearte databases for rat and mouse
rat_db = sqldf("select * from merged_ts_mdb_valid where 
               merged_ts_mdb_valid.miRNA like '%rno%'")
mouse_db = sqldf("select * from merged_ts_mdb_valid where 
                 merged_ts_mdb_valid.miRNA like '%mmu%'")


##########################################################################################
#uploading rno-miR fold change data of epileptogenesis for dentate gyrus region only
dg_sele_hm <- read.delim("D:/mirnas94/dg_sele_hm.txt", stringsAsFactors=FALSE)
# preprocessd file no need to follow the below steps # goto quant score
#upload Fold change data opf mirnas for Dentate gyrus region 
dg_sele <- read.delim("D:/mirnas94/dg_sele.txt", stringsAsFactors=FALSE)
#convert FC to log2values # using FC_16day readings only - before 1st Seizure
colnames(dg_sele)  = c("rno_mir",	"P-value",	"1hr",	"24hr",	"72hr",	"10day",	"16day",
                       "1stSeiz",	"Chronic",	"Control")
dg_sele[c("log2_FC_1hr","log2_FC_24hr","log2_FC_72hr","log2_FC_10day",
          "log2_FC_16day","log2_FC_1Siez","log2_FC_Chronic")] = NA

dg_sele$log2_FC_16day = log2(as.numeric(as.character(dg_sele$`16day`)))
########################## get fold change values for mmu miRNAs ########################
#don't have mmmu-miRNA expression data yet
#equate same rat FC values if mmu and rno seq are same 
rno_mirSeq = sqldf("select dg_sele.*,mirbase21.Mature_Seq_21 
                   from dg_sele,mirbase21
                   where dg_sele.rno_mir = mirbase21.Mature_ID_21")

mmu_mirna = sqldf("select mirbase21.Mature_ID_21,rno_mirSeq.Mature_Seq_21
                  from rno_mirSeq, mirbase21
                  where (rno_mirSeq.Mature_Seq_21 = mirbase21.Mature_Seq_21) 
                  and mirbase21.Mature_ID_21 like '%mmu%'")
mmu_mirna = mmu_mirna[!duplicated(mmu_mirna),]

dg_sele_hm=merge(rno_mirSeq, mmu_mirna, by.x = 'Mature_Seq_21',
                 by.y = 'Mature_Seq_21', suffixes = c('_rno','_mmu'), 
                 all.x = TRUE, all.y = TRUE)
dg_sele_hm=dg_sele_hm[!duplicated(dg_sele_hm),]
rm(rno_mirSeq,mmu_mirna,dg_sele)
##########################################################################################
