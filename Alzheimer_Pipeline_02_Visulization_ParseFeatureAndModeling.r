##############################
#Created on Tue Feb 23 09:56:23 2016

#@author: Chester
##############################

##############################
# import libraries
##############################
library(hash)
library(sqldf)

##############################
# set parameters
##############################
str_inputFilePath = "D:\\Phd\\Grade_04\\Alzheimer\\Log\\0212"

##############################
# define functions
##############################
# get batchs by scanning the directory of input path
list_inputFilePath = list.files(path=str_inputFilePath, pattern=".csv", full.names = FALSE)
hash_inputFileBatch = hash()
for(idx in 1:length(list_inputFilePath)){
    str_inputFileBatch = paste(strsplit(list_inputFilePath[idx],"_")[[1]][2],strsplit(list_inputFilePath[idx],"_")[[1]][3],sep="_")
    if (!(has.key(str_inputFileBatch,hash_inputFileBatch))){
        if (nchar(toString(strsplit(str_inputFileBatch,"_")[[1]][1])) <= 2){
            hash_inputFileBatch[str_inputFileBatch] <- 1
        }else{
            hash_inputFileBatch[str_inputFileBatch] <- 2
        }
    }
}

list_inputFileBatch = keys(hash_inputFileBatch)
for(idx in 1:length(list_inputFileBatch)){
    if(values(hash_inputFileBatch[list_inputFileBatch[idx]]) == 1){
        # this batch is a classification problem
        str_thisBatch = list_inputFileBatch[idx]
        df_rawDataAll = data.frame(Chr=integer(), GeneSymbol=character(), Filter=character(), Accuracy=double())
        for(idxc in 1:22){
            df_rawData = data.frame(read.csv(paste0(str_inputFilePath, "\\chr", toString(idxc), "_", str_thisBatch)))
            df_rawData = sqldf(paste0("SELECT '", toString(idxc), "' AS chr,* FROM df_rawData"))
            colnames(df_rawData) = c("Chr", "GeneSymbol", "Filter",  "Accuracy")
            df_rawDataAll = rbind(df_rawDataAll, df_rawData)
        }
        df_rawDataAll = df_rawDataAll[complete.cases(df_rawDataAll),]
        df_rawDataAll = sqldf("SELECT * FROM df_rawDataAll WHERE Accuracy>=0")
        df_rawDataAll$Chr = factor(df_rawDataAll$Chr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"))
        
        # draw boxplot
        #=============================
        # 1. ADNI
        str_validationType = "Accuracy"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='ADNI'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_ADNI_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(Accuracy~Chr,data=df_selectedData, main=paste0(str_validationType, " of each Gene with ADNI filter"), xlab="Chromosome No.", ylab=paste0(str_validationType), ylim=c(0.2,0.5))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
              #text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.505, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
        #=============================
        # 2. dbSNP
        str_validationType = "Accuracy"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='dbSNP'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_dbSNP_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(Accuracy~Chr,data=df_selectedData, main=paste0(str_validationType, " of each Gene with dbSNP filter"), xlab="Chromosome No.", ylab=paste0(str_validationType), ylim=c(0.2,0.5))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            #text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.505, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
        #=============================
        # 3. HapMapPhased
        str_validationType = "Accuracy"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='HapMapPhased'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_HapMapPhased_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(Accuracy~Chr,data=df_selectedData, main=paste0(str_validationType, " of each Gene with HapMapPhased filter"), xlab="Chromosome No.", ylab=paste0(str_validationType), ylim=c(0.2,0.5))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            #text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.505, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
    }else{
        # otherwise, the batch is a regression problem
        str_thisBatch = list_inputFileBatch[idx]
        df_rawDataAll = data.frame(Chr=integer(), GeneSymbol=character(), Filter=character(), Pearson=double(), Spearman=double(), AVG_P_S=double())
        for(idxc in 1:22){
            df_rawData = data.frame(read.csv(paste0(str_inputFilePath, "\\chr", toString(idxc), "_", str_thisBatch)))
            df_rawData = sqldf(paste0("SELECT '", toString(idxc), "' AS chr,* FROM df_rawData"))
            colnames(df_rawData) = c("Chr", "GeneSymbol", "Filter",  "Pearson", "Spearman", "AVG_P_S")
            df_rawDataAll = rbind(df_rawDataAll, df_rawData)
        }
        df_rawDataAll = df_rawDataAll[complete.cases(df_rawDataAll),]
        df_rawDataAll = sqldf("SELECT * FROM df_rawDataAll WHERE Pearson>=0 AND Spearman>=0 AND AVG_P_S>=0 ORDER BY Chr")
        df_rawDataAll$Chr = factor(df_rawDataAll$Chr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"))
      
        # draw boxplot
        #=============================
        # 1. ADNI
        str_validationType = "Pearson"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='ADNI'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_ADNI_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(Pearson~Chr,data=df_selectedData, main=paste0(str_validationType, " Cor. of each Gene with ADNI filter"), xlab="Chromosome No.", ylab=paste0(str_validationType, " Cor."), ylim=c(0,0.3))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.305, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
        #-----------------------------
        str_validationType = "Spearman"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='ADNI'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_ADNI_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(Spearman~Chr,data=df_selectedData, main=paste0(str_validationType, " Cor. of each Gene with ADNI filter"), xlab="Chromosome No.", ylab=paste0(str_validationType, " Cor."), ylim=c(0,0.3))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.305, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
        #-----------------------------
        str_validationType = "AVG_P_S"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='ADNI'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_ADNI_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(AVG_P_S~Chr,data=df_selectedData, main=paste0(str_validationType, " Cor. of each Gene with ADNI filter"), xlab="Chromosome No.", ylab=paste0(str_validationType, " Cor."), ylim=c(0,0.3))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.305, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
        #=============================
        # 2. dbSNP
        str_validationType = "Pearson"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='dbSNP'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_dbSNP_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(Pearson~Chr,data=df_selectedData, main=paste0(str_validationType, " Cor. of each Gene with dbSNP filter"), xlab="Chromosome No.", ylab=paste0(str_validationType, " Cor."), ylim=c(0,0.3))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.305, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
        #-----------------------------
        str_validationType = "Spearman"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='dbSNP'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_dbSNP_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(Spearman~Chr,data=df_selectedData, main=paste0(str_validationType, " Cor. of each Gene with dbSNP filter"), xlab="Chromosome No.", ylab=paste0(str_validationType, " Cor."), ylim=c(0,0.3))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.305, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
        #-----------------------------
        str_validationType = "AVG_P_S"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='dbSNP'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_dbSNP_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(AVG_P_S~Chr,data=df_selectedData, main=paste0(str_validationType, " Cor. of each Gene with dbSNP filter"), xlab="Chromosome No.", ylab=paste0(str_validationType, " Cor."), ylim=c(0,0.3))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.305, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
        #=============================
        # 3. HapMapPhased
        str_validationType = "Pearson"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='HapMapPhased'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_HapMapPhased_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(Pearson~Chr,data=df_selectedData, main=paste0(str_validationType, " Cor. of each Gene with HapMapPhased filter"), xlab="Chromosome No.", ylab=paste0(str_validationType, " Cor."), ylim=c(0,0.3))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.305, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
        #-----------------------------
        str_validationType = "Spearman"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='HapMapPhased'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_HapMapPhased_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(Spearman~Chr,data=df_selectedData, main=paste0(str_validationType, " Cor. of each Gene with HapMapPhased filter"), xlab="Chromosome No.", ylab=paste0(str_validationType, " Cor."), ylim=c(0,0.3))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.305, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
        #-----------------------------
        str_validationType = "AVG_P_S"
        df_selectedData = sqldf("SELECT * FROM df_rawDataAll WHERE Filter='HapMapPhased'")
        str_resultOfAPOE = sqldf("SELECT * FROM df_selectedData WHERE GeneSymbol=='APOE'")
        str_countBetterThanAPOE = sqldf(paste0("SELECT COUNT(*) FROM df_selectedData WHERE ", str_validationType, ">", str_resultOfAPOE[str_validationType]))
        png(paste0(str_inputFilePath, "\\boxplot_", tolower(str_validationType), "_HapMapPhased_", gsub(".csv", "", str_thisBatch), ".png"), width = 1024, height = 950, pointsize=15)
        bp_plot = boxplot(AVG_P_S~Chr,data=df_selectedData, main=paste0(str_validationType, " Cor. of each Gene with HapMapPhased filter"), xlab="Chromosome No.", ylab=paste0(str_validationType, " Cor."), ylim=c(0,0.3))
        if (length(bp_plot$out) > 0){
            vector_labelOfOutlier = c()
            for(idxo in 1:length(bp_plot$out)){
                int_idxoOfTarget = as.integer(rownames(df_selectedData)[which(df_selectedData["Chr"] == bp_plot$group[idxo] & df_selectedData[str_validationType] == bp_plot$out[idxo])])
                vector_labelOfOutlier = c(vector_labelOfOutlier, ifelse(df_selectedData[int_idxoOfTarget,str_validationType]>=as.numeric(str_resultOfAPOE[str_validationType]),toString(df_selectedData[int_idxoOfTarget,"GeneSymbol"]),""))
            }
            text(bp_plot$group, bp_plot$out, vector_labelOfOutlier, pos=4, col="darkred")
        }
        text(0.1, as.numeric(str_resultOfAPOE[str_validationType])+0.005, paste0(str_validationType, " of APOE: ",toString(round(str_resultOfAPOE[str_validationType],4))), pos=4)
        abline(h = str_resultOfAPOE[str_validationType], col = "gray60")
        legend(0.1, 0.305, paste0("Select ", str_countBetterThanAPOE, " Genes    "))
        dev.off()
    }
}






