library(randomForest)
library(caret)
library(e1071)
library(caTools)

setwd("/Users/adigenova/Projects/IARC/bioinfo/mesomics/SVs/MERGE/random_forest/svaba")
svaba_sv=read.table("training/svaba-training-rf.txt",h=T)
#we convert some values of the
#note on reads_supporting variants calls
#https://github.com/walaj/svaba/issues/33

svaba_sv$RAF<-ifelse(svaba_sv$RAF>1, 1, svaba_sv$RAF)
sapply(svaba_sv, class)
#filter variants matching PON or GNOMAD

svaba_so=svaba_sv[svaba_sv$SOMATIC >0,]
svaba_germ=svaba_sv[svaba_sv$SOMATIC ==0,]
  
#pick useful values for random forest
cmodel=c("PON_BC1","PON_BC2","GNOMAD_AC",
           "GNOMAD_BC1","GNOMAD_BC2",
           "PCAWG_BC1","PCAWG_BC2","SOMATIC","COSMIC_GENE",
           "CENTROMER","EXON","CONSERVATION_BC1",
           "CONSERVATION_BC2","CNV_TCN1", "CNV_TCN2", 
           "CNV_CF1", "CNV_CF2", "SVLEN","RAF","RFS")

  
svaba_so=svaba_so[,cmodel]
svaba_germ=svaba_germ[,cmodel]
svaba_so$SOMATIC=1
svaba_so$SOMATIC=as.factor(svaba_so$SOMATIC)
svaba_germ$SOMATIC=as.factor(svaba_germ$SOMATIC)
  
  #we split the data into trainin and test
d_germ=svaba_germ[sample(nrow(svaba_germ), dim(svaba_so)*1),]
d_soma=svaba_so[sample(nrow(svaba_so), dim(svaba_so)*1),]
dim(d_germ)
dim(d_soma)
svaba_data=rbind(d_soma,d_germ)
head(svaba_data)
svaba_data$COSMIC_GENE=as.factor(svaba_data$COSMIC_GENE)                   
svaba_data$CENTROMER=as.factor(svaba_data$CENTROMER) 
svaba_data$EXON=as.factor(svaba_data$EXON) 
sapply(svaba_data, class)
#we drop unused columnns (vars)
#df <- subset(svaba_data, select = -c(FILE, CHROM, POS,CODING))
df <- svaba_data
#we split the data for testing and evaluation
sample = sample.split(df$SOMATIC, SplitRatio = 0.75)
train = subset(df, sample == TRUE) # 75% training
test  = subset(df, sample == FALSE) # 25% for evaluation
# fitting datasset to avoid overfitting.
dim(train)
dim(test)
#test with regular training 
rf <- randomForest(
  SOMATIC ~ .,
  #importance=T,
  data=train,
  keep.forest = TRUE
)
#we make the predictions
pred = predict(rf, newdata=test)
#we evaluate the prediction
confusionMatrix(pred, test$SOMATIC)
#we plot the relevant of variables
varImp(rf)
varImpPlot(rf,main="SVaba")

  #train the forest with the whole data 1:1
rf <- randomForest(
  SOMATIC ~ .,
  #importance=T,
  data=df,
  keep.forest = TRUE
)

cat('Model building complete. Saving model...\n')
saveRDS(rf, file="svaba_rf_model_v1.rds")
cat('Model saved.\n')
save.image(file=paste0("svaba_training_v1",".RData"))

rf <- readRDS("svaba_rf_model_v1.rds")    
  
#we load data to use the model and predict
svaba_tonly=read.table("prediction/svaba_tonly_pred.txt",h=T)
#https://github.com/walaj/svaba/issues/33
#the above affect 1/1 genotypes, usually related to homozygous variants.
svaba_tonly$RAF<-ifelse(svaba_tonly$RAF>1, 1, svaba_tonly$RAF)
svaba_tonly$COSMIC_GENE=as.factor(svaba_tonly$COSMIC_GENE)                   
svaba_tonly$CENTROMER=as.factor(svaba_tonly$CENTROMER) 
svaba_tonly$EXON=as.factor(svaba_tonly$EXON) 
svaba_pto=svaba_tonly[,cmodel]  
svaba_pto$SOMATIC=0
svaba_pto$SOMATIC=as.factor(svaba_pto$SOMATIC)  
sample.features<-subset(svaba_pto, select = names(rf$forest$xlevels))
any.missing.features <- which(!names(rf$forest$xlevels) %in% colnames(sample.features));
#we predict using the model
predicted.class <- predict(rf, sample.features) 
predicted.prob <- predict(rf, sample.features, type="prob")
predictions <- data.frame(
  svaba_tonly,
  Predicted.Class = predicted.class,
  Predicted.prob = predicted.prob
)  

#black listed zones for SVaba, the others callers use a similar approach

#no match PON, hit an exon , read support > 5,and have somatic prob >0.5 & 
ea=predictions[predictions$PON==0 & predictions$Predicted.prob.1 > 0.5 & predictions$PE_SR >=5 & predictions$EXON==1 & predictions$RAF < 1 & predictions$CENTROMER == 0,]
#no match PON, do not hit an exon, read support > 5 , and have somatic prob > 0.75
ia=predictions[predictions$PON==0 & predictions$Predicted.prob.1 > 0.75 & predictions$PE_SR >=5 & predictions$EXON==0 & predictions$RAF < 1 & predictions$CENTROMER == 0,]
a=rbind(ea,ia)



#sort(summary(a$TUMOR_ID))

#sort(summary(a$Tumor_ID))
write.table(
  x = predictions,
  file = paste0("svaba", "_predictions_unfiltered.txt"),
  sep = '\t',
  row.names = TRUE,
  quote = FALSE
);

#write the filtered predictions
write.table(
  x = a,
  file = paste0("svaba", "_predictions_filtered.txt"),
  sep = '\t',
  row.names = TRUE,
  quote = FALSE
);

#we save the sesion information
writeLines(capture.output(sessionInfo()), paste("svaba_R_sessionInfo.txt",sep=""))
