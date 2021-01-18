library(randomForest)
library(caret)
library(e1071)
library(caTools)

setwd("/Users/adigenova/Projects/IARC/bioinfo/mesomics/SVs/MERGE/random_forest/manta/")
manta_sv=read.table("training/manta_training_v1.txt",h=T)
sapply(manta_sv, class)
#filter variants matching PON or GNOMAD

manta_so=manta_sv[manta_sv$SOMATIC >0,]
manta_germ=manta_sv[manta_sv$SOMATIC ==0,]
  
#pick useful values for random forest
cmodel=c("PON_BC1","PON_BC2","GNOMAD_AC",
           "GNOMAD_BC1","GNOMAD_BC2",
           "PCAWG_BC1","PCAWG_BC2","SOMATIC","COSMIC_GENE",
           "CENTROMER","EXON","CONSERVATION_BC1",
           "CONSERVATION_BC2","CNV_TCN1", "CNV_TCN2", 
           "CNV_CF1", "CNV_CF2", "SVLEN","RAF","RFS")
  
manta_so=manta_so[,cmodel]
manta_germ=manta_germ[,cmodel]
manta_so$SOMATIC=1
manta_so$SOMATIC=as.factor(manta_so$SOMATIC)
manta_germ$SOMATIC=as.factor(manta_germ$SOMATIC)
  
  #we split the data into trainin and test
d_germ=manta_germ[sample(nrow(manta_germ), dim(manta_so)*1),]
d_soma=manta_so[sample(nrow(manta_so), dim(manta_so)*1),]
dim(d_germ)
dim(d_soma)
manta_data=rbind(d_soma,d_germ)
head(manta_data)
manta_data$COSMIC_GENE=as.factor(manta_data$COSMIC_GENE)                   
manta_data$CENTROMER=as.factor(manta_data$CENTROMER) 
manta_data$EXON=as.factor(manta_data$EXON) 
sapply(manta_data, class)
#we drop unused columnns (vars)
#df <- subset(manta_data, select = -c(FILE, CHROM, POS,CODING))
df <- manta_data
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
varImp(rf)
#we plot the relevant of variables
varImpPlot(rf,main="Manta")
  #A second random forest without matching PON 
  #train the forest with the whole data 1:1
rf <- randomForest(
  SOMATIC ~ .,
  #importance=T,
  data=df,
  keep.forest = TRUE
)

cat('Model building complete. Saving model...\n')
saveRDS(rf, file="manta_rf_model_v1.rds")
cat('Model saved.\n')
save.image(file=paste0("manta_training_v1",".RData"))

rf <- readRDS("manta_rf_model_v1.rds")  

#we load data to use the model and predict
manta_tonly=read.table("prediction/manta_prediction_v1.txt",h=T)
manta_tonly$COSMIC_GENE=as.factor(manta_tonly$COSMIC_GENE)                   
manta_tonly$CENTROMER=as.factor(manta_tonly$CENTROMER) 
manta_tonly$EXON=as.factor(manta_tonly$EXON) 
manta_pto=manta_tonly[,cmodel]  
manta_pto$SOMATIC=0
manta_pto$SOMATIC=as.factor(manta_pto$SOMATIC)  
sample.features<-subset(manta_pto, select = names(rf$forest$xlevels))
any.missing.features <- which(!names(rf$forest$xlevels) %in% colnames(sample.features));
#we predict using the model
predicted.class <- predict(rf, sample.features) 
predicted.prob <- predict(rf, sample.features, type="prob")
predictions <- data.frame(
  manta_tonly,
  Predicted.Class = predicted.class,
  Predicted.prob = predicted.prob
)  

#no match PON, hit an exon , read support > 5,and have somatic prob >0.5
ea=predictions[predictions$PON==0 & predictions$Predicted.prob.1 > 0.5 & predictions$PE_SR >=5 & predictions$EXON==1,]
#no match PON, do not hit an exon, read support > 5 , and have somatic prob > 0.75
ia=predictions[predictions$PON==0 & predictions$Predicted.prob.1 > 0.75 & predictions$PE_SR >=5 & predictions$EXON==0,]
a=rbind(ea,ia)
#sort(summary(a$TUMOR_ID))

#sort(summary(a$Tumor_ID))
write.table(
  x = predictions,
  file = paste0("manta", "_predictions_unfiltered.txt"),
  sep = '\t',
  row.names = TRUE,
  quote = FALSE
);

#write the filtered predictions
write.table(
  x = a,
  file = paste0("manta", "_predictions_filtered.txt"),
  sep = '\t',
  row.names = TRUE,
  quote = FALSE
);


#we save the sesion information
writeLines(capture.output(sessionInfo()), paste("manta_R_sessionInfo.txt",sep=""))

