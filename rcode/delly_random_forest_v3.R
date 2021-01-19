library(randomForest)
library(caret)
library(e1071)
library(caTools)

setwd("/Users/adigenova/Projects/IARC/bioinfo/mesomics/SVs/MERGE/random_forest/delly")
delly_sv=read.table("raw_data3/rf_delly_v3.matrix.txt",h=T)
sapply(delly_sv, class)
#filter variants matching PON or GNOMAD
#delly_sv=delly_sv[delly_sv$PON==0,]
#delly_sv=delly_sv[delly_sv$GNOMAD==0,]
#we split somatic and germline variants
#delly_so=delly_sv[delly_sv$SOMATIC >0 & delly_sv$PON ==0,]
#delly_germ=delly_sv[delly_sv$SOMATIC ==0 & delly_sv$PON ==0,]

delly_so=delly_sv[delly_sv$SOMATIC >0,]
delly_germ=delly_sv[delly_sv$SOMATIC ==0,]
  
  #pick useful values for random forest
  cmodel=c("PON_BC1","PON_BC2","GNOMAD_AC",
           "GNOMAD_BC1","GNOMAD_BC2",
           "PCAWG_BC1","PCAWG_BC2","SOMATIC","COSMIC_GENE",
           "CENTROMER","EXON","CONSERVATION_BC1",
           "CONSERVATION_BC2","CNV_TCN1", "CNV_TCN2", 
           "CNV_CF1", "CNV_CF2", "SVLEN","RAF","RFS")
  
  delly_so=delly_so[,cmodel]
  delly_germ=delly_germ[,cmodel]
  delly_so$SOMATIC=1
  delly_so$SOMATIC=as.factor(delly_so$SOMATIC)
  delly_germ$SOMATIC=as.factor(delly_germ$SOMATIC)
  
  #we split the data into trainin and test
  d_germ=delly_germ[sample(nrow(delly_germ), dim(delly_so)*1),]
  d_soma=delly_so[sample(nrow(delly_so), dim(delly_so)*1),]
  dim(d_germ)
  dim(d_soma)
  delly_data=rbind(d_soma,d_germ)
  head(delly_data)
  delly_data$COSMIC_GENE=as.factor(delly_data$COSMIC_GENE)                   
  delly_data$CENTROMER=as.factor(delly_data$CENTROMER) 
  delly_data$EXON=as.factor(delly_data$EXON) 
  sapply(delly_data, class)
  #we drop unused columnns (vars)
  #df <- subset(delly_data, select = -c(FILE, CHROM, POS,CODING))
  df <- delly_data
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
  varImpPlot(rf,main="Delly")
  #A second random forest without matching PON 
  #train the forest with the whole data 1:1
  rf <- randomForest(
    SOMATIC ~ .,
    #importance=T,
    data=df,
    keep.forest = TRUE
  )
  cat('Model building complete. Saving model...\n')
  saveRDS(rf, file="delly_rf_model_v3.rds")
  cat('Model saved.\n')
  save.image(file=paste0("delly_training_v3",".RData"))
  
#we load the trained model to made the classification
rf <- readRDS("delly_rf_model_v3.rds")

#we load data to use the model and predict
delly_tonly=read.table("t-only/t-only-delly.txt",h=T)
delly_tonly$COSMIC_GENE=as.factor(delly_tonly$COSMIC_GENE)                   
delly_tonly$CENTROMER=as.factor(delly_tonly$CENTROMER) 
delly_tonly$EXON=as.factor(delly_tonly$EXON) 
delly_pto=delly_tonly[,cmodel]  
delly_pto$SOMATIC=0
delly_pto$SOMATIC=as.factor(delly_pto$SOMATIC)  
sample.features<-subset(delly_pto, select = names(rf$forest$xlevels))
any.missing.features <- which(!names(rf$forest$xlevels) %in% colnames(sample.features));
#we predict using the model
predicted.class <- predict(rf, sample.features) 
predicted.prob <- predict(rf, sample.features, type="prob")
predictions <- data.frame(
  delly_tonly,
  Predicted.Class = predicted.class,
  Predicted.prob = predicted.prob
)  
  
#no match PON, hit an exon , read support > 5,and have somatic prob >0.5
ea=predictions[predictions$PON==0 & predictions$Predicted.prob.1 > 0.5 & predictions$PE_SR >=5 & predictions$EXON==1 & predictions$CENTROMER == 0,]
#no match PON, do not hit an exon, read support > 5 , and have somatic prob > 0.75
ia=predictions[predictions$PON==0 & predictions$Predicted.prob.1 > 0.75 & predictions$PE_SR >=5 & predictions$EXON==0 & predictions$CENTROMER == 0,]
a=rbind(ea,ia)
#sort(summary(a$Tumor_ID))
write.table(
  x = predictions,
  file = paste0("delly", "_predictions_unfiltered.txt"),
  sep = '\t',
  row.names = TRUE,
  quote = FALSE
);

#write the filtered predictions
write.table(
  x = a,
  file = paste0("delly", "_predictions_filtered.txt"),
  sep = '\t',
  row.names = TRUE,
  quote = FALSE
);


#we save the sesion information
writeLines(capture.output(sessionInfo()), paste("delly_R_sessionInfo.txt",sep=""))