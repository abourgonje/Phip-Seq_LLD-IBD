# =========================================
# By: R.Gacesa (UMCG, 2021)
# 
# prediction of UC, CD or CD vs UC 
# using antibody panels
# =========================================
#
# load libraries
# =========================================
library(caret)
library(pROC)
library(plotROC)
library(foreach)
library(plyr)
library(gbm)
library(doSNOW)

# set WD and load helper functions
# =========================================
setwd('D:/UMCG/Arno_Antibody_prediction/')
source('R_ML_helperscripts.R')
source('R_ML_scripts_v4.R')

# prep paralelization
registerDoSNOW(makeCluster(12, type = "SOCK"))

# PREDICTIVE MODELLING RUN FOR ANTIBODIES
# =========================================================
# PART I: training (on training set), test set validation, 
#   segata data external (negative) validation
# =========================================================

# RUN PARAMETERS
# =========================================================
# runNameMain has to be set to one of:
#  - "CD" for Crohn's disease prediction
#  - "UC" for Ulcerative colitis prediction
#  - "CD_vs_UC" for separation between CD and UC
 
runNameMain <- "CD"
pC <- "Y"; if (runNameMain == "CD_vs_UC") {pC <- "CD"}

# mlMethods: caret-implemented ML algorithms to use 
#   NOTE: code is implemented & tested for glmnet, gbm, avNNet, svmRadial
#         and might not work with other algorithms
mlMethods <- c("glmnet", "gbm", "avNNet", "svmRadial")

# datasets to consider (all = all antibodies, agilent & twist are sub-sets)
dataSubsets <- c("all","agilent","twist")

# training set parameters
trPerc <- 0.80 # proportion of data used for training

# feature pre-selection parameters (will be calculated from training set)
pVF <- 0.005   # feature selection p-value cutoff
minP <- 0.01   # minimal presence below which features are removed
maxP <- 0.99   # maximal presence above which features are removed

# colors for plots
if (runNameMain == "CD") {
  myCol = "red3"; myCol2 = "red3"; myCol3 = "red1"; myCol4="red4"
} else if (runNameMain == "UC") {
  myCol = "blue3"; myCol2 = "blue3"; myCol3 = "blue1"; myCol4 = "blue4"
} else if (runNameMain == "CD_vs_UC") {
  myCol = "purple3"; myCol2 = "purple3"; myCol3 = "purple1"; myCol4 = "purple4"
}

# MODEL TRAINING AND TESTING
# ==========================================
# loop iterates over algorithms and data-subsets
for (dType in dataSubsets) {
  for (mlMethod in mlMethods) {
    # set run name
    runName <- paste0(runNameMain,'_base_',dType,'_',mlMethod)
    print(paste0(' >> STARTING ML RUN FOR ',runName))
    
    # LOAD DATA
    # =====================
    #  main data
    inDF <- read.table(paste0('Datasets/',runNameMain,'_matched_antibodies.csv'),sep=',',header=T,row.names = 1,stringsAsFactors = F)
    #  negative external test-set
    inDFs <- read.table('Datasets/Israel_cohort_filtered_prevalences.csv',sep=',',header=T)
    
    # DATA PREP    
    # ==============================================================================
    # i) basic prep (variable types, subset antibody panel if needed)
    inDFpr <- prepAbData(inDF,subSet=dType)
    
    # ii) prep training & test sets
    # ==============================================================================
    set.seed(123897) # fixed seed ensures that all training / test splits are identical
    inTrain <- createDataPartition(y=inDFpr[["Cohort"]],p=trPerc,list=F)
    inDFtr <- inDFpr[inTrain,]
    inDFtst  <- inDFpr[-inTrain,]
    
    # iii) pre-process (using only training set to avoid data leakage)
    # ==============================================================================
    pt <- prepProcessAbData(inDFtr=inDFtr,pVF = pVF,maxP = maxP,minP=minP)
    inDFtrpp <- pt[[2]]; inDFprep <- pt[[1]]
    
    #  > save datasets
    write.table(inDFtrpp,paste0('Model_data_',runNameMain,'/',runName,'_training.csv'),sep=',',row.names = F)
    write.table(inDFtst,paste0('Model_data_',runNameMain,'/',runName,'_test.csv'),sep=',',row.names = F)
    saveRDS(inDFprep,paste0('Model_data_',runNameMain,'/',runName,'_preprocessor.RDS'))
    
    # MODEL TRAINING
    # ==============================================================================
    # iv) model training
    # ==============================================================================
    #  > set model training scheme
    trainC <- trainControl(method="repeatedcv",number=5,repeats = 5,savePredictions = T,classProbs = T,allowParallel = T,
                           verboseIter = F,returnData = F,preProcOptions = NULL,trim = T)
    set.seed(123899) # seed is fixed for reproducibility / consistency [this is not necessary]
    mdlFit <- caret::train(Cohort ~ ., data=inDFtrpp, method = mlMethod, metric = "Kappa", trControl = trainC, 
                           tuneLength=10) #  tune length size defines size of optimization grid
    
    # v) report performance on training set / optimization using cross-validation
    # =============================================================================
    #  > confusion matrix
    mdlMetrics <- getMdlFitXVresults(mdlFit=mdlFit,posClass = pC,mdName = runName)
    #  > ROC
    mdlROC <- compareModelsTrainingCV(fittedMdls = list(mdlFit),
                                      modelNames = c(runName),roc.conf.boot = 100,
                                      posClass = pC,annotateAUConly = T,roc.conf = 0.95,
                                      tit = paste0(runName, " model"),roc.smooth = F,
                                      annotS = 5,diagonalLine = F,textOffSetY = +0.05,textOffSetX = -0.195)
    #   - re-style ROC plot for publication
    mdlRocS <- mdlROC + scale_color_manual(values=c(myCol)) + theme_classic() +
      scale_fill_manual(values = c(myCol2) ) + theme(legend.position = "none") + ggtitle("") + 
      geom_abline(intercept = c(1), slope = 1,col="darkblue",linetype="longdash",size=1.05) + 
      coord_cartesian(xlim=c(1.005,-0.00),ylim=c(0,1.02),expand = F) +
      theme(axis.line = element_line(size = 1.05),axis.ticks = element_line(size = 0.9))
    #print(mdlRocS) # output on screen for debugging
    
    #  > extract model betas (GLM) / variable importance (non-GLM models)
    varImpTable <- getVarImpTbl(mdlFit)
    
    #   >> save model, metrics, ROC, variable importance
    saveRDS(mdlFit,paste0('Model_data_',runNameMain,'/',runName,'_ModelFit.RDS'))
    write.table(mdlMetrics,paste0('Model_metrics_',runNameMain,'/',runName,'_resultsXV_metrics.csv'),sep=',',row.names = F)
    ggsave(plot = mdlRocS,filename = paste0('Model_metrics_',runNameMain,'/',runName,'_resultsXV_ROC.png'),width = 6,height = 6,units = "in",dpi = 320)
    write.table(varImpTable,paste0('Model_metrics_',runNameMain,'/',runName,'_variables.csv'),sep=',',row.names = F)
    
    # MODEL TEST (using Test set and external negative control)
    # =============================================================================
    # v) model test on test set
    # =============================================================================
    #  > preprocess test set (using pre-processing scheme from training set)
    inDFtstpp <- predict(inDFprep,inDFtst)
    
    #  > predict test set and report prediction metrics
    mdlMetricsTest <- getMdlTestResults(mdlFit=mdlFit,mdName = runName,testSet=inDFtstpp,posClass = pC)
    
    #  > generate ROC curve
    mdlRocTest <- compareMdlsDatasets(mdls = list(mdlFit),
                                      dataSets = list(inDFtstpp),
                                      posClass = pC,
                                      mdNames = c(runName),
                                      response = "Cohort",
                                      removeLegend = T,
                                      roc.conf.boot = 100,roc.conf = 0.95,roc.smooth = F,
                                      tit = paste0(runName, " model"),
                                      annotateROC = T,
                                      annotateROCbaseName = F,
                                      annotS = 5,
                                      diagonalLine = F,
                                      textOffSetX = -0.195,
                                      textOffSetY = +0.05)[[1]]
    #print(mdlRocTest) # debug
    #  - ROC styling for publication
    mdlRocTestS <- mdlRocTest + scale_color_manual(values=c(myCol)) + theme_classic() +
      scale_fill_manual(values = c(myCol2) ) + theme(legend.position = "none") + ggtitle("") + 
      geom_abline(intercept = c(1), slope = 1,col="darkblue",linetype="longdash",size=1.05) + 
      coord_cartesian(xlim=c(1.005,-0.00),ylim=c(0,1.02),expand = F) +
      theme(axis.line = element_line(size = 1.05),axis.ticks = element_line(size = 0.9))
    #print(mdlRocTestS) # debug
    # >> save results
    write.table(mdlMetricsTest,paste0('Model_metrics_',runNameMain,'/',runName,'_resultsTest_metrics.csv'),sep=',',row.names = F)
    ggsave(plot = mdlRocTestS,filename = paste0('Model_metrics_',runNameMain,'/',runName,'_resultsTest_ROC.png'),width = 6,height = 6,units = "in",dpi = 320)
    
    # vi) model test on external negative dataset
    # ================================================
    # - note: it includes ONLY NEGATIVES, so we cannot make ROC curves
    #         and not all prediction metrics can be calculated
    if (runNameMain != "CD_vs_UC") {
      inDFs$filename <- NULL
      inDFs$Cohort <- "N"
      #   - merge with our test set positive cases
      inDFstst <- rbind.fill(inDFs,inDFtst[inDFtst$Cohort=="Y",])
      inDFstst[is.na(inDFstst)] <- 0
      inDFstst$Cohort <- as.factor(inDFstst$Cohort)
      inDFststneg <- inDFstst[inDFstst$Cohort == "N",]
      #   - apply pre-processing scheme
      inDFststnegpp <- predict(inDFprep,inDFststneg)
      #   - predict
      mdlMetricsTestExt <- getMdlTestResults(mdlFit=mdlFit,mdName = runName,testSet=inDFststnegpp,dataSetName="Test.set.externalneg")
      #   > save results of testing on external set
      write.table(mdlMetricsTestExt,paste0('Model_metrics_',runNameMain,'/',runName,'_resultsTestExt_metrics.csv'),sep=',',row.names = F)
    }
    print(paste0(' >> DONE WITH ML RUN FOR ',runName))
  }
}

# ==============================================================================
#  DATA COLLECTION
# ==============================================================================
#   - code collects results of individual models and puts it into one table and
#     moves main results to separate folder
# ==============================================================================
runs <- c("CD","UC","CD_vs_UC")

# make folder for results if it does not exist
if (!dir.exists('Results_ROCs_main')) {dir.create('Results_ROCs_main')}

# collect data from each run
for (run in runs) {
  fldr <- paste0('Model_metrics_',run)
  for (dt in c("all","agilent","twist")) {
    res <- NULL
    toMerge <- list.files(pattern = paste0('.*base_',dt,'_.*_metrics.csv'), fldr)
    for (f in toMerge) {
      res <- rbind.data.frame(res,read.table(paste0(fldr,'/',f),sep=',',header=T))
    }
    write.table(res,paste0('Results_merged/Model_results_',run,'_',dt,'.csv'),sep=',',row.names = F)
  }
}
file.copy(from = 'Model_metrics_CD/CD_base_all_glmnet_resultsTest_ROC.png',to = 'Results_ROCs_main/',copy.mode = TRUE,overwrite = T)
file.copy(from = 'Model_metrics_CD/CD_base_all_glmnet_resultsXV_ROC.png',to = 'Results_ROCs_main/',copy.mode = TRUE,overwrite = T)
file.copy(from = 'Model_metrics_UC/UC_base_all_glmnet_resultsTest_ROC.png',to = 'Results_ROCs_main/',copy.mode = TRUE,overwrite = T)
file.copy(from = 'Model_metrics_UC/UC_base_all_glmnet_resultsXV_ROC.png',to = 'Results_ROCs_main/',copy.mode = TRUE,overwrite = T)
file.copy(from = 'Model_metrics_CD_vs_UC/CD_vs_UC_base_all_glmnet_resultsTest_ROC.png',to = 'Results_ROCs_main/',copy.mode = TRUE,overwrite = T)
file.copy(from = 'Model_metrics_CD_vs_UC/CD_vs_UC_base_all_glmnet_resultsXV_ROC.png',to = 'Results_ROCs_main/',copy.mode = TRUE,overwrite = T)

# ==============================================================================
# MODEL OPTIMIZATION using recursive feature selection (RFE)
# ==============================================================================

# which models to optimize?
dType <- "all" # dataset (all, twist or agilent)
mlMethod <- "glmnet" # algorithm

for (runNameMain in c("CD","UC","CD_vs_UC"))  {
  # set colors for plots
  if (runNameMain == "CD") {
    myCol = "red3"; myCol2 = "red3"; myCol3 = "red1"; myCol4="red4"; pC = "Y"
  } else if (runNameMain == "UC") {
    myCol = "blue3"; myCol2 = "blue3"; myCol3 = "blue1"; myCol4 = "blue4"; pC = "Y"
  } else if (runNameMain == "CD_vs_UC") {
    myCol = "purple3"; myCol2 = "purple3"; myCol3 = "purple1"; myCol4 = "purple4"; pC = "CD"
  }
  # set run names
  runNameL <- paste0(runNameMain,'_base_',dType,'_',mlMethod)
  runNameS <- paste0(runNameMain,'_opt_',dType,'_',mlMethod)
  print(paste0(' >> STARTING OPTIMIZATION RUN FOR ',runNameL))
  
  # LOAD DATASETS (Tr, Test, Ext, pre-processor)
  # ===========================================
  # training
  inDFtrpp <- read.table(paste0('Model_data_',runNameMain,'/',runNameL,'_training.csv'),sep=',',header = T)
  inDFtrpp$Cohort <- as.factor(inDFtrpp$Cohort)
  # test
  inDFtst <- read.table(paste0('Model_data_',runNameMain,'/',runNameL,'_test.csv'),sep=',',header=T)
  inDFtst$Cohort <- as.factor(inDFtst$Cohort)
  inDFprep <- readRDS(paste0('Model_data_',runNameMain,'/',runNameL,'_preprocessor.RDS'))
  inDFtstpp <- predict(inDFprep,inDFtst)
  # external
  if (runNameMain != "CD_vs_UC") {
    inDFs <- read.table('Datasets/Israel_cohort_filtered_prevalences.csv',sep=',',header=T)
    #inDFs$Cohort <- as.factor(inDFs$Cohort)
    inDFs$filename <- NULL
    inDFs$Cohort <- "N"
    #   > merge with our test set positive cases
    inDFstst <- rbind.fill(inDFs,inDFtst[inDFtst$Cohort=="Y",])
    inDFstst[is.na(inDFstst)] <- 0
    inDFstst$Cohort <- as.factor(inDFstst$Cohort)
    inDFststneg <- inDFstst[inDFstst$Cohort == "N",]
    #   - apply pre-processing scheme
    inDFststnegpp <- predict(inDFprep,inDFststneg)
  }
  # run RFE if not already done
  # ==========================================
  if (!file.exists(paste0('Model_RFE_',runNameMain,'/',runNameS,'_RFE.RDS'))) {
    # define RFE setup
    rfeMethod <- mlMethod
    # train controller for RFE model fit
    trC <- trainControl(method="repeatedcv",
                        repeats=1,
                        number=5,
                        savePredictions = T,
                        classProbs = T,
                        allowParallel = T)
    # train controller for RFE algorithm
    rfeCtrl <- rfeControl(functions = rfFuncs,
                          method = "repeatedcv",
                          repeats = 5,
                          number=50,
                          verbose = F,
                          allowParallel = T)
    # set steps to test
    if (ncol(inDFtrpp) <= 100) {szs=c(seq(1,ncol(inDFtrpp)-1))
    } else if (ncol(inDFtrpp) <= 210) {szs=c(seq(1,50,1),seq(50,100,1),seq(105,ncol(inDFtrpp)-1,5))
    } else if (ncol(inDFtrpp) > 210) {szs=c(seq(1,50,1),seq(50,100,1),seq(105,200,5), seq(200,ncol(inDFtrpp)-1,10))  }
    szs <- unique(szs)
    # run RFE
    print ('  >> doing RFE profile'); time1 <- Sys.time()
    rfeProfile <- rfe(x=inDFtrpp[,-grep("Cohort",colnames(inDFtrpp))],
                      y=inDFtrpp[["Cohort"]],
                      sizes=szs,
                      rfeControl = rfeCtrl,
                      metric="Kappa",
                      method=rfeMethod,
                      maximize = T,
                      trControl = trC)
    time2 <- Sys.time(); print ('    >>> DONE!'); print(time2 - time1)
    # save it
    saveRDS(rfeProfile,file=paste0('Model_RFE_',runNameMain,'/',runNameS,'_RFE.RDS'))
  }
  
  # (re)load RFE profile
  rfeProfile <- readRDS(paste0('Model_RFE_',runNameMain,'/',runNameS,'_RFE.RDS'))
  
  # plot it
  nMax <- pickSizeTolerance(rfeProfile$results, metric="Kappa",maximize = T,tol=0)
  nT5 <- pickSizeTolerance(rfeProfile$results, metric="Kappa",maximize = T,tol=5)
  nT10 <- pickSizeTolerance(rfeProfile$results, metric="Kappa",maximize = T,tol=10) 
  nT20 <- pickSizeTolerance(rfeProfile$results, metric="Kappa",maximize = T,tol=20) 
  rfeplot <- ggplot(rfeProfile$results) +
    aes(x=Variables,y=Kappa) + 
    geom_line(col=myCol,size=1.05) + geom_point(col=myCol2) +
    geom_errorbar(aes(ymax = Kappa+KappaSD,ymin=Kappa-KappaSD),alpha=0.15,col=myCol2) +
    ggtitle(paste0("Recursive feature elimination (",runNameMain,")")) + 
    theme_classic() + xlab("Number of Variables") + 
    theme(axis.line = element_line(size = 1.05),axis.ticks = element_line(size = 0.9)) + ggtitle("") + 
    xlab('Number of antibody-bound peptides') + ylab("Cohen's Kappa")
      
  #print(rfeplot) # debug
  # save it
  ggsave(plot = rfeplot,filename = paste0('Model_RFE_',runNameMain,'/',runNameS,'_RFEplot.png'),dpi = 600,width = 6,height = 4,scale = 1)
    
  # # REFIT OPTIMIZED MODEL {5 variables}
  # # ===========================================
  # variables
  varT <- rfeProfile$optVariables
  # - keep only optimized vars
  inDFtrppp5 <- inDFtrpp[,colnames(inDFtrpp) %in% c("Cohort",varT[1:5])]
  inDFtrppp10 <- inDFtrpp[,colnames(inDFtrpp) %in% c("Cohort",varT[1:10])]
  # # - train it
  trainC <- trainControl(method="repeatedcv",number=5,repeats = 10,savePredictions = T,classProbs = T,allowParallel = T,
                         verboseIter = F,returnData = F,preProcOptions = NULL,trim = T)
  set.seed(123899)
  mdlFitOpt5 <- train(Cohort ~ ., data=inDFtrppp5, method = mlMethod, metric = "Kappa", trControl = trainC, 
                         tuneLength=20)
  mdlFitOpt10 <- train(Cohort ~ ., data=inDFtrppp10, method = mlMethod, metric = "Kappa", trControl = trainC, 
                               tuneLength=20)
  # > save variable betas
  varImpTable5 <- getVarImpTbl(mdlFitOpt5)
  varImpTable10 <- getVarImpTbl(mdlFitOpt10)
  write.table(varImpTable5,paste0('Model_RFE_',runNameMain,'/',runNameS,'_model5_betas.csv'),sep=',',row.names = F)
  write.table(varImpTable10,paste0('Model_RFE_',runNameMain,'/',runNameS,'_model10_betas.csv'),sep=',',row.names = F)
    
  # > report performance on XV set
  mdlMetricsXV5 <- getMdlFitXVresults(mdlFit=mdlFitOpt5,mdName = runNameS,posClass = pC)
  mdlMetricsXV5$Model <- paste0(mdlMetricsXV5$Model,'_top5')
  mdlMetricsXV10 <- getMdlFitXVresults(mdlFit=mdlFitOpt10,mdName = runNameS,posClass = pC)
  mdlMetricsXV10$Model <- paste0(mdlMetricsXV10$Model,'_top10')
  # > report performance on test set
  mdlMetricsTest5 <- getMdlTestResults(mdlFit=mdlFitOpt5,mdName = runNameS,testSet=inDFtstpp,posClass = pC)
  mdlMetricsTest5$Model <- paste0(mdlMetricsTest5$Model,'_top5')
  mdlMetricsTest10 <- getMdlTestResults(mdlFit=mdlFitOpt10,mdName = runNameS,testSet=inDFtstpp,posClass = pC)
  mdlMetricsTest10$Model <- paste0(mdlMetricsTest10$Model,'_top10')
  # > report performance on external test set (if not doing CD vs UC)
  if (runNameMain != "CD_vs_UC") {
    mdlMetricsTestExt5 <- getMdlTestResults(mdlFit=mdlFitOpt5,mdName = runNameS,testSet=inDFststnegpp,posClass = pC,dataSetName = "Test.set.externalneg")
    mdlMetricsTestExt5$Model <- paste0(mdlMetricsTest5$Model,'_top5')
    mdlMetricsTestExt10 <- getMdlTestResults(mdlFit=mdlFitOpt10,mdName = runNameS,testSet=inDFststnegpp,posClass = pC,dataSetName = "Test.set.externalneg")
    mdlMetricsTestExt10$Model <- paste0(mdlMetricsTest10$Model,'_top10')
  }
  #    - merge tables
  mdlMetricsMrg <- rbind.data.frame(mdlMetricsXV5,mdlMetricsXV10,mdlMetricsTest5,mdlMetricsTest10)
  if (runNameMain != "CD_vs_UC") {
    mdlMetricsMrg <- rbind.data.frame(mdlMetricsMrg,mdlMetricsTestExt5,mdlMetricsTestExt10)
  }
  row.names(mdlMetricsMrg) <- NULL
  #    - save
  write.table(mdlMetricsMrg,paste0('Model_RFE_',runNameMain,'/',runNameS,'_resultsTest_metrics.csv'),sep=',',row.names = F)
  #  > ROC (XV)
  mdlROCxv <- compareModelsTrainingCV(fittedMdls = list(mdlFitOpt5,mdlFitOpt10),
                                      modelNames = c("M1(5 abp)","M2(10 abp)"),roc.conf.boot = 100,
                                      posClass = pC,annotateAUConly = T,roc.conf = 0.95,
                                      tit = paste0(runNameS, " model"),
                                      annotS = 5, diagonalLine = F, textOffSetY = +0.05,textOffSetX = +0.005)
  #  - style it
  mdlRocTestXVS <- mdlROCxv + scale_color_manual(values=c(myCol3,myCol4)) + theme_classic() +
    scale_fill_manual(values = c(myCol3,myCol4) ) + theme(legend.position = "bottom") + ggtitle("") + 
    geom_abline(intercept = c(1), slope = 1,col="darkblue",linetype="longdash",size=1.05) + 
    coord_cartesian(xlim=c(1.005,-0.00),ylim=c(0,1.02),expand = F) +
    theme(axis.line = element_line(size = 1.05),axis.ticks = element_line(size = 0.9))
  
  #print(mdlRocTestXVS) # debug
  ggsave(plot = mdlRocTestXVS,filename = paste0('Model_RFE_',runNameMain,'/',runNameS,'_ROC_XV.png'),dpi = 600,width = 6,height = 6,scale = 1.0)
  
  #  > ROC (test sets)
  mdlRocTest <- compareMdlsDatasets(mdls = list(mdlFitOpt5,mdlFitOpt10),
                                    dataSets = list(inDFtstpp),
                                    posClass = pC,
                                    mdNames = c("M1(5 abp)","M2(10 abp)"),
                                    response = "Cohort",
                                    removeLegend = T,
                                    roc.conf.boot = 100,roc.conf = 0.95,roc.smooth = F,
                                    tit = paste0(runNameS, " model"),
                                    annotateROC = T,
                                    diagonalLine = F,
                                    annotateROCbaseName = T,
                                    annotS = 5,
                                    textOffSetY = +0.05,textOffSetX = +0.005)[[1]]
  mdlRocTestS <- mdlRocTest + scale_color_manual(values=c(myCol3,myCol4)) + theme_classic() +
    scale_fill_manual(values = c(myCol3,myCol4) ) + theme(legend.position = "bottom") + ggtitle("") + 
    geom_abline(intercept = c(1), slope = 1,col="darkblue",linetype="longdash",size=1.05) + 
    coord_cartesian(xlim=c(1.005,-0.00),ylim=c(0,1.02),expand = F) +
    theme(axis.line = element_line(size = 1.05),axis.ticks = element_line(size = 0.9))
  
  print(mdlRocTestS)
  ggsave(plot = mdlRocTestS,filename = paste0('Model_RFE_',runNameMain,'/',runNameS,'_ROC_TestSet.png'),dpi = 600,width = 6,height = 6,scale = 1.0)
  
  # DeLong tests
  # - compare ROC of optimized models to original model, compare top-5 and top-10 optimizations
  resDeLong <- NULL
  # prep ROC
  mdlFitOrig <- readRDS(file = paste0('Model_data_',runNameMain,'/',runNameMain,'_base_all_glmnet_ModelFit.RDS'))
  roc5  <- roc(inDFtstpp$Cohort,predict(mdlFitOpt5,newdata = inDFtstpp,type="prob")[[pC]],auc=T,percent=F)
  roc10 <- roc(inDFtstpp$Cohort,predict(mdlFitOpt10,newdata = inDFtstpp,type="prob")[[pC]],auc=T,percent=F)
  rocOrig <- roc(inDFtstpp$Cohort,predict(mdlFitOrig,newdata = inDFtstpp,type="prob")[[pC]],auc=T,percent=F)
  
  # compare ROCs
  #  > orig vs top5
  t <- roc.test(rocOrig,roc5,method="delong")
  resDeLong <- rbind.data.frame(resDeLong,
                                data.frame(dataset=runNameMain,
                                           method="glmnet",
                                           ROC1=paste0("no_opt"),
                                           ROC1_AUC=t$roc1$auc,
                                           ROC2=paste0("top5"),
                                           ROC2_AUC=t$roc2$auc,
                                           Zstat=t$statistic,
                                           pvalue=t$p.value))
  #  > orig vs top10
  t <- roc.test(rocOrig,roc10,method="delong")
  resDeLong <- rbind.data.frame(resDeLong,
                                data.frame(dataset=runNameMain,
                                           method="glmnet",
                                           ROC1=paste0("no_opt"),
                                           ROC1_AUC=t$roc1$auc,
                                           ROC2=paste0("top10"),
                                           ROC2_AUC=t$roc2$auc,
                                           Zstat=t$statistic,
                                           pvalue=t$p.value))
  #  > top5 vs top10
  t <- roc.test(roc5,roc10,method="delong")
  resDeLong <- rbind.data.frame(resDeLong,
                                data.frame(dataset=runNameMain,
                                           method="glmnet",
                                           ROC1=paste0("top5"),
                                           ROC1_AUC=t$roc1$auc,
                                           ROC2=paste0("top10"),
                                           ROC2_AUC=t$roc2$auc,
                                           Zstat=t$statistic,
                                           pvalue=t$p.value))
  
  # #  > ROC curve with all 3 models (test sets)
  mdlRocTest <- compareMdlsDatasets(mdls = list(mdlFitOpt5,mdlFitOpt10,mdlFitOrig),
                                    dataSets = list(inDFtstpp),
                                    posClass = pC,
                                    mdNames = c("M1(5 abp)","M2(10 abp)","M3(all)"),
                                    response = "Cohort",
                                    removeLegend = T,
                                    roc.conf.boot = 100,roc.conf = 0.95,roc.smooth = F,
                                    tit = paste0(runNameS, " model"),
                                    annotateROC = T,
                                    diagonalLine = F,
                                    annotateROCbaseName = T,
                                    annotS = 5,
                                    textOffSetY = +0.05,textOffSetX = +0.005)[[1]]
  #  - style it
  mdlRocTestS <- mdlRocTest + scale_color_manual(values=c(myCol3,myCol4,"orange")) + theme_classic() +
    scale_fill_manual(values = c(myCol3,myCol4,"orange") ) + theme(legend.position = "bottom") + ggtitle("") +
    geom_abline(intercept = c(1), slope = 1,col="darkblue",linetype="longdash",size=1.05) +
    coord_cartesian(xlim=c(1.005,-0.00),ylim=c(0,1.02),expand = F) +
    theme(axis.line = element_line(size = 1.05),axis.ticks = element_line(size = 0.9))
  
  #print(mdlRocTestS) # debug
  
  write.table(resDeLong,paste0('Results_merged/','test_DeLong_',runNameMain,'.csv'),sep=',',row.names = F)
}
# merge delong tests and do FDR correction
inDL <- read.table('Results_merged/test_DeLong_CD.csv',sep=',',header=T)
inDL <- rbind.data.frame(inDL,read.table('Results_merged/test_DeLong_UC.csv',sep=',',header=T))
inDL <- rbind.data.frame(inDL,read.table('Results_merged/test_DeLong_CD_vs_UC.csv',sep=',',header=T))
inDL$FDR <- p.adjust(inDL$pvalue)

write.table(inDL,paste0('Results_merged/tests_DeLong_merged_FDRcorr.csv'),sep=',',row.names = F)

