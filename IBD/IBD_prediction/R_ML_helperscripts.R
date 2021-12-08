# ===========================================================================
# by R.Gacesa (UMCG)
# HELPER FUNCTIONS FOR CARET
# ===========================================================================

styleROC <- function(rocGG) {
  rocGG <- rocGG + theme_classic() + ggtitle('') + 
    geom_abline(intercept = c(1), slope = 1,col="darkblue",size=1.05)
  rocGG
}


# prepprocess data for antibody models
# PRE-PROCESS TRAINING SET
# - note: we do not do any pre-processing on training set
#         to avoid data leakage
# - returns list(preprocessing scheme, pre-processed data)
# =====================================
# > remove features with prevalence < X% & > Y%
prepProcessAbData <- function(inDFtr=inDFtr,pVF=0.05,maxP=0.99,minP=0.01,verbose=T) {
  if (verbose) {print(' > pre-processing run')}
  vPrev <- colSums(inDFtr[,colnames(inDFtr) != c("Cohort")])/nrow(inDFtr)
  toRemP <- vPrev[vPrev > maxP | vPrev < minP]
  if (verbose) {print(paste0('  >> dropping ',length(toRemP),' vars with prevalence > ',maxP,' or < ',minP))}
  inDFtr <- inDFtr[,!(colnames(inDFtr) %in% names(toRemP) )]
  if (verbose) {print(paste0('   >> keeping ',ncol(inDFtr)-1,' features'))}
  if (ncol(inDFtr) <= 1) {
    quit('ERROR: ALL FEATURES DROPPED!!')
  }
  # > statistical filter preprocess
  # - do simple chi-squared for differences between classes
  #   (basically keep everything that shows some difference)
  resT <- NULL
  for (cc in colnames(inDFtr)) {
    if (cc != "Cohort") {
      tst <- fisher.test(table(inDFtr[[cc]],inDFtr$Cohort))
      resT <- rbind.data.frame(resT,data.frame(Var=cc,pV=tst$p.value))
    }
  }
  resT <- resT[order(resT$pV),]
  res <- resT$Var[resT$pV < pVF]
  if (verbose) {print(paste0(' >> keeping ',length(res),' predictors after chi-sq filter [pV < ',pVF,']'))}
  toKeep <- c("Cohort",res)
  inDFtrp <- inDFtr[,colnames(inDFtr) %in% toKeep]
  # > "automagical" Caret preprocess
  #  - remove zero-variance variables (zv)
  #  - remove highly correlated variables ( corr > cutoff)
  prepP <- preProcess(inDFtrp,
                      c("zv","corr","conditionalX"),
                      cutoff = 0.99 # correlation coefficient
  )
  if (verbose) {print(prepP)}
  inDFtrpp <- predict(prepP,inDFtrp)
  list(prepP,inDFtrpp)
}

# prep data for antibody models
# ==================================
prepAbData <- function(inDF=inDF,subSet = "all") {
  # subset features as needed
  if (subSet == "agilent") {
    # agilent only
    inDF <- inDF[,c(1,grep('agilent_',colnames(inDF)))]
  } else if (subSet == "twist") {
    # twist only
    inDF <- inDF[,c(1,grep('twist_',colnames(inDF)))]
  } else if (subSet != "all") {
    quit("ERROR: subset must be one of <'agilent','twist','all'>")
  }
  # > convert response var to factor for 2-way prediction
  cn <- c("Cohort")
  #cn <- colnames(inDF)
  for (cc in cn) {
    inDF[[cc]] <- as.character(inDF[[cc]])
    inDF[[cc]][inDF[[cc]]=="0"] <- "N"
    inDF[[cc]][inDF[[cc]]=="1"] <- "Y"
    inDF[[cc]] <- as.factor(inDF[[cc]])
  }
  inDF
}

# get model results on crossvalidation
getMdlFitXVresults <- function(mdlFit,mdName="test",posClass="Y",negClass="N") {
  fittedMdl <- extractBestTunePred(mdlFit)
  print(paste0(' >> analyzing x-validation performance'))
  confM <- confusionMatrix(fittedMdl$pred,fittedMdl$obs)
  # acc, kappa
  acc <- confM[[3]][1]
  kappa <- confM[[3]][2]
  # sens spec ppv npv prec recall F1 prev bacc
  sens <- confM[[4]][1]
  spec <- confM[[4]][2]
  ppv <- confM[[4]][3]
  npv <- confM[[4]][4]
  f1 <- confM[[4]][7]
  bacc <- confM[[4]][11]
  # get auc
  r <- roc(predictor = fittedMdl[[posClass]], response = fittedMdl$obs, auc=T,percent=T,smooth=F)
  auc <- as.numeric(r$auc)/100
  nrVars <- length(mdlFit$coefnames)
  # get nr vars which are non0
  varsCoef <- getVarImpTbl(mdlFit)
  if (!is.null(varsCoef)) {
    nrVarsNZ <- sum(varsCoef$RelImp != 0)
  } else {
    nrVarsNZ <- nrVars
  }
  # turn confusion matrix to wide form for one-row output
  cmTab <- as.data.frame(confM$table)
  p1r1name <- paste0('CM.Pr.',cmTab$Prediction[1],'_Ref.',cmTab$Reference[1]) 
  p1r1nmbr <- cmTab$Freq[1]
  p2r2name <- paste0('CM.Pr.',cmTab$Prediction[2],'_Ref.',cmTab$Reference[2]) 
  p2r2nmbr <- cmTab$Freq[2]
  p3r3name <- paste0('CM.Pr.',cmTab$Prediction[3],'_Ref.',cmTab$Reference[3]) 
  p3r3nmbr <- cmTab$Freq[3]
  p4r4name <- paste0('CM.Pr.',cmTab$Prediction[4],'_Ref.',cmTab$Reference[4]) 
  p4r4nmbr <- cmTab$Freq[4]
  cmNR = p1r1nmbr+p2r2nmbr+p3r3nmbr+p4r4nmbr
  oneRow <- data.frame(Model=mdName,Method=mdlFit$method,Dataset="TrSet.XV",NrVars=nrVars,NrVarsNonZero=nrVarsNZ,ACC=acc,Kappa=kappa,Sensitivity=sens,Specificity=spec,
                       PPV=ppv,NPV=npv,B.ACC=bacc,F1=f1,AUC=auc,
                       p1r1name=p1r1nmbr,p2r2name=p2r2nmbr,p3r3name=p3r3nmbr,p4r4name=p4r4nmbr,CM.NR=cmNR)
  colnames(oneRow)[15:18] <- c(p1r1name,p2r2name,p3r3name,p4r4name)
  oneRow
}

# report metrics for model on test set
# ====================================
getMdlTestResults <- function(mdlFit,testSet,mdName="test",responseVar="Cohort",posClass="Y",negClass="N",
                              dataSetName="Test.set") {
  print(paste0(' >> analyzing test-set performance'))
  # #   - predict on test set
  predTst <- predict(mdlFit,testSet)
  confM <- confusionMatrix(predTst,testSet[[responseVar]])
  # acc, kappa
  acc <- confM[[3]][1]
  kappa <- confM[[3]][2]
  # sens spec ppv npv prec recall F1 prev bacc
  sens <- confM[[4]][1]
  spec <- confM[[4]][2]
  ppv <- confM[[4]][3]
  npv <- confM[[4]][4]
  f1 <- confM[[4]][7]
  bacc <- confM[[4]][11]
  # get auc
  if (sum(testSet[[responseVar]]==posClass) == 0 | sum(testSet[[responseVar]]!=posClass) == 0) {
    auc <- NA
  } else {
    r <- roc(predictor = predict(mdlFit,testSet,type = "prob")[[posClass]],response = testSet[[responseVar]],auc=T,percent=T)
    auc <- as.numeric(r$auc)/100
  }
  # NR Vars
  nrVars <- length(mdlFit$coefnames)
  #   nr vars which are non0
  varsCoef <- getVarImpTbl(mdlFit)
  if (!is.null(varsCoef)) {
    nrVarsNZ <- sum(varsCoef$RelImp != 0)
  } else {
    nrVarsNZ <- nrVars
  }
  # turn confusion matrix to wide form for one-row output
  cmTab <- as.data.frame(confM$table)
  p1r1name <- paste0('CM.Pr.',cmTab$Prediction[1],'_Ref.',cmTab$Reference[1]) 
  p1r1nmbr <- cmTab$Freq[1]
  p2r2name <- paste0('CM.Pr.',cmTab$Prediction[2],'_Ref.',cmTab$Reference[2]) 
  p2r2nmbr <- cmTab$Freq[2]
  p3r3name <- paste0('CM.Pr.',cmTab$Prediction[3],'_Ref.',cmTab$Reference[3]) 
  p3r3nmbr <- cmTab$Freq[3]
  p4r4name <- paste0('CM.Pr.',cmTab$Prediction[4],'_Ref.',cmTab$Reference[4]) 
  p4r4nmbr <- cmTab$Freq[4]
  cmNR = p1r1nmbr+p2r2nmbr+p3r3nmbr+p4r4nmbr
  oneRow <- data.frame(Model=mdName,Method=mdlFit$method,Dataset=dataSetName,NrVars=nrVars,NrVarsNonZero=nrVarsNZ,ACC=acc,Kappa=kappa,Sensitivity=sens,Specificity=spec,
                       PPV=ppv,NPV=npv,B.ACC=bacc,F1=f1,AUC=auc,
                       p1r1name=p1r1nmbr,p2r2name=p2r2nmbr,p3r3name=p3r3nmbr,p4r4name=p4r4nmbr,CM.NR=cmNR)
  colnames(oneRow)[15:18] <- c(p1r1name,p2r2name,p3r3name,p4r4name)
  oneRow
}

# get variable importance (in pretty table)
getVarImpTbl <- function(mdlFit) {
  if (mdlFit$method == "glmnet") {
    betas <- getGlmnetBetas(mdlFit)
    betas$RelImpScaled <- abs(betas$Beta/betas$Beta[1])
    viTbl <- betas
  } else if (mdlFit$method == "gbm") {
    viTbl <- as.data.frame(summary.gbm(mdlFit$finalModel))
    colnames(viTbl) <- c("Var","Importance")
    viTbl$RelImpScaled <- abs(viTbl$Importance/viTbl$Importance[1])
    rownames(viTbl) <- NULL
  } else if (mdlFit$method == "avNNet") {
    viTbl <- NULL
    #viTbl <- as.data.frame(varImp(mdlFit$finalModel))
    #viTbl$Var <- mdlFit$finalModel$coefnames
    #viTbl <- viTbl[order(viTbl$Overall,decreasing = T),]
  } else if (mdlFit$method == "svmRadial") {
    viTbl <- NULL
  }
  else {
    viTbl <- as.data.frame(varImp(mdlFit,scale = F)$importance)
    viTbl$Var <- rownames(viTbl)
    viTbl <- viTbl[order(viTbl$Overall,decreasing = T),]
    viTbl <- viTbl[,c("Var","Overall")]
    colnames(viTbl) <- c("Var","Importance")
    viTbl$RelImp <- viTbl$Importance/sum(viTbl$Importance)*100
    viTbl$RelImpScaled <- abs(viTbl$Importance/viTbl$Importance[1])
    rownames(viTbl) <- NULL
  }
  viTbl
}



# helper function for extracting best tuning prediction data from Caret model
# ===========================================================================
extractBestTunePred <- function(inCaret) {
  bt <- inCaret$bestTune
  if (bt[1] == "none") {
    res <- inCaret$pred
  } else {
    res <- inCaret$pred
    for (v in colnames(bt)) {
      #print(v)
      res <- res[res[[v]] == bt[[v]][1], ]
    }
  }
  res
}

plotROC <- function (rocIn,ggTitle="ROC",textSize=16,AUCsize=5) {
  g <- pROC::ggroc(rocIn,col="darkblue",size=1.5,alpha=0.75)
  g <- g + 
    #annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    annotate("text", x = .35, y = .13, label = paste("AUC =",round(rocIn$auc,2)),hjust = 0,size=AUCsize) +
#    annotate("text", x = .45, y = .30, label = paste("Test Set AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    geom_abline(slope = 1,intercept = 1,col="grey") +
    ggtitle(ggTitle) + xlab("Specificity") + ylab("Sensitivity") +
    labs(color='Dataset')  + theme(text=element_text(size=textSize))
  g
}

draw_confusion_matrix <- function(cmtrx,colP="blue",colN="orange") {
  total <- sum(cmtrx$table)
  res <- as.numeric(cmtrx$table)
  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- brewer.pal(9,'Greens');
  redPalette <- brewer.pal(9,'Reds');
  bluePalette <- brewer.pal(9,'Blues');
  orangePalette <- brewer.pal(9,'Oranges');
  
  getColor <- function (col = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    if (col == "green") {
      palette <- greenPalette
    } else if (col == "red") {
      palette <- redPalette
    } else if (col == "blue") {
      palette <- bluePalette
    } else if (col == "orange") {
      palette <- orangePalette
    }
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }
  
  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  # create the matrix
  classes = colnames(cmtrx$table)
  rect(150, 430, 240, 370, col=getColor(colP, res[1]))
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor(colN, res[3]))
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Reference', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor(colN, res[2]))
  rect(250, 305, 340, 365, col=getColor(colP, res[4]))
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)
  # add in the cmtrx results
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  # add in the specifics
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cmtrx$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cmtrx$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cmtrx$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cmtrx$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cmtrx$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cmtrx$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cmtrx$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cmtrx$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cmtrx$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cmtrx$byClass[7]), 3), cex=1.2)
  # add in the accuracy information
  text(30, 35, names(cmtrx$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cmtrx$overall[1]), 3), cex=1.4)
  text(70, 35, names(cmtrx$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cmtrx$overall[2]), 3), cex=1.4)
}

getGlmnetBetas <- function(mdlFit) {
   mdlBetas <- as.data.frame(as.matrix(coef(mdlFit$finalModel, mdlFit$bestTune$lambda)))
   mdlBetas$Var <- row.names(mdlBetas); colnames(mdlBetas)[1] <- "Beta"
   mdlBetas <- mdlBetas[mdlBetas$Var != "(Intercept)",]
   mdlBetas <- mdlBetas[order(abs(mdlBetas$Beta),decreasing = T),]
   mdlBetas$RelImp <- abs(mdlBetas$Beta)/sum(abs(mdlBetas$Beta))*100
   mdlBetas <- mdlBetas[,c("Var","Beta","RelImp")]; row.names(mdlBetas) <- NULL
   mdlBetas
}
