FDA_PCAmodels <- function(dataset, ClassVar, TranInd = NA, NrPCs = NA){
  matrix.pca <- prcomp(dataset$NIR)
  PCAscores <- predict(matrix.pca)       
  i <- which(colnames(dataset) == ClassVar)
  if (any(is.na(NrPCs))){
	D <- cbind(dataset[,i], data.frame(PCAscores))
  } else {
	D <- cbind(dataset[,i], data.frame(PCAscores[,NrPCs]))
  }
  
  colnames(D)[1]<-'group'
  if (any(is.na(TranInd))){
    FDA_PCAmodel <- fda(group ~.,data = D, na.action = "na.omit")
    p <- predict(FDA_PCAmodel, newdata = D, type='variates')
    confusion <- confusion(FDA_PCAmodel)
    out <- list(FDA_PCAmodel = FDA_PCAmodel, scoors = p, confusion = confusion)
  } else {
    FDA_PCAmodel <- fda(group ~.,data = D[TranInd, ], na.action = "na.omit")
	pTran <- predict(FDA_PCAmodel, D[TranInd, ], type='variates')
    pVal <- predict(FDA_PCAmodel, D[-TranInd, ], type='variates')
    confusionTran <- confusion(FDA_PCAmodel)
    confusionVal <- confusion(FDA_PCAmodel, D[-TranInd, ])
    out <- list(FDA_PCAmodel = FDA_PCAmodel, scoors = pTran, scoorsVal = pVal, confusion = confusionTran, confusionVal = confusionVal, matrix.pca = matrix.pca)
  }
}

calcAveConfMatrix <- function(FDAmodel1, FDAmodel2, FDAmodel3) {
  ave_c <- c1 <- FDAmodel1$confusion
  c2 <- FDAmodel2$confusion
  c3 <- FDAmodel3$confusion
  for (i in 1:ncol(c1)) {
    for (j in 1:nrow(c1)) {
      ave_c[j, i] <- round(mean(c(c1[j, i], c2[j, i], c3[j, i])), 2)
    }
    ave_c[, i] <- round(ave_c[, i]/sum(ave_c[, i])*100, 2)
  }
  ave_cv <- c1 <- FDAmodel1$confusionVal
  c2 <- FDAmodel2$confusionVal
  c3 <- FDAmodel3$confusionVal
  for (i in 1:ncol(c1)) {
    for (j in 1:nrow(c1)) {
      ave_cv[j, i] <- round(mean(c(c1[j, i], c2[j, i], c3[j, i])), 2)
    }
    ave_cv[, i] <- round(ave_cv[, i]/sum(ave_cv[, i])*100, 2)
    
    
  }
  
  accuracy_c <- sum(diag(ave_c))/sum(ave_c)*100
  
  accuracy_cv <-sum(diag(ave_cv))/sum(ave_cv)*100
  
  return(list(ave_c = ave_c, ave_cv = ave_cv, recognition = accuracy_c, prediction = accuracy_cv))
 
  
} #Eof

getFdaVars <- function(FDA_PCAmodel, round = 2){
  #Vars <- round((LDAmodel$svd)^2/sum((LDAmodel$svd)^2)*100, round)
  round(c(FDA_PCAmodel$percent.explained[1], 
          (FDA_PCAmodel$percent.explained[-1] - FDA_PCAmodel$percent.explained[-length(FDA_PCAmodel$percent.explained)])), round)
}

plotFDAmodel <- function(FDA_PCAmodel, Roots = 1:2, col = 1, pch = 16, main = "", sub = "", projCVres = FALSE, col2 = 1, pch2 = "x", ...){
  p <- FDA_PCAmodel$scoors
  vars <- getFdaVars(FDA_PCAmodel$FDA_PCAmodel)
  if (projCVres) {
	ylims <- apply(rbind(FDA_PCAmodel$scoors, FDA_PCAmodel$scoorsVal), 2, range)
	pCv <- FDA_PCAmodel$scoorsVal
  } else {
	ylims <- apply(FDA_PCAmodel$scoors, 2, range)
  }
	if (ncol(p)>1) {
		xlab=paste0('root ', Roots[1],' - ', vars[Roots[1]], "%")
		ylab=paste0('root ', Roots[2],' - ', vars[Roots[2]], "%")
  
		plot(p[,Roots[1]],p[,Roots[2]], xlab=xlab, ylab=ylab, col = col, pch = pch, main = main, sub = sub, xlim = ylims[,Roots[1]], ylim = ylims[,Roots[2]], ...)
		if (projCVres) {
			points(pCv[,Roots[1]], pCv[,Roots[2]], col = col2, pch = pch2, ...)
		}
	} else {
		plot(p, col = col, pch = pch, main = main, sub = sub, ...)
	}
  }

  
aveConfMatrix <- function(colConf){
	TrInd <- which(names(colConf) == "Training")
	VaInd <- which(names(colConf) == "Valid")
	
	AveTr <- 0
	for (i in TrInd){
		AveTr <- AveTr + (colConf[i]$Training)/colSums(colConf[i]$Training)*100
		#AveTr <- AveTr + t(t(colConf[i]$Training)/colSums(t(colConf[i]$Training))*100)
		
	}
	AveTr <- round(AveTr/3, 2)
	
	AveVa <- 0
	for (i in VaInd){
		 AveVa <- AveVa + (colConf[i]$Valid)/colSums(colConf[i]$Valid)*100
		#AveVa <- AveVa + t(t(colConf[i]$Valid)/colSums(t(colConf[i]$Valid))*100)
	}
	AveVa <- round(AveVa/3, 2)
	out <- list(Training = AveTr, Valid = AveVa)
}
