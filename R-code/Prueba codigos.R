library(aquap2)
library(openxlsx)
library(mda)

source("C:/Users/Dell/Desktop/R_analysis/Codes/LDAsum.R")


#fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = "custom@importTRH_VoltcraftFizAut_Lew.R", slType = "xls", multiplyRows = F, dol = F)
fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = F, slType = "xls", multiplyRows = F, dol = F)
fullData$header$C_all<-as.factor(fullData$header$C_all)



petreatments <- list(c("sgol@2-13-0"), c("sgol@2-17-0"), c("sgol@2-21-0"), c("snv"), c("msc"), c("deTr"), c("deTr", "snv"), c("deTr", "msc"),
                     c("sgol@2-13-0", "snv"), c("sgol@2-17-0",  "snv"), c("sgol@2-21-0","snv"), c("sgol@2-13-0", "msc"), c("sgol@2-17-0",  "msc"), c("sgol@2-21-0","msc"),
                     c("sgol@2-13-0", "deTr"), c("sgol@2-17-0",  "deTr"), c("sgol@2-21-0","deTr"),
                     c("sgol@2-13-0", "deTr","snv"), c("sgol@2-17-0", "deTr", "snv"), c("sgol@2-21-0","deTr", "snv"),
                     c("sgol@2-13-0", "deTr", "msc"), c("sgol@2-17-0",  "deTr", "msc"), c("sgol@2-21-0","deTr", "msc"),
                     
                     c("sgol@2-21-0", "sgol@2-21-1"), c("sgol@2-21-0", "sgol@2-21-2"),
                     c("sgol@2-21-0", "sgol@2-13-1"), c("sgol@2-21-0", "sgol@2-13-2"),
                     c("sgol@2-21-0", "sgol@2-17-1"), c("sgol@2-21-0", "sgol@2-17-2"),
                     
                     c("sgol@2-13-0", "sgol@2-21-1"),c("sgol@2-13-0", "sgol@2-21-2"),
                     c("sgol@2-13-0", "sgol@2-13-1"),c("sgol@2-13-0", "sgol@2-13-2"),
                     c("sgol@2-13-0", "sgol@2-17-1"),c("sgol@2-13-0", "sgol@2-17-2"),
                     
                     c("sgol@2-17-0", "sgol@2-21-1"),c("sgol@2-17-0", "sgol@2-21-2"),
                     c("sgol@2-17-0", "sgol@2-17-1"),  c("sgol@2-17-0", "sgol@2-17-2"),
                     c("sgol@2-17-0", "sgol@2-13-1"),  c("sgol@2-17-0", "sgol@2-13-2"))



fullDataLDA<-fullData
fullDataLDA$header$C_typecod
fullDataLDA<-ssc(fullDataLDA, C_typecod %in% c("N","P","A"), include = T)


petreatments <- list(c("sgol@2-13-0"), c("sgol@2-17-0"), c("sgol@2-21-0"), c("snv"), c("msc"), c("deTr"), c("deTr", "snv"), c("deTr", "msc"))

i=1

#LDA group CV

for(i in 1:length(petreatments)){
  
  
  pret <- petreatments[[i]]
  
  cu <- gdmm(fullDataLDA, getap(dpt.pre = pret))
  
  #plot_pca(cu, pg.fns = paste0(petreatments[[i]], collapse = "_"))
  data  <- getcd(cu, 1)
  
  plot_spectra(data, colorBy = c("C_conc", "C_tempLev","C_typecode"), pg.fns = paste0(petreatments[[i]], collapse = "_", "_fullData") )
  
  
  source("R-code/FDAnirDeveloper.R")
  
  #source("D:/egyetem/PhD/honey/Rproject_elemzesek/kodok/FDAnirDeveloper.R")
  fullDatadout2 <- data
  
  fullData <- 1
  fullData$header <- fullDatadout2$header
  fullData$colRep <- fullDatadout2$colRep
  fullData$NIR <- fullDatadout2$NIR
  
  GR <- "C_conc"
  a <- daOptima(fullData,GR,N=10)
  dn<- which(a[,2]==max(a[,2]))[1]
  
  valid <- "C_samplename" 
  
  nrG <- nrow(unique(fullData$header[GR])) #how many groups do you have in your classification model
  
  nvG <- nrow(unique(fullData$header[valid])) #how many groups do you have within your validation group
  
  conf1 <- confval <- array(NA, c(nrG, nrG, nvG))
  
  
  colnames(conf1) <- levels(fullData$header[[GR]])
  rownames(conf1) <- levels(fullData$header[[GR]])
  colnames(confval) <- levels(fullData$header[[GR]])
  rownames(confval) <- levels(fullData$header[[GR]])
  
  
  sampleColor<- fullData$colRep[[GR]]
  pdf(paste0("results/Prueba", paste0(petreatments[[i]], collapse = "_"),  "_fullData_groupval.pdf"))
# dev.off()
   tranind <- levels(fullData$header[[valid]])
  
  #dn=50
  j=1
  for (j in 1:nvG) {
    
    ind  <-tranind[j] ##chooses which group to leave out forvalidation
    tri <- -which(fullData$header[[valid]] == ind) #gets the rows that belongs to de ind (tri=training index)
    da1 <- FDA_PCAmodels(fullData, ClassVar = GR, NrPCs = 1:dn,  TranInd = tri)
    plotFDAmodel2(da1, groups= fullData$header[[GR]][tri], projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], cex = 1.2)
    
    #legend("topleft", legend= unique(dadata$P_Type), col = unique(dadata$col$P_Type), pch = 16, cex = 1)
    conf1[,,j] <- da1$confusion
    confval[,,j] <- da1$confusionVal
    #plotFDAmodel(da1, projCVres = TRUE, col = dadata$col$Sample)
    #legend("bottomright", legend= unique(dadata$Sample), col = unique(dadata$col$Sample), pch = 16, cex = 0.8
    conf1_ave <- apply(conf1, MARGIN = c(1,2), mean)
    #(conf1_ave2 <- round(conf1_ave/colSums(conf1_ave)*100, 1))
    confval_ave <- apply(confval, c(1,2), mean)
    #(confval_ave <- round(confval_ave/colSums(confval_ave)*100, 1))
    
    
    
    for(k in 1:nrG){ 
      conf1_ave[, k] <- round(conf1_ave[, k]/sum(conf1_ave[, k])*100, 2)
      
      confval_ave[, k] <- round(confval_ave[, k]/sum(confval_ave[, k])*100, 2)
      
    }
    
    
    
    accuracy_c <- sum(diag(conf1_ave))/sum(conf1_ave)*100
    
    accuracy_cv <-sum(diag(confval_ave))/sum(confval_ave)*100
  }
  
  
  
  conf <- list(ave_c =conf1_ave, ave_cv = confval_ave, recognition = accuracy_c, prediction= accuracy_cv)
  
  write.xlsx(conf, file = paste0("results/Prueba", paste0(petreatments[[i]], collapse = "_"),  "_fullData.xlsx"))
  dev.off()
}




impLDAres("results/Sunflower/all/GCV", name = "Sunflower_all_GCV")


#LDA 3cv



for(i in 1:length(petreatments)){
  
  
  
  pret <- petreatments[[i]]
  
  cu <- gdmm(Sunflower2, getap(dpt.pre = pret))
  
  #plot_pca(cu, pg.fns = paste0(petreatments[[i]], collapse = "_"))
  data  <- getcd(cu, 1)
  
  plot_spectra(data, colorBy = c("C_level", "C_Syrup", "C_aqua"), pg.fns = paste0(petreatments[[i]], collapse = "_", "_Sunflower_all") )
  
  
  source("D:/egyetem/PhD/honey/Rproject_elemzesek/kodok/FDAnirDeveloper.R")
  fullDatadout2 <- data
  
  fullData <- 1
  fullData$header <- fullDatadout2$header
  fullData$colRep <- fullDatadout2$colRep
  fullData$NIR <- fullDatadout2$NIR
  
  GR <- "C_level"
  a <- daOptima(fullData,GR,N=40)
  dn<- which(a[,2]==max(a[,2]))[1]
  
  
  conf1 <- confval <- array(NA, c(10, 10, 3))
  
  colnames(conf1) <- levels(fullData$header$C_level)
  rownames(conf1) <- levels(fullData$header$C_level)
  colnames(confval) <- levels(fullData$header$C_level)
  rownames(confval) <- levels(fullData$header$C_level)
  
  
  sampleColor<- fullData$colRep$C_level
  pdf(paste0("results/Sunflower/all/3CV/NIR_lda_", paste0(petreatments[[i]], collapse = "_"),  "_Sunflower_all_3cv.pdf"))
  tranind <- levels(fullData$header$C_Group)
  
  for (j in 1:3) {
    
    
    tri <- -seq(j, nrow(fullData$header), 3)
    da1 <- FDA_PCAmodels(fullData, ClassVar = "C_level", NrPCs = 1:dn,  TranInd = tri)
    plotFDAmodel2(da1, groups= fullData$header$C_level[tri], projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], cex = 1.2)
    
    #legend("topleft", legend= unique(dadata$P_Type), col = unique(dadata$col$P_Type), pch = 16, cex = 1)
    conf1[,,j] <- da1$confusion
    confval[,,j] <- da1$confusionVal
    #plotFDAmodel(da1, projCVres = TRUE, col = dadata$col$Sample)
    #legend("bottomright", legend= unique(dadata$Sample), col = unique(dadata$col$Sample), pch = 16, cex = 0.8)
    conf1_ave <- apply(conf1, MARGIN = c(1,2), mean)
    #(conf1_ave2 <- round(conf1_ave/colSums(conf1_ave)*100, 1))
    confval_ave <- apply(confval, c(1,2), mean)
    #(confval_ave <- round(confval_ave/colSums(confval_ave)*100, 1))
    
    
    
    for(k in 1:10){ 
      conf1_ave[, k] <- round(conf1_ave[, k]/sum(conf1_ave[, k])*100, 2)
      
      confval_ave[, k] <- round(confval_ave[, k]/sum(confval_ave[, k])*100, 2)
      
    }
    
    
    
    accuracy_c <- sum(diag(conf1_ave))/sum(conf1_ave)*100
    
    accuracy_cv <-sum(diag(confval_ave))/sum(confval_ave)*100
  }
  
  
  
  conf <- list(ave_c =conf1_ave, ave_cv = confval_ave, recognition = accuracy_c, prediction= accuracy_cv)
  
  write.xlsx(conf, file = paste0("results/Sunflower/all/3CV/NIR_lda_", paste0(petreatments[[i]], collapse = "_"),  "_Sunflower_all_CV3.xlsx"))
  dev.off()
}


impLDAres("results/Sunflower/all/3CV", name = "Sunflower_all_3CV")




plssum <- data.frame(matrix(NA, nrow = length(petreatments)*2*2, ncol = 9))
colnames(plssum) <- c("variable", "validtype", "pretreat", "NrLV", "NrObs","R2Tr", "RMSEC", "R2CV", "RMSECV")
rownames(plssum) <- 1:(length(petreatments)*2*2)
source("D:/egyetem/PhD/honey/Rproject_elemzesek/plot_plsr.R")  #source the code I sent!!!

l <- 0
for(k  in c("Y_ttemp", "Y_ttime")){
  
  for(j in c("C_sample", "C_SampleID")){
    
    for(i in 1:length(petreatments)){
      
      #load("R-data/fulldataout2")
      
      
      pret <- petreatments[[i]]
      
      valid <- j
      
      var <- k  #rewrite this
      
      cu <- gdmm(fullData2, getap(dpt.pre = pret , pls.colorBy = "C_ttemp", do.pls = TRUE, pls.ncomp = NULL, pls.regOn = var, pls.valid =  valid, pls.exOut = TRUE))
      
      #save(cu, file = paste0("exports/cu_",var,"_", paste0(pret,collapse = "_")))
      
      plot(cu, pg.fns = paste0("_", k, j, "_",  paste0(pret,collapse = "_"), "_PLS_all"))
      
      
      plsModel <- cu@.Data[[1]]@plsr$model[[1]]
      R2tr <-  getR2C(plsModel)
      RMSE <-  getRMSEC(plsModel)
      R2CV <-  getR2CV(plsModel)
      RMSECV <-  getRMSECV(plsModel)
      LV <- plsModel$ncomp
      Nr <- nrow(plsModel$scores)
      
      
      
      plssum[i+l,1]<- k
      plssum[i+l,2]<- j
      plssum[i+l,3]<-  paste0(pret,collapse = "_")
      plssum[i+l,4]<- LV
      plssum[i+l,5]<- Nr
      plssum[i+l,6]<- R2tr
      plssum[i+l,7]<- RMSE
      plssum[i+l,8]<- R2CV
      plssum[i+l,9]<- RMSECV
      print(plssum) 
      
      write.xlsx(plssum, "results/plssum_Y_temp_all.xlsx")  #rewrite this
      
      
    } 
    l <- l + i
  }
}


impLDAres("C:/Users/Dell/Desktop/R_analysis/20210921_Probiotics_metri/results/Treatmets_excel", name = "prueba")



#for 3 fold cross valid
#valid<- 3
#and change in: conf1 <- confval <- array(NA, c(GR, GR, 3))
