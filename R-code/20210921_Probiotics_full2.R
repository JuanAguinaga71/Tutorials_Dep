
#Neccesary libraries and supporting codes
library(aquap2) 
library(openxlsx)
#ibrary(mda)
#source("R-code/FDAnirDeveloper.R")
source("C:/Users/Dell/Desktop/R_analysis/Codes/LDAsum.R")
#source("R-code/forpls.r")  



###Import the Data

#fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = "custom@importTRH_VoltcraftFizAut_Lew.R", slType = "xls", multiplyRows = F, dol = F)
fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = F, slType = "xls", multiplyRows = F, dol = F)
fullData$header$C_all<-as.factor(fullData$header$C_all) #converts character data to factor

#fullData_average<-do_avg(fullData,c("C_samplename", "Y_avgS"))

###Select Data of interest
fullData_nowater<- ssc(fullData, C_typecod %in% c("N","P","A"), include = T)


###Plot  SPECTRA
fullData1<-selectWls(fullData_nowater,from=950, to=1630)##select wavelenght of interest
fullData1<- reColor(fullData1)
cb<- c("C_conc", "C_tempLev","C_typecode","C_samplename")
plot_spectra(fullData1, colorBy=  cb, pg.fns="_rawspectra")
#Note: Minimal apply sgolay pretreatment

save (fullData1, file= "R-data/fullData1")
load (fullData1, file= "R-data/fullData1")

###PCA and PLS

cu1 <- gdmm(fullData1, getap(do.pca= T, pca.colorBy =cb, do.pls= T, pls.valid= "C_tempLev",
                            pls.regOn = c("Y_conc", "Y_tempLev")))
plot_pca(cu1, pg.fns= "PCA__fullData")
plot_pls(cu1, pg.fns= "PLS__validbytempLev_fullData")

#Outliers detection
#Options 1) Manual detection in PCA 2) Manual detection in DA

#1) Manual detection in PCA
cu <- gdmm(fullData1, getap(do.pca= T, pca.colorBy =cb, do.pls= F, pls.valid="C_smplID", 
                                    pls.regOn =c("Y_conc", "Y_tempLev")))
sm<-getcm(cu, 1)
plot(sm$model$scores)
identify(sm$model$scores)
#datos_for_outs<-sm$model$scores
outs<-c(361, 362,363)
fullData2<-fullData1[-outs]


#################################################################################
###################LDA para 950-1630 #####################
#Option 1: LDA with sensor contributing arrows (only can be performed 1 pretreatment at the time) uses DA optima
#Option 2: LDA with DA optima (can perform several pretreatment at the time)

fullDataLDA<-fullData2
fullDataLDA<-ssc(fullDataLDA, C_typecod %in% c("N","P","A"), include = T)
fullDataLDA$header$C_ID <- seq(1, nrow(fullDataLDA), 1)
fullDataLDA <- reColor(fullDataLDA)
fullDataLDA$header$C_typecod

fullData_forLDA1<-fullData2
fullData_forLDA2<-fullData2
#################################################################################

####LDA con with sensor contributing arrows(Option 1)####

#For this DA is used: 
library(mda)
source("R-code/FDAnirDeveloper.R")

daData<- 1
daData$header<- fullData_forLDA1$header
daData$colRep<- fullData_forLDA1$colRep
daData$NIR<- fullData_forLDA1$NIR

#GR by "C_typecode" (probiotic), "C_tempLev", "C_conc"
GR <- "C_conc"
a <- daOptima(fullData_forLDA1,GR,N=30)
dn<- which(a[,2]==max(a[,2]))[1]

Title0<- paste0("DA of","C_conc")
Title<- paste0(Title0, "\n", "Correction:", "\n", "sample num: ", 
               nrow(fullData_forLDA1),", LV:", dn,"(100)")
#valid <- "C_Repl"

daObjectList<- list()

pdf("results/test.pdf")
for (i in 1:3) {
  TranInd1 = -seq(i+1, nrow(daData$header), 3)
  daObjectList[[i]] <- FDA_PCAmodels(daData,GR, NrPCs = 1:dn,  TranInd = TranInd1)
  sampleColor = daData$colRep[[GR]]
  plotFDAmodel(daObjectList[[i]], projCVres = TRUE, col = sampleColor[TranInd1], col2 = sampleColor[-TranInd1], main = Title)
  legend("bottomright", legend= unique(daData$heade[[GR]]), col = unique(sampleColor), pch = 16, cex = 0.8)
  plotFDAmodel2(daObjectList[[i]], groups= daData$header[[GR]][TranInd1], projCVres = TRUE, col = sampleColor[TranInd1], col2 = sampleColor[-TranInd1], main = Title, sub=sub)
  d <- SensorContributionTongue(fullData = daData, daObjectList[[i]], corner = "")
}
dev.off()

####################################################################
#For this DA is used: 
library(mda)
source("R-code/FDAnirDeveloper.R")


P_forLDA2<-ssc(fullData_forLDA2, C_typecod %in% "P", include = T)
P_T25_forLDA2<- ssc(P_forLDA2, C_tempLev %in% "T1", include = T)

#####
######

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

petreatments <- list(c("sgol@2-13-0"), c("sgol@2-17-0"))

for(i in 1:length(petreatments)){
  
  
  pret <- petreatments[[i]]
  
  cu <- gdmm(P_T25_forLDA2, getap(dpt.pre = pret))
  #plot_pca(cu, pg.fns = paste0(petreatments[[i]], collapse = "_"))
  
  data  <- getcd(cu, 1)  
  
  plot_spectra(data, colorBy = c("C_conc", "C_tempLev","C_typecode"), 
               pg.fns = paste0(petreatments[[i]], collapse = "_", "_PT25") )
             
  fullDatadout2 <- data
  
  fullData <- 1
  fullData$header <- fullDatadout2$header
  fullData$colRep <- fullDatadout2$colRep
  fullData$NIR <- fullDatadout2$NIR
  
  #GR by "C_typecode" (probiotic), "C_tempLev", "C_conc"
  GR <- "C_conc"
  a <- daOptima(fullData,GR,N=30)
  dn<- which(a[,2]==max(a[,2]))[1]
  
  valid <- "C_Repl" 
  
  nrG <- nrow(unique(fullData$header[GR])) #how many groups do you have in your classification model
  
  nvG <- nrow(unique(fullData$header[valid])) #how many groups do you have within your validation group
  
  conf1 <- confval <- array(NA, c(nrG, nrG, nvG))
  
  
  colnames(conf1) <- levels(fullData$header[[GR]])
  rownames(conf1) <- levels(fullData$header[[GR]])
  colnames(confval) <- levels(fullData$header[[GR]])
  rownames(confval) <- levels(fullData$header[[GR]])
  
  
  sampleColor<- fullData$colRep[[GR]]
  
  pdf(paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "_PT25.pdf"))

  tranind <- levels(fullData$header[[valid]])

  #j=1
  for (j in 1:nvG) {
    
    ind  <-tranind[j] ##chooses which group to leave out forvalidation
    tri <- -which(fullData$header[[valid]] == ind) #gets the rows that belongs to de ind (tri=training index)
    da1 <- FDA_PCAmodels(fullData, ClassVar = GR, NrPCs = 1:dn,  TranInd = tri)
    plotFDAmodel2(da1, groups= fullData$header[[GR]][tri], projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], cex = 1.2)
    
    
  #  d <- SensorContributionTongue(fullData = daData, da1, corner = "") ####esta sirve para crear contribucion de wavelengths
    
    
    #Ir al apartado de LDA outlayers y correr las lineas:
    #plotFDAmodel(da1, projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], pch = 16 , pch2 = 1)
    #identify(da1$scoorsVal)
    
    conf1[,,j] <- da1$confusion
    confval[,,j] <- da1$confusionVal
    conf1_ave <- apply(conf1, MARGIN = c(1,2), mean)
    confval_ave <- apply(confval, c(1,2), mean)
    
    
    
    for(k in 1:nrG){ 
      conf1_ave[, k] <- round(conf1_ave[, k]/sum(conf1_ave[, k])*100, 2)
      
      confval_ave[, k] <- round(confval_ave[, k]/sum(confval_ave[, k])*100, 2)
      
    }
    
    
    
    accuracy_c <- sum(diag(conf1_ave))/sum(conf1_ave)*100
    
    accuracy_cv <-sum(diag(confval_ave))/sum(confval_ave)*100
  }
  
  
  
  conf <- list(ave_c =conf1_ave, ave_cv = confval_ave, recognition = accuracy_c, prediction= accuracy_cv)
  
  
  #write.xlsx(conf, file = paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "PT25.xlsx"))
  write.xlsx(conf, file = paste0("results/LDA_results/", paste0(petreatments[[i]], collapse = "_"), "-" , "PT25.xlsx"))
  dev.off()
}
source("C:/Users/Dell/Desktop/R_analysis/Codes/LDAsum.R")
#impLDAres("results/LDA_950_1630/PT25", name = "_resume_LDA_PT25")
impLDAres("results/LDA_results", name = "_resume_LDA_PT25")

#########PLSR_MICROBIOLOGÍA#################################

#Option 0: PLSR normal
#Option 1: PLSR_independentData (can perform one pretreatment at the time)
#Option 2: PLSR (can perform several pretreatment at the time)


fullDataLDA<-fullData2
fullDataLDA<-ssc(fullDataLDA, C_typecod %in% c("N","P","A"), include = T)
fullDataLDA$header$C_ID <- seq(1, nrow(fullDataLDA), 1)
fullDataLDA <- reColor(fullDataLDA)
fullDataLDA$header$C_typecod

fullData_forPLS0<-fullData2
fullData_forPLS1<-fullData2
fullData_forPLS2<-fullData2
fullData_forPLS0<-fullData_forLDA1

#Option 0: PLSR normal

fullData_forPLS0<-do_avg(fullData_forPLS0,c("C_samplename", "Y_Repl"))
fullData_forPLS0<- reColor(fullData_forPLS0)
cu_pls1 <- gdmm(fullData_forPLS1, getap(do.pca= F, pca.colorBy =cb, do.pls= T, pls.valid= "C_microb",
                            pls.regOn = c("Y_conc", "Y_tempLev", "Y_microb")))
plot_pls(cu_pls1, pg.fns= "PLS__validbymicrob_fullData")

#Option 1: PLSR_independentData (can perform one pretreatment at the time)
#indepPLSR #Select your best LDA pretreatment 

fullData_forPLS2<-fullData_forLDA1
fullData_forPLS2_avg<-do_avg(fullData_forPLS2,c("C_samplename", "Y_Repl"))
#fullData_forPLS1<-do_avg(fullData_forPLS1,c( "C_samplename", "Y_samplename"), inParallel = FALSE) #new code from Flora
fullData_forPLS2_avg<- reColor(fullData_forPLS2_avg)

cu <- gdmm(fullData_forPLS2_avg, getap(dpt.pre = c("sgol@2-13-0", "sgol@2-13-2")))
data  <- getcd(cu, 1)

#plot_pls(cu_pls1, pg.fns= "PLS__validbymicrob_fullData")

fullData_cal<- ssc(data, C_Repl %in% c("R3"), include=FALSE)
fullData_test<- ssc(data, C_Repl %in% c("R3"), include=TRUE)

cu<-gdmm(fullData_cal, getap(do.pca = F, do.pls=T, pca.colorBy =cb, pls.regOn="Y_microb_log", pls.ncomp= 9, pls.valid="C_samplename"))
plot_pls_indepPred(fullData_test, cu, pg.fns= "fullData_PLSind_by_Prob")



#Option 2: PLSR (can perform several pretreatment at the time)
#validation by "C_Repl" 

source("R-code/forpls.r")


fullData_forPLS1r<-fullData_forPLS1
fullData_forPLS1r<- reColor(fullData_forPLS1r)

regon <- list(c("Y_microb_log"))

#validate <- list(c(3), "LOO", "C_Group");list("C_samplename")
validate <- list("C_samplename")
petreatments <- list(c("sgol@2-13-0", "sgol@2-13-2"))

plssum <- data.frame(matrix(NA, nrow = length(petreatments)*length(regon)*length(validate), ncol = 9))
colnames(plssum) <- c("variable", "validtype", "pretreat", "NrLV", "NrObs","R2Tr", "RMSEC", "R2CV", "RMSECV")
rownames(plssum) <- 1:(length(petreatments)*length(regon)*length(validate))

l <- 0
for(k  in 1:length(regon)){
  
  for(j in 1:length(validate)){
    
    for(i in 1:length(petreatments)){
      
      pret <- petreatments[[i]]
      
      valid <- validate[[j]]
      
      var <- regon[[k]]  
      
      cu <- gdmm(fullData_forPLS1, getap(dpt.pre = pret , do.pca=F, pls.colorBy = "C_samplename", do.pls = TRUE, pls.ncomp = 9, pls.regOn = var, pls.valid =  valid, pls.exOut = FALSE))
      
      plot(cu, pg.fns = paste0("_","_", regon[[k]], "_",valid, "_",  paste0(pret,collapse = "_"), "_PLS_full"))
      
      plsModel <- cu@.Data[[1]]@plsr$model[[1]]
      R2tr <-  getR2C(plsModel)
      RMSE <-  getRMSEC(plsModel)
      R2CV <-  getR2CV(plsModel)
      RMSECV <-  getRMSECV(plsModel)
      LV <- plsModel$ncomp
      Nr <- nrow(plsModel$scores)
      
      
      
      plssum[i+l,1]<- var
      plssum[i+l,2]<- valid
      plssum[i+l,3]<-  paste0(pret,collapse = "_")
      plssum[i+l,4]<- LV
      plssum[i+l,5]<- Nr
      plssum[i+l,6]<- R2tr
      plssum[i+l,7]<- RMSE
      plssum[i+l,8]<- R2CV
      plssum[i+l,9]<- RMSECV
      print(plssum) 
      
      write.xlsx(plssum, "results/plssum_all.xlsx")  #rewrite this
      
      
    } 
    l <- l + i
  }
}

pg.fns = paste0("", subs, "",valid, "", "TEMP_detailed", paste0(pret,collapse = ""), "_PLS_all")


