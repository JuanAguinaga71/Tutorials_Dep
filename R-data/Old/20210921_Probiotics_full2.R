library(aquap2)
library(openxlsx)
library(mda)

source("R-code/FDAnirDeveloper.R")
source("C:/Users/Dell/Desktop/R_analysis/Codes/LDAsum.R")
source("R-code/forpls.r")  #source the code I sent!!!


#fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = "custom@importTRH_VoltcraftFizAut_Lew.R", slType = "xls", multiplyRows = F, dol = F)
fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = F, slType = "xls", multiplyRows = F, dol = F)
fullData$header$C_all<-as.factor(fullData$header$C_all)
cb<- c("C_conc", "C_tempLev","C_typecode","C_samplename")
levels(fullData$header$C_typecod)
fullData$header$C_microb_log
fullData_average<-do_avg(fullData,c("C_samplename", "Y_avgS"))

fullData_nowater<- ssc(fullData, C_typecod %in% c("N","P","A"), include = T)


########################### SPECTRA, pca, pls #########################################
########################Sin aplicar pretreatments#####################################################
fullData<-selectWls(fullData,from=950, to=1630)
fullData<- reColor(fullData)
#fullData_forLDA<-fullData

plot_spectra(fullData_average, colorBy=  cb, pg.fns="_rawspectra")
fullDatatrunc <- selectWls(fullData,from=1300, to=1600)
fullDatatrunc<- reColor(fullDatatrunc)

plot_spectra(fullDatatrunc, colorBy=  cb, pg.fns="_rawspectra_1300-1600")

save (fullData, file= "R-data/fullData")
load (fullData, file= "R-data/fullData")
save (fullDatatrunc, file= "R-data/fullDatatrunc")
load (fullDatatrunc, file= "R-data/fullDatatrunc")

cu1 <- gdmm(fullData, getap(do.pca= T, pca.colorBy =cb, do.pls= T, pls.valid= "C_tempLev",
                            pls.regOn = c("Y_conc", "Y_tempLev")))
plot_pca(cu1, pg.fns= "PCA__fullData_950-1630")
plot_pls(cu1, pg.fns= "PLS__validbytempLev_fullData_950-1630")

cu2 <- gdmm(fullDatatrunc, getap(do.pca= T, pca.colorBy =cb, do.pls= T, pls.valid= "C_tempLev",
                                 pls.regOn = c("Y_conc", "Y_tempLev")))
plot_pca(cu2, pg.fns= "PCA__probiotic_1300-1600")
plot_pls(cu2, pg.fns= "PLS__validbytempLev_fullData_1300-1600")


#Desde aqui se trabajará primero con Wavelenght 950-1630

#Primera detección de outlayers de acuerdo a PCA
cu_fullData <- gdmm(fullData, getap(do.pca= T, pca.colorBy =cb, do.pls= F, pls.valid="C_smplID", 
                                    pls.regOn =c("Y_conc", "Y_tempLev")))
dev.off()
sm<-getcm(cu_fullData, 1)
plot(sm$model$scores)
identify(sm$model$scores)
datos_for_outs<-sm$model$scores

outs_prueba<-c(37, 38,39, 415, 416, 417, 547, 548, 549)
fullData_outs_prueba$header$C_all<-fullData[outs]
fullData_outs_prueba$header$C_all
fullData$header[outs_prueba,] #para ver a que muestras corresponden los outlayers
#[c(1:3),]

#################################################################################
###################LDA para 950-1630 #####################
#################################################################################
load (fullData, file= "R-data/fullData")
load (fullDatatrunc, file= "R-data/fullDatatrunc")

fullDataLDA<-fullData
fullDataLDA<-ssc(fullDataLDA, C_typecod %in% c("N","P","A"), include = T)


fullDataLDA$header$C_ID <- seq(1, nrow(fullDataLDA), 1)
fullDataLDA <- reColor(fullDataLDA)

fullDataLDA$header$C_typecod
#################################################################################
################################################################################
#OUTLIER DETECTION for LDA models - solo correr una vez (para fullData)
#1) Run the model and visually detect if the outliers are in the training or validation set 
#2) Create an identifier for the data to be analized. e.g. C_ID 
fullDataLDA$header$C_ID <- seq(1, nrow(fullDataLDA), 1)
fullDataLDA <- reColor(fullDataLDA)
#3)Identify the outliers(clicking)
#correr j=(1, 2, 3, etc) hasta plotFDA, hasta que se cargue el da1 correcto 
plotFDAmodel(da1, projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], pch = 16 , pch2 = 1) #se requiere de esta manera porque no reconoce con elipsis (plotFDAmodel2) 
identify(da1$scoorsVal)

##########j=1
#4)Identify the outlayers in the data (depending if it is training or validation)
fullData$header[-tri,][c(55,  56,  57, 121, 122, 123, 136, 137, 138, 169,
                         170, 171, 196, 197, 198, 202, 203, 204),] ##to identify the 
#5)kick out the outlayers
indexout <- which(fullDataLDA$header$C_ID %in% c(178:180, 349:351, 409:411, 514:516,
                                                 595:597, 622:624))
fullDataLDA2 <- fullDataLDA[-indexout,]

###########j=2
#4)Identify the outlayers in the data (depending if it is training or validation)
fullData$header[-tri,][c(151, 152, 153, 184, 185, 186),] ##to identify the 
#5)kick out the outlayers

#j=2
indexout2 <- which(fullDataLDA2$header$C_ID %in% c(433:435, 541:543))
fullDataLDA3 <- fullDataLDA2[-indexout2,]       

###########j=3
#4)Identify the outlayers in the data (depending if it is training or validation)
fullData$header[-tri,][c(163, 164, 165),] ##to identify the 
#5)kick out the outlayers

#j=3
indexout3 <- which(fullDataLDA3$header$C_ID %in% c(487:489))
fullDataLDA4 <- fullDataLDA3[-indexout3,]   


#ultima correccion, se reviso y se encontro 3 outlayers en j=2 (211 212 213).
#Además se elimina los outlayers del PCA inicial: fullData$header[outs_prueba,]

#4)Identify the outlayers in the data (depending if it is training or validation)
fullData$header[-tri,][c(211, 212, 213),] ##to identify the 
#5)kick out the outlayers

#j=3
indexout3 <- which(fullDataLDA3$header$C_ID %in% c(487:489, 649:651)) #j3        j2 
indexout4 <- which(fullDataLDA3$header$C_chron %in% c(37:39, 415:417, 547:549)) # PCA inicial

fullDataLDA4 <- fullDataLDA3[-indexout3,]
fullDataLDA4 <- fullDataLDA4[-indexout4,]

###############ultima-ultima correccion
#j=2(otra vez)
#4)Identify the outlayers in the data (depending if it is training or validation)
fullData$header[-tri,][c(97:99, 124:126),] ##to identify the 
#5)kick out the outlayers

indexout5 <- which(fullDataLDA4$header$C_ID %in% c(298:300, 388:390))
fullDataLDA5 <- fullDataLDA4[-indexout5,]

#############################################################################
#############################################################################
#este paso es solo para partir con la data original y sacar los outlayers 

fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = F, slType = "xls", multiplyRows = F, dol = F)
fullData$header$C_all<-as.factor(fullData$header$C_all)
cb<- c("C_conc", "C_tempLev","C_typecode")
fullData<-selectWls(fullData,from=950, to=1630)
fullData<- reColor(fullData)


fullData_forLDA<-fullData
fullData_forLDA1<-ssc(fullData_forLDA, C_typecod %in% c("N","P","A"), include = T)

fullData_forLDA1$header$C_ID <- seq(1, nrow(fullData_forLDA1), 1)

fullData_forLDA1<-fullData_forLDA1[-indexout,]
fullData_forLDA1<-fullData_forLDA1[-indexout2,]
fullData_forLDA1<-fullData_forLDA1[-indexout3,]
fullData_forLDA1<-fullData_forLDA1[-indexout4,]
fullData_forLDA1<-fullData_forLDA1[-indexout5,]
fullData_forLDA1 <- reColor(fullData_forLDA1) #con este me quedo


############################################################################
###########################################################################

petreatments <- list(c("sgol@2-13-0","deTr", "snv")) #pretreatment de prueba

###############################################################################

#LDA group CV para detectar outlayers

i=1                #se refiere al primer petreatment
for(i in 1:length(petreatments)){
  
  
  pret <- petreatments[[i]]
  
  cu <- gdmm(fullDataLDA5, getap(dpt.pre = pret))
    #plot_pca(cu, pg.fns = paste0(petreatments[[i]], collapse = "_"))
  
  data  <- getcd(cu, 1)  
  
  plot_spectra(data, colorBy = c("C_conc", "C_tempLev","C_typecode"), 
               pg.fns = paste0(petreatments[[i]], collapse = "_", "_fullData") )
 
  fullDatadout2 <- data
  
  fullData <- 1
  fullData$header <- fullDatadout2$header
  fullData$colRep <- fullDatadout2$colRep
  fullData$NIR <- fullDatadout2$NIR
  
  GR <- "C_typecode"
  
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
  pdf(paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "_fullData_950-1630_groupval.pdf"))
   dev.off()
  tranind <- levels(fullData$header[[valid]])
  
  #dn=50
  j=2    #se refiere a la primera validacion
  for (j in 1:nvG) {
    
    ind  <-tranind[j] ##chooses which group to leave out forvalidation
    tri <- -which(fullData$header[[valid]] == ind) #gets the rows that belongs to de ind (tri=training index)
    da1 <- FDA_PCAmodels(fullData, ClassVar = GR, NrPCs = 1:dn,  TranInd = tri)
    plotFDAmodel2(da1, groups= fullData$header[[GR]][tri], projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], cex = 1.2)
    
    #Ir al apartado de LDA outlayers y correr las lineas:
    #plotFDAmodel(da1, projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], pch = 16 , pch2 = 1)
    #identify(da1$scoorsVal)
    
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
  
  write.xlsx(conf, file = paste0("results/LDA_1350_1630", paste0(petreatments[[i]], collapse = "_"),  "_LDAfullData_1350_1630.xlsx"))
  dev.off()
}

impLDAres("results/LDA_1350_1630", name = "_resume_LDAfullData_1350_1630")




###############################################################
################################################################
########################LDA 950-1630 sin outlayers######################

fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = F, slType = "xls", multiplyRows = F, dol = F)
fullData$header$C_all<-as.factor(fullData$header$C_all)
cb<- c("C_conc", "C_tempLev","C_typecode")
fullData<-selectWls(fullData,from=950, to=1630)
fullData<- reColor(fullData)

fullData_forLDA<-fullData
fullData_forLDA1<-ssc(fullData_forLDA, C_typecod %in% c("N","P","A"), include = T)

fullData_forLDA1$header$C_ID <- seq(1, nrow(fullData_forLDA1), 1)

fullData_forLDA1<-fullData_forLDA1[-indexout,]
fullData_forLDA1<-fullData_forLDA1[-indexout2,]
fullData_forLDA1<-fullData_forLDA1[-indexout3,]
fullData_forLDA1<-fullData_forLDA1[-indexout4,]
fullData_forLDA1<-fullData_forLDA1[-indexout5,]
fullData_forLDA1 <- reColor(fullData_forLDA1) #con este me quedo

#save (fullData_forLDA1, file= "R-data/fullData_forLDA1")

load(fullData_forLDA1, file= "R-data/fullData_forLDA1")
#######################################################################

#LDA group CV 

#i=1                #se refiere al primer petreatment
fullData_noouts <- fullData_forLDA1   #
#N_noouts<-ssc(fullData_noouts, C_typecod %in% "N", include = T)
#A_noouts<-ssc(fullData_noouts, C_typecod %in% "A", include = T)
#P_noouts<-ssc(fullData_noouts, C_typecod %in% "P", include = T)
A_noouts$NIR

#save (fullData_noouts, file= "R-data/fullData_noouts")
#load(fullData_noouts, file= "R-data/fullData_noouts")

#petreatments <- list(c("sgol@2-13-0"), c("sgol@2-17-0"), c("sgol@2-21-0"), c("snv"), c("msc"), c("deTr"), c("deTr", "snv"), c("deTr", "msc"),
#                     c("sgol@2-13-0", "snv"), c("sgol@2-17-0",  "snv"), c("sgol@2-21-0","snv"), c("sgol@2-13-0", "msc"), c("sgol@2-17-0",  "msc"), c("sgol@2-21-0","msc"),
#                     c("sgol@2-13-0", "deTr"), c("sgol@2-17-0",  "deTr"), c("sgol@2-21-0","deTr"),
#                     c("sgol@2-13-0", "deTr","snv"), c("sgol@2-17-0", "deTr", "snv"), c("sgol@2-21-0","deTr", "snv"),
#                     c("sgol@2-13-0", "deTr", "msc"), c("sgol@2-17-0",  "deTr", "msc"), c("sgol@2-21-0","deTr", "msc"))

#petreatments <- list(c("sgol@2-13-0"), c("sgol@2-17-0"), c("sgol@2-21-0"))

#####
#este use hasta ahora# petreatments <- list(c("sgol@2-13-0"), c("snv"), c("msc"), c("deTr"), c("deTr", "snv"), 
                     #c("deTr", "msc"), c("sgol@2-13-0", "snv"), c("sgol@2-13-0", "msc"), 
                     #c("sgol@2-13-0", "deTr"),  c("sgol@2-13-0", "deTr","snv"),  
                     #c("sgol@2-13-0", "deTr", "msc"))
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
#petreatments <- list( c("deTr","snv"))
for(i in 1:length(petreatments)){
  
  
  pret <- petreatments[[i]]
  
  cu <- gdmm(fullData_noouts, getap(dpt.pre = pret))
  #plot_pca(cu, pg.fns = paste0(petreatments[[i]], collapse = "_"))
  
  data  <- getcd(cu, 1)  
  
  plot_spectra(data, colorBy = c("C_conc", "C_tempLev","C_typecode"), 
               pg.fns = paste0(petreatments[[i]], collapse = "_", "_fullData_Nouts") )
               #pg.fns = paste0(petreatments[[i]], collapse = "_", "_N_NOuts") )
               #pg.fns = paste0(petreatments[[i]], collapse = "_", "_A_NOuts") )
               #pg.fns = paste0(petreatments[[i]], collapse = "_", "_P_NOuts") )
  
  fullDatadout2 <- data
  
  fullData <- 1
  fullData$header <- fullDatadout2$header
  fullData$colRep <- fullDatadout2$colRep
  fullData$NIR <- fullDatadout2$NIR
  
  #GR by "C_typecode" (probiotic), "C_tempLev", "C_conc"
  GR <- "C_typecode"
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
  pdf(paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "_fullData_Nouts_950-1630_groupval.pdf"))
  #pdf(paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "_N_Nouts_950-1630_groupval.pdf"))
  #pdf(paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "_A_Nouts_950-1630_groupval.pdf"))
  #pdf(paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "_P_Nouts_950-1630_groupval.pdf"))
  
  
  tranind <- levels(fullData$header[[valid]])
  
  #dn=50
  j=3
  for (j in 1:nvG) {
    
    ind  <-tranind[j] ##chooses which group to leave out forvalidation
    tri <- -which(fullData$header[[valid]] == ind) #gets the rows that belongs to de ind (tri=training index)
    da1 <- FDA_PCAmodels(fullData, ClassVar = GR, NrPCs = 1:dn,  TranInd = tri)
    plotFDAmodel2(da1, groups= fullData$header[[GR]][tri], projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], cex = 1.2)
    
    #Ir al apartado de LDA outlayers y correr las lineas:
    #plotFDAmodel(da1, projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], pch = 16 , pch2 = 1)
    #identify(da1$scoorsVal)
    
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
  
  write.xlsx(conf, file = paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "_LDAfullData_NOuts_950_1630.xlsx"))
  #write.xlsx(conf, file = paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "_LDA_N_NOuts_950_1630.xlsx"))
  #write.xlsx(conf, file = paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "_LDA_A_NOuts_950_1630.xlsx"))  
  #write.xlsx(conf, file = paste0("results/LDA_950_1630", paste0(petreatments[[i]], collapse = "_"),  "_LDA_P_NOuts_950_1630.xlsx")) 
    dev.off()
}

impLDAres("results/LDA_950_1630/fullData (Prob_Repl)", name = "_resume_LDAfullData_NOuts_950_1630")
#impLDAres("results/LDA_950_1630/N (conc_Repl)", name = "_resume_LDA_N_NOuts_950_1630")
#impLDAres("results/LDA_950_1630/A (tempLev_Repl)", name = "_resume_LDA_A_NOuts_950_1630")
#impLDAres("results/LDA_950_1630/P (conc_Repl)", name = "_resume_LDA_P_NOuts_950_1630")

###############################################################
################################################################
########################LDA 1300-1600 sin outlayers######################
#LDA group CV 

#i=1                #se refiere al primer petreatment
fullData_noouts <- fullData_forLDA1
fullData_noouts <-selectWls(fullData_noouts,from=1300, to=1600)
fullData_noouts <-reColor(fullData_noouts)
#
N_noouts<-ssc(fullData_noouts, C_typecod %in% "N", include = T)
#A_noouts<-ssc(fullData_noouts, C_typecod %in% "A", include = T)
#P_noouts<-ssc(fullData_noouts, C_typecod %in% "P", include = T)
N_noouts <-reColor(N_noouts)

#save (fullData_noouts, file= "R-data/fullData_noouts")
#load(fullData_noouts, file= "R-data/fullData_noouts")

#petreatments <- list(c("sgol@2-13-0"), c("sgol@2-17-0"), c("sgol@2-21-0"), c("snv"), c("msc"), c("deTr"), c("deTr", "snv"), c("deTr", "msc"),
#                     c("sgol@2-13-0", "snv"), c("sgol@2-17-0",  "snv"), c("sgol@2-21-0","snv"), c("sgol@2-13-0", "msc"), c("sgol@2-17-0",  "msc"), c("sgol@2-21-0","msc"),
#                     c("sgol@2-13-0", "deTr"), c("sgol@2-17-0",  "deTr"), c("sgol@2-21-0","deTr"),
#                     c("sgol@2-13-0", "deTr","snv"), c("sgol@2-17-0", "deTr", "snv"), c("sgol@2-21-0","deTr", "snv"),
#                     c("sgol@2-13-0", "deTr", "msc"), c("sgol@2-17-0",  "deTr", "msc"), c("sgol@2-21-0","deTr", "msc"))

#petreatments <- list(c("sgol@2-13-0"), c("sgol@2-17-0"), c("sgol@2-21-0"))

petreatments <- list(c("sgol@2-13-0"), c("snv"), c("msc"), c("deTr"), c("deTr", "snv"), 
                     c("deTr", "msc"), c("sgol@2-13-0", "snv"), c("sgol@2-13-0", "msc"), 
                     c("sgol@2-13-0", "deTr"),  c("sgol@2-13-0", "deTr","snv"),  
                     c("sgol@2-13-0", "deTr", "msc"))

#petreatments <- list(c("sgol@2-13-0"))

for(i in 1:length(petreatments)){
  
  
  pret <- petreatments[[i]]
  
  cu <- gdmm(N_noouts, getap(dpt.pre = pret))
  #plot_pca(cu, pg.fns = paste0(petreatments[[i]], collapse = "_"))
  
  data  <- getcd(cu, 1)  
  
  plot_spectra(data, colorBy = c("C_conc", "C_tempLev","C_typecode"), 
               #pg.fns = paste0(petreatments[[i]], collapse = "_", "_fullData_Nouts") )
               pg.fns = paste0(petreatments[[i]], collapse = "_", "_N_NOuts") )
               #pg.fns = paste0(petreatments[[i]], collapse = "_", "_A_NOuts") )
               #pg.fns = paste0(petreatments[[i]], collapse = "_", "_P_NOuts") )
  
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
  #pdf(paste0("results/LDA_1300_1600", paste0(petreatments[[i]], collapse = "_"),  "_fullData_Nouts_1300_1600_groupval.pdf"))
  pdf(paste0("results/LDA_1300_1600", paste0(petreatments[[i]], collapse = "_"),  "_N_Nouts_1300_1600_groupval.pdf"))
  #pdf(paste0("results/LDA_1300_1600", paste0(petreatments[[i]], collapse = "_"),  "_A_Nouts_1300_1600_groupval.pdf"))
  #pdf(paste0("results/LDA_1300_1600", paste0(petreatments[[i]], collapse = "_"),  "_P_Nouts_1300_1600_groupval.pdf"))
  
  
  tranind <- levels(fullData$header[[valid]])
  
  #dn=50
  j=3
  for (j in 1:nvG) {
    
    ind  <-tranind[j] ##chooses which group to leave out forvalidation
    tri <- -which(fullData$header[[valid]] == ind) #gets the rows that belongs to de ind (tri=training index)
    da1 <- FDA_PCAmodels(fullData, ClassVar = GR, NrPCs = 1:dn,  TranInd = tri)
    plotFDAmodel2(da1, groups= fullData$header[[GR]][tri], projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], cex = 1.2)
    
    #Ir al apartado de LDA outlayers y correr las lineas:
    #plotFDAmodel(da1, projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], pch = 16 , pch2 = 1)
    #identify(da1$scoorsVal)
    
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
  
  #write.xlsx(conf, file = paste0("results/LDA_1300_1600", paste0(petreatments[[i]], collapse = "_"),  "_LDAfullData_NOuts_1300_1600.xlsx"))
  write.xlsx(conf, file = paste0("results/LDA_1300_1600", paste0(petreatments[[i]], collapse = "_"),  "_LDA_N_NOuts_1300_1600.xlsx"))
  #write.xlsx(conf, file = paste0("results/LDA_1300_1600", paste0(petreatments[[i]], collapse = "_"),  "_LDA_A_NOuts_1300_1600.xlsx"))  
  #write.xlsx(conf, file = paste0("results/LDA_1300_1600", paste0(petreatments[[i]], collapse = "_"),  "_LDA_P_NOuts_1300_1600.xlsx")) 
  dev.off()
}

#impLDAres("results/LDA_1300_1600/fullData (conc_Repl)", name = "_resume_LDAfullData_NOuts_1300_1600")
impLDAres("results/LDA_1300_1600/N (conc_Repl)", name = "_resume_LDA_N_NOuts_1300_1600")
#impLDAres("results/LDA_1300_1600/A (conc_Repl)", name = "_resume_LDA_A_NOuts_1300_1600")
#impLDAres("results/LDA_1300_1600/P (tempLev_Repl)", name = "_resume_LDA_P_NOuts_1300_1600")

##################Aquaphotomics 950-1630###############################################
##############################################################################

#i=1                #se refiere al primer petreatment
fullData_noouts <- fullData_forLDA1
#fullData_noouts <-selectWls(fullData_noouts,from=1300, to=1600)
fullData_noouts <-reColor(fullData_noouts)
#
#N_noouts<-ssc(fullData_noouts, C_typecod %in% "N", include = T)
#A_noouts<-ssc(fullData_noouts, C_typecod %in% "A", include = T)
#P_noouts<-ssc(fullData_noouts, C_typecod %in% "P", include = T)
#N_noouts <-reColor(N_noouts)




findClosest <- function(what, where){ # often problem the non integer wls, this fun can find the closest to what (num vect) in where (num vect) & gives back the indexes
  inexes <- NA
  for (i in what){
    inexes <- c(inexes, order(abs(where-i))[1])
  }
  out <- inexes[-1]
} # this is a function to find the closest wavelengths

aqWls <- .ap2$stn$aqg_wlsAquagram # 12 predefined water matrics coordinates
wls <- getWavelengths(fullData_noouts) # wavelength available in the dataset
ind <- findClosest(aqWls, wls) # indexes of the closest wavelengths

selWls <- wls[ind]


cu_Aq1 <- gdmm(fullData_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
                         aqg.vars= c("C_typecode", "C_tempLev", "C_conc"), 
                         aqg.TCalib= "def", do.pca=F, do.pls=F))

plot_aqg(cu_Aq1, pg.fns="_fullData_950-1630")

#### N

N_noouts<-ssc(fullData_noouts, C_typecod %in% "N", include = T)
N_noouts <-reColor(N_noouts)


cu_Aq2 <- gdmm(N_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
                                  aqg.vars= c("C_tempLev","C_conc"),aqg.minus="T1", 
                                  aqg.TCalib= "def", do.pca=F, do.pls=F))

plot_aqg(cu_Aq2, pg.fns="_N_950-1630")

#by Temp
NT1_noouts<-ssc(N_noouts, C_tempLev %in% "T1", include = T)
#NT2_noouts<-ssc(N_noouts, C_tempLev %in% "T2", include = T)
#NT3_noouts<-ssc(N_noouts, C_tempLev %in% "T3", include = T)


#aqg.selWls = c(1342, 1364, 1374, 1384, 1412, 1426, 1440, 1452, 1462, 1476, 1488, 1512),

cu_Aq3 <- gdmm(NT1_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
                                  aqg.vars= c("C_conc"), 
                                  aqg.TCalib= "def", do.pca=F, do.pls=F))

plot_aqg(cu_Aq3, pg.fns="_NT1_950-1630")
#plot_aqg(cu_Aq3, pg.fns="_NT2_950-1630")
#plot_aqg(cu_Aq3, pg.fns="_NT3_950-1630")

#by conc

NC1_noouts<-ssc(N_noouts, C_conc %in% "C1", include = T)
NC2_noouts<-ssc(N_noouts, C_conc %in% "C2", include = T)
NC3_noouts<-ssc(N_noouts, C_conc %in% "C3", include = T)



cu_Aq4 <- gdmm(NC3_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
                                    aqg.vars= c("C_tempLev"), 
                                    aqg.TCalib= "def", do.pca=F, do.pls=F))

plot_aqg(cu_Aq3, pg.fns="_NC1_950-1630")
plot_aqg(cu_Aq3, pg.fns="_NC2_950-1630")
plot_aqg(cu_Aq3, pg.fns="_NC3_950-1630")

##################Aquaphotomics 1300-1600###############################################
##############################################################################

#i=1                #se refiere al primer petreatment
fullData_noouts <- fullData_forLDA1
fullData_noouts <-selectWls(fullData_noouts,from=1300, to=1600)
fullData_noouts <-reColor(fullData_noouts)
#
N_noouts<-ssc(fullData_noouts, C_typecod %in% "N", include = T)
#A_noouts<-ssc(fullData_noouts, C_typecod %in% "A", include = T)
#P_noouts<-ssc(fullData_noouts, C_typecod %in% "P", include = T)
#N_noouts <-reColor(N_noouts)

plot_spectra(N_noouts,colorBy = cb)



findClosest <- function(what, where){ # often problem the non integer wls, this fun can find the closest to what (num vect) in where (num vect) & gives back the indexes
  inexes <- NA
  for (i in what){
    inexes <- c(inexes, order(abs(where-i))[1])
  }
  out <- inexes[-1]
} # this is a function to find the closest wavelengths

aqWls <- .ap2$stn$aqg_wlsAquagram # 12 predefined water matrics coordinates
wls <- getWavelengths(fullData_noouts) # wavelength available in the dataset
ind <- findClosest(aqWls, wls) # indexes of the closest wavelengths

selWls <- wls[ind]


cu_Aq1 <- gdmm(fullData_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
                                         aqg.vars= c("C_typecode", "C_tempLev", "C_conc"), 
                                         aqg.TCalib= "def", do.pca=F, do.pls=F))

plot_aqg(cu_Aq1, pg.fns="_fullData_1300-1600")

#### N

N_noouts<-ssc(fullData_noouts, C_typecod %in% "N", include = T)
N_noouts <-reColor(N_noouts)


#cu_Aq2 <- gdmm(N_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
#                                  aqg.vars= c("C_tempLev"), aqg.minus = "T1" 
#                                  aqg.TCalib= "def", do.pca=F, do.pls=F))
#cu_Aq2 <- gdmm(N_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
#                                  aqg.vars= c("C_conc"), aqg.minus = "C3", 
#                                  aqg.TCalib= "def", do.pca=F, do.pls=F))

cu_Aq2 <- gdmm(N_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
                                  aqg.vars= c("C_tempLev"), 
                                  aqg.TCalib= "def", do.pca=F, do.pls=F))
cu_Aq2 <- gdmm(N_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
                                  aqg.vars= c("C_conc"), 
                                  aqg.TCalib= "def", do.pca=F, do.pls=F))

plot_aqg(cu_Aq2, pg.fns="_NT_1300-1600")
plot_aqg(cu_Aq2, pg.fns="_NC_1300-1600")
#by Temp
#NT1_noouts<-ssc(N_noouts, C_tempLev %in% "T1", include = T)
#NT2_noouts<-ssc(N_noouts, C_tempLev %in% "T2", include = T)
NT3_noouts<-ssc(N_noouts, C_tempLev %in% "T3", include = T)



#cu_Aq3 <- gdmm(NT1_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
#                                    aqg.vars= c("C_conc"), aqg.minus = "C3", 
#                                    aqg.TCalib= "def", do.pca=F, do.pls=F))

cu_Aq3 <- gdmm(NT3_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
                                    aqg.vars= c("C_conc"), 
                                    aqg.TCalib= "def", do.pca=F, do.pls=F))

#plot_spectra(NT1_noouts,colorBy = cb)
NT1_noouts$NIR


#plot_aqg(cu_Aq3, pg.fns="_NT1_1300-1600")
#plot_aqg(cu_Aq3, pg.fns="_NT2_1300-1600")
plot_aqg(cu_Aq3, pg.fns="_NT3_1300-1600")

#by conc

NC1_noouts<-ssc(N_noouts, C_conc %in% "C1", include = T)
NC2_noouts<-ssc(N_noouts, C_conc %in% "C2", include = T)
NC3_noouts<-ssc(N_noouts, C_conc %in% "C3", include = T)



#cu_Aq4 <- gdmm(NC3_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
#                                    aqg.vars= c("C_tempLev"), aqg.minus = "T1", 
#                                    aqg.TCalib= "def", do.pca=F, do.pls=F))

cu_Aq4 <- gdmm(NC3_noouts, ap=getap(do.aqg=T,aqg.selWls = selWls,
                                    aqg.vars= c("C_tempLev"), 
                                    aqg.TCalib= "def", do.pca=F, do.pls=F))

plot_aqg(cu_Aq4, pg.fns="_NC1_1300-1600")
plot_aqg(cu_Aq4, pg.fns="_NC2_1300-1600")
plot_aqg(cu_Aq4, pg.fns="_NC3_1300-1600")





#prueba
#aqg_wlsAquagram = c(1342, 1364, 1374, 1384, 1412, 1426, 1440, 1452, 1462, 1476, 1488, 1512), 	## the wavelengths for the classic aquagram (argument aqg.selWls)

cu_Aq5 <- gdmm(NC3_noouts, ap=getap(do.aqg=T,aqg.selWls = c(1342, 1364, 1374, 1384, 1412, 1426, 1440, 1452, 1462, 1476, 1488, 1512),
                                    aqg.vars= c("C_tempLev"), 
                                    aqg.TCalib= "def", do.pca=F, do.pls=F))

plot_aqg(cu_Aq3, pg.fns="_prueba NC1_1300-1600")


#########PLSR_MICROBIOLOGÍA#################################
#############################################################
load(fullData_forLDA1, file= "R-data/fullData_forLDA1")
fullData_forPLS1<-fullData_forLDA1

fullData_forPLS1<-do_avg(fullData_forPLS1,c("C_samplename", "Y_avgS"))

fullData_forPLS1<-selectWls(fullData_forPLS1,from=950, to=1630)
fullData_forPLS1<- reColor(fullData_forPLS1)

cu_pls1 <- gdmm(fullData_forPLS1, getap(do.pca= F, pca.colorBy =cb, do.pls= T, pls.valid= "C_microb",
                            pls.regOn = c("Y_conc", "Y_tempLev", "Y_microb")))
plot_pls(cu_pls1, pg.fns= "PLS__validbymicrob_fullData_950-1630")


#validation by "C_Repl" 
###
fullData_forPLS1<-selectWls(fullData_forPLS1,from=950, to=1630)
fullData_forPLS1r<- reColor(fullData_forPLS1)
fullData_forPLS1r$header$Y_microb_log

cu_pls1 <- gdmm(fullData_forPLS1r, getap(do.pca= F, pca.colorBy =cb, do.pls= T, pls.valid= "C_Repl", pls.colorBy= "C_samplename",
                                         pls.regOn = c("Y_conc", "Y_tempLev", "Y_microb", "Y_microb_log")))
plot_pls(cu_pls1, pg.fns= "PLSr__validbymicrob_fullData_950-1630")
###N##
fullData_forPLS1_Nr<-ssc(fullData_forPLS1, C_typecod %in% c("N"), include = T)
fullData_forPLS1_Nr<- reColor(fullData_forPLS1_Nr)

###A
fullData_forPLS1_Ar<-ssc(fullData_forPLS1, C_typecod %in% c("A"), include = T)
fullData_forPLS1_Ar<- reColor(fullData_forPLS1_Ar)

###P
fullData_forPLS1_Pr<-ssc(fullData_forPLS1, C_typecod %in% c("P"), include = T)
fullData_forPLS1_Pr<- reColor(fullData_forPLS1_Pr)

###################
##################




regon <- list(c("Y_microb_log"))

#validate <- list(c(3), "LOO", "C_Group")
validate <- list("C_Repl")
petreatments <- list(c("sgol@2-13-0"), c("deTr", "msc"))

petreatments <- list(c("sgol@2-13-0"), c("snv"), c("msc"), c("deTr"), c("deTr", "snv"), 
                     c("deTr", "msc"), c("sgol@2-13-0", "snv"), c("sgol@2-13-0", "msc"), 
                     c("sgol@2-13-0", "deTr"),  c("sgol@2-13-0", "deTr","snv"),  
                     c("sgol@2-13-0", "deTr", "msc"))


plssum <- data.frame(matrix(NA, nrow = length(petreatments)*length(regon)*length(validate), ncol = 9))
colnames(plssum) <- c("variable", "validtype", "pretreat", "NrLV", "NrObs","R2Tr", "RMSEC", "R2CV", "RMSECV")
rownames(plssum) <- 1:(length(petreatments)*length(regon)*length(validate))
#source("C:/Users/Dell/Desktop/R_analysis/20210921_Probiotics_metri/R-code/forpls.r")  #source the code I sent!!!


library(pls)

l <- 0
for(k  in 1:length(regon)){
  
  for(j in 1:length(validate)){
    
    for(i in 1:length(petreatments)){
      
      #load("R-data/fulldataout2")
      
      
      pret <- petreatments[[i]]
      
      valid <- validate[[j]]
      
      var <- regon[[k]]  
      
      cu <- gdmm(fullData_forPLS1, getap(dpt.pre = pret , do.pca=F, pls.colorBy = "C_samplename", do.pls = TRUE, pls.ncomp = NULL, pls.regOn = var, pls.valid =  valid, pls.exOut = FALSE))
      
      #save(cu, file = paste0("exports/cu_",var,"", paste0(pret,collapse = "")))
      
      plot(cu, pg.fns = paste0("_","_", regon[[k]], "_",valid, "_",  paste0(pret,collapse = "_"), "_PLS_all4"))
      
     
      
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

#"PLS_fullData_950-1630" ; "PLS_N_950-1630"; 

############PCAs_950_1630
fullData_noouts <- fullData_forLDA1
fullData_for_PCA1<-fullData_noouts


cu <- gdmm(fullData_for_PCA1, getap(do.pca= T, pca.colorBy =cb, do.pls= F, pls.valid= "C_tempLev",
                            pls.regOn = c("Y_conc", "Y_tempLev")))
plot_pca(cu, pg.fns= "PCA__fullData_950-1630")


fullData_noouts <- fullData_forLDA1
fullData_for_PCA1<-fullData_noouts

fullData_for_PCA1<-ssc(fullData_for_PCA1, C_typecod %in% c("P"), include = T)
cu <- gdmm(fullData_for_PCA1, getap(do.pca= T, pca.colorBy =cb, do.pls= F, pls.valid= "C_tempLev",
                                    pls.regOn = c("Y_conc", "Y_tempLev")))
plot_pca(cu, pg.fns= "PCA__P_950-1630")


############PCAs_950_1630_con pretreatments seleccionados

fullData_noouts <- fullData_forLDA1
fullData_for_PCA1<-fullData_noouts


cu <- gdmm(fullData_for_PCA1, getap(do.pca= T, pca.colorBy =cb, dpt.pre =c("sgol", "msc"), do.pls= F, pls.valid= "C_tempLev",
                                    pls.regOn = c("Y_conc", "Y_tempLev")))
plot_pca(cu, pg.fns= "PCA__fullData_sgol_msc_950-1630")


fullData_noouts <- fullData_forLDA1
fullData_for_PCA1<-fullData_noouts

fullData_for_PCA1<-ssc(fullData_for_PCA1, C_typecod %in% c("N"), include = T)
cu <- gdmm(fullData_for_PCA1, getap(do.pca= T, pca.colorBy =cb, dpt.pre =c("sgol", "msc"), do.pls= F, pls.valid= "C_tempLev",
                                    pls.regOn = c("Y_conc", "Y_tempLev")))
plot_pca(cu, pg.fns= "PCA__N_sgol_msc_950-1630")


