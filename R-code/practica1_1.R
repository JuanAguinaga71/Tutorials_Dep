#PRACTICA 1 JP
###Primero ir a more y click en "set as working directory" (escoge la carpeta que estarás trabajando)
###Abre las siguientes librerias
library(aquap2)
library(openxlsx)
library(mda)
source("R-code/FDAnir.R")

###Llama a los datos (gfd es get full data)
#fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = "custom@importTRH_VoltcraftFizAut_Lew.R", slType = "xls", multiplyRows = F, dol = F)
fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = F, slType = "xls", multiplyRows = F, dol = F)
fullData <- reColor(fullData)
cb<- c("C_samplename","C_type" ,"C_mixture", "C_variety", "C_adultlev", "C_adultNr","C_conSNr","C_Repl")
plot_spectra(fullData, colorBy = cb, pg.fns= "rawSpec1") ###el numero de gráficos esta dado por cb; pg.fns es el nombre que le damos al pdf

fullData$header$C_all<-as.factor(fullData$header$C_all)
cu <- gdmm(fullData, getap(do.pca = T,pca.colorBy = cb,do.pls = F)) ###gdmm es generate datasets and make models;getap es (get analysis procedure) en este caso PCA
plot_pca(cu,pg.fns = "_fullData1", pg.where = "pdf")


#AP (anato powder)   ###ssc es seleccionar observaciones, y aqui la seleccion es segun C_mixture y luego se especifica cuales son las observaciones deseadas "AP","Control_TU", "Control_Anatto","Control_NA","Control_TY".
fullData_AP_Controls <- ssc(fullData,C_mixture %in% c("AP","Control_TU", "Control_Anatto","Control_NA","Control_TY"), include = TRUE)
plot_spectra(fullData_AP_Controls, colorBy = cb, pg.fns= "rawSpecAP_Controls") ###el numero de gráficos esta dado por cb

#TOMATO_AP_RAW 
fullData_tom <- ssc(fullData_AP_Controls, C_type %in% c("tomato","adulterant"), include = TRUE)###ssc es seleccionar observaciones
plot_spectra(fullData_tom, colorBy = cb, pg.fns= "raw_tomSpec")
fullData_tom <- selectWls(fullData_tom, 1300,1650) #selectWls que se desee analizar
plot_spectra(fullData_tom, colorBy = cb, pg.fns= "raw_tom_wls")
fullData_tom_AP <- fullData_tom 
fullData_tom_AP  <- reColor(fullData_tom_AP) ### recolor es importante redefinir colores en un data set mas pequeño 
fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE) 
fullData_tom_AP <-reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants"))) ###aqg es acuagrama
plot_pca(cu,pg.fns = "_tom_AP", pg.where = "pdf")
plot_pls(cu,pg.fns = "RAW_tom_AP", pg.where = "pdf")
plot_aqg(cu,pg.fns = "_tom_AP", pg.where = "pdf")  
#fullData_tom_AP <- fullData_tom 
#cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
#plot_aqg(cu,pg.fns = "_tom_AP_ByVariety", pg.where = "pdf")  

########################OUTLIERS###############
sm <- getcm(cu, 1) ###getcm es get a single model from the cube por ejemplo antes generado por gdmm 
plot(sm$model$scores, col= fullData_tom_AP$colRep$C_adultlev) ### col significa column indexes
y<- identify(sm$model$scores) ###identify es identificar puntos en un grafico disperso (en este caso para identificar los outlayers)
outs <- c(y)
fullData_tom_AP<- fullData_tom_AP[-outs,]
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pca(cu,pg.fns = "_tom_AP_NoOuts", pg.where = "pdf")
plot_pls(cu,pg.fns = "_tom_AP_NoOuts", pg.where = "pdf")
plot_aqg(cu,pg.fns = "_tom_AP_NoOuts", pg.where = "pdf")  
fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE) 
fullData_tom_AP  <- reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "_tom_AP_ByAdlev_NoOuts", pg.where = "pdf")  



###23/11/2020 llegue hasta aqui

#LDA Tom_AP_AdLev
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR   ####y esto????
dadata$col <- fullData_tom_AP$colRep  ####y esto????
col = dadata$col$C_lev_AP  ####y esto????
pdf("results/fulldata_tom_AP_AdLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd) *###FDA es flexible discriminant analysis
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("bottomleft", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("bottomleft", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("bottomleft", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/Tom_AP_AdLev.xlsx", colNames=TRUE)
dev.off()


#LDA fulldata_tom_AP variety
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_variety

pdf("results/fulldata_tom_AP_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_tom_AP_variety.xlsx", colNames=TRUE)
dev.off()
 
#LDA tom_AP_lowestLev-ByVariety
fullData_tom_AP <-ssc(fullData_tom_AP, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_variety

pdf("results/fulldata_tom_AP_LowesLev_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_tom_AP_Lowestlev_variety.xlsx", colNames=TRUE)
dev.off()


 
# ############################################################
# ############################################################
# ############################################################
# ########### MSC PRETREATED #######################
# ############################################################
# ############################################################
# ############################################################

fullData_tom <- ssc(fullData_AP_Controls, C_type %in% c("tomato","adulterant"), include = TRUE)
fullData_tom  <- reColor(fullData_tom)
fullData_tom  <- do_msc(fullData_tom)
plot_spectra(fullData_tom, colorBy = cb, pg.fns= "MSC_tomSpec")
fullData_tom <- selectWls(fullData_tom, 1300,1650)
plot_spectra(fullData_tom, colorBy = cb, pg.fns= "MSC_tom_wls")
fullData_tom_AP <- fullData_tom 
fullData_tom_AP  <- reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pca(cu,pg.fns = "MSC_tom_AP", pg.where = "pdf")
plot_pls(cu,pg.fns = "MSC_tom_AP", pg.where = "pdf")
plot_aqg(cu,pg.fns = "MSC_tom_AP", pg.where = "pdf")  

fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE) 
fullData_tom_AP <-reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "MSC_tom_AP_level8Notincluded", pg.where = "pdf")  

fullData_tom_AP <- fullData_tom 
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars = "C_variety",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "MSC_tom_AP_ByVariety", pg.where = "pdf")  

########################OUTLIERS###############
sm <- getcm(cu, 1)
plot(sm$model$scores, col= fullData_tom_AP$colRep$C_adultlev)
y<- identify(sm$model$scores)
outs <- c(y)
fullData_tom_AP<- fullData_tom_AP[-outs,]
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.mod = "classic",aqg.vars = "C_variety",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pca(cu,pg.fns = "MSC_tom_AP_NoOuts", pg.where = "pdf")
plot_pls(cu,pg.fns = "MSC_tom_AP_NoOuts", pg.where = "pdf")
plot_aqg(cu,pg.fns = "MSC_tom_AP_ByVarietyNoOuts", pg.where = "pdf")  
###############
fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE) 
fullData_tom_AP<- fullData_tom_AP[-outs,]
fullData_tom_AP  <- reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "MSC_tom_AP_ByAdlev_NoOuts_level8Noincluded", pg.where = "pdf") 
###############
fullData_tom_AP <- fullData_tom 
fullData_tom_AP<- fullData_tom_AP[-outs,]
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "MSC_tom_AP_ByAdlev_NoOuts", pg.where = "pdf")  


#LDA Tom_AP_AdLev
fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE)
fullData_tom_AP<- fullData_tom_AP[-outs,]
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_lev_AP
pdf("results/fulldata_MSC_tom_AP_AdLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("topright", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("topright", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("bottomleft", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/MSC_Tom_AP_AdLev.xlsx", colNames=TRUE)
dev.off()


#LDA fulldata_MSC_tom_AP variety
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_variety

pdf("results/fulldata_MSC_tom_AP_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_MSC_tom_AP_variety.xlsx", colNames=TRUE)
dev.off()

#LDA MSC_tom_AP_lowestLev-ByVariety
fullData_tom_AP <-ssc(fullData_tom_AP, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_variety

pdf("results/fulldata_MSC_tom_AP_LowesLev_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_MSC_tom_AP_Lowestlev_variety.xlsx", colNames=TRUE)
dev.off()

# ############################################################
# ############################################################
# ############################################################
# ########### Sgolay pretreated #######################
# ############################################################
# ############################################################
# ############################################################

fullData_tom <- ssc(fullData_AP_Controls, C_type %in% c("tomato","adulterant"), include = TRUE)
fullData_tom  <- reColor(fullData_tom)
fullData_tom  <- do_sgolay(fullData_tom,p = 2, n = 21)
plot_spectra(fullData_tom, colorBy = cb, pg.fns= "SGOL_tomSpec")
fullData_tom <- selectWls(fullData_tom, 1300,1650)
plot_spectra(fullData_tom, colorBy = cb, pg.fns= "SGOL_tom_wls")
fullData_tom_AP <- fullData_tom 
fullData_tom_AP  <- reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pca(cu,pg.fns = "SGOL_tom_AP", pg.where = "pdf")
plot_pls(cu,pg.fns = "SGOL_tom_AP", pg.where = "pdf")
plot_aqg(cu,pg.fns = "SGOL_tom_AP", pg.where = "pdf")  

fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE) 
fullData_tom_AP <-reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SGOL_tom_AP_level8Notincluded", pg.where = "pdf")  

fullData_tom_AP <- fullData_tom 
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars = "C_variety",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SGOL_tom_AP_ByVariety", pg.where = "pdf")  

########################OUTLIERS###############
fullData_tom_AP <- fullData_tom 
sm <- getcm(cu, 1)
plot(sm$model$scores, col= fullData_tom_AP$colRep$C_adultlev)
y<- identify(sm$model$scores)
outs <- c(y)
fullData_tom_AP<- fullData_tom_AP[-outs,]
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.mod = "classic",aqg.vars = "C_variety",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pca(cu,pg.fns = "SGOL_tom_AP_NoOuts", pg.where = "pdf")
plot_pls(cu,pg.fns = "SGOL_tom_AP_NoOuts", pg.where = "pdf")
plot_aqg(cu,pg.fns = "SGOL_tom_AP_ByVarietyNoOuts", pg.where = "pdf")  
###############
fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE) 
fullData_tom_AP<- fullData_tom_AP[-outs,]
fullData_tom_AP  <- reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SGOL_tom_AP_ByAdlev_NoOuts_level8Noincluded", pg.where = "pdf") 
###############
fullData_tom_AP <- fullData_tom 
fullData_tom_AP<- fullData_tom_AP[-outs,]
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SGOL_tom_AP_ByAdlev_NoOuts", pg.where = "pdf")  


#LDA SGOL_Tom_AP_AdLev
fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE)
fullData_tom_AP<- fullData_tom_AP[-outs,]
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_lev_AP

pdf("results/fulldata_SGOL_tom_AP_AdLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("topleft", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("topright", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("bottomleft", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/SGOL_Tom_AP_AdLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata_SGOL_tom_AP variety
fullData_tom_AP <- fullData_tom 
fullData_tom_AP<- fullData_tom_AP[-outs,]
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_variety

pdf("results/fulldata_SGOL_tom_AP_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_SGOL_tom_AP_variety.xlsx", colNames=TRUE)
dev.off()

#LDA SGOL_tom_AP_lowestLev-ByVariety
fullData_tom_AP <-ssc(fullData_tom_AP, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_variety

pdf("results/fulldata_SGOL_tom_AP_LowesLev_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_SGOL_tom_AP_Lowestlev_variety.xlsx", colNames=TRUE)
dev.off()

 
# ############################################################
# ############################################################
# ############################################################
# ########### Sgolay and MSC pretreated #######################
# ############################################################
# ############################################################
# ############################################################


fullData_tom <- ssc(fullData_AP_Controls, C_type %in% c("tomato","adulterant"), include = TRUE)
fullData_tom  <- reColor(fullData_tom)
fullData_tom  <- do_sgolay(fullData_tom,p = 2, n = 21)
fullData_tom  <- do_msc(fullData_tom)
plot_spectra(fullData_tom, colorBy = cb, pg.fns= "SGOL_tomSpec")
fullData_tom <- selectWls(fullData_tom, 1300,1650)
plot_spectra(fullData_tom, colorBy = cb, pg.fns= "SGOL_MSC_tom_wls")
fullData_tom_AP <- fullData_tom 
fullData_tom_AP  <- reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pca(cu,pg.fns = "SGOL_MSC_tom_AP", pg.where = "pdf")
plot_pls(cu,pg.fns = "SGOL_MSC_tom_AP", pg.where = "pdf")
plot_aqg(cu,pg.fns = "SGOL_MSC_tom_AP", pg.where = "pdf")  

fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE) 
fullData_tom_AP <-reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SGOL_MSC_tom_AP_level8Notincluded", pg.where = "pdf")  

fullData_tom_AP <- fullData_tom 
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars = "C_variety",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SGOL_MSC_tom_AP_ByVariety", pg.where = "pdf")  

########################OUTLIERS###############
fullData_tom_AP <- fullData_tom 
sm <- getcm(cu, 1)
plot(sm$model$scores, col= fullData_tom_AP$colRep$C_adultlev)
y<- identify(sm$model$scores)
outs <- c(y)
fullData_tom_AP<- fullData_tom_AP[-outs,]
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.mod = "classic",aqg.vars = "C_variety",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pca(cu,pg.fns = "SGOL_MSC_tom_AP_NoOuts", pg.where = "pdf")
plot_pls(cu,pg.fns = "SGOL_MSC_tom_AP_NoOuts", pg.where = "pdf")
plot_aqg(cu,pg.fns = "SGOL_MSC_tom_AP_ByVarietyNoOuts", pg.where = "pdf")  
###############
fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE) 
fullData_tom_AP<- fullData_tom_AP[-outs,]
fullData_tom_AP  <- reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SGOL_MSC_tom_AP_ByAdlev_NoOuts_level8Noincluded", pg.where = "pdf") 
###############
fullData_tom_AP <- fullData_tom 
fullData_tom_AP<- fullData_tom_AP[-outs,]
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SGOL_MSC_tom_AP_ByAdlev_NoOuts", pg.where = "pdf")  


#LDA MSC_SGOL_Tom_AP_AdLev
fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE)
fullData_tom_AP<- fullData_tom_AP[-outs,]
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_lev_AP

pdf("results/fulldata_SGOL_MSC_tom_AP_AdLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("topleft", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("topright", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("bottomleft", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/SGOL_MSC_Tom_AP_AdLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata_SGOL_MSC_tom_AP variety
fullData_tom_AP <- fullData_tom 
fullData_tom_AP<- fullData_tom_AP[-outs,]
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_variety

pdf("results/fulldata_SGOL_MSC_tom_AP_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_SGOL_MSC_tom_AP_variety.xlsx", colNames=TRUE)
dev.off()

#LDA SGOL_MSC_tom_AP_lowestLev-ByVariety
fullData_tom_AP <-ssc(fullData_tom_AP, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_variety

pdf("results/fulldata_SGOL_MSC_tom_AP_LowesLev_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomleft", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_SGOL_MSC_tom_AP_Lowestlev_variety.xlsx", colNames=TRUE)
dev.off()



# ############################################################
# ############################################################
# ############################################################
# ########### SNV pretreated #######################
# ############################################################
# ############################################################
# ############################################################

fullData_tom <- ssc(fullData_AP_Controls, C_type %in% c("tomato","adulterant"), include = TRUE)
fullData_tom  <- reColor(fullData_tom)
fullData_tom  <- do_snv(fullData_tom)
plot_spectra(fullData_tom, colorBy = cb, pg.fns= "SNV_tomSpec")
fullData_tom <- selectWls(fullData_tom, 1300,1650)
plot_spectra(fullData_tom, colorBy = cb, pg.fns= "SNV_tom_wls")
fullData_tom_AP <- fullData_tom 
fullData_tom_AP  <- reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = T,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pca(cu,pg.fns = "SNV_tom_AP", pg.where = "pdf")
plot_pls(cu,pg.fns = "SNV_tom_AP", pg.where = "pdf")
plot_aqg(cu,pg.fns = "SNV_tom_AP", pg.where = "pdf")  

fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE) 
fullData_tom_AP <-reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = F,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SNV_tom_AP_level8Notincluded", pg.where = "pdf")  

fullData_tom_AP <- fullData_tom 
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = F,aqg.vars = "C_variety",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SNV_tom_AP_ByVariety", pg.where = "pdf")  

########################OUTLIERS###############
fullData_tom_AP <- fullData_tom 
sm <- getcm(cu, 1)
plot(sm$model$scores, col= fullData_tom_AP$colRep$C_adultlev)
y<- identify(sm$model$scores)
outs <- c(y)
fullData_tom_AP<- fullData_tom_AP[-outs,]
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = F,aqg.mod = "classic",aqg.vars = "C_variety",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pca(cu,pg.fns = "SNV_tom_AP_NoOuts", pg.where = "pdf")
plot_pls(cu,pg.fns = "SNV_tom_AP_NoOuts", pg.where = "pdf")
plot_aqg(cu,pg.fns = "SNV_tom_AP_ByVarietyNoOuts", pg.where = "pdf")  
###############
fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE) 
fullData_tom_AP<- fullData_tom_AP[-outs,]
fullData_tom_AP  <- reColor(fullData_tom_AP)
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = F,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SNV__tom_AP_ByAdlev_NoOuts_level8Noincluded", pg.where = "pdf") 
###############
fullData_tom_AP <- fullData_tom 
fullData_tom_AP<- fullData_tom_AP[-outs,]
cu <- gdmm(fullData_tom_AP, getap(do.pca = T,do.pls = F,aqg.vars="C_lev_AP",aqg.mod = "classic",do.aqg = T,pca.colorBy = cb,pls.regOn=c("Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_aqg(cu,pg.fns = "SGOL_MSC_tom_AP_ByAdlev_NoOuts", pg.where = "pdf")  


#LDA SNV_Tom_AP_AdLev
fullData_tom_AP <- ssc(fullData_tom, C_lev_AP%in% c("0","1","2","3","4","5","6","7"), include = TRUE)
fullData_tom_AP<- fullData_tom_AP[-outs,]
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_lev_AP

pdf("results/fulldata_SNV_tom_AP_AdLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("bottomleft", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("topright", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_lev_AP", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of AP_AdLev")
legend("bottomleft", legend= unique(dadata$C_lev_AP), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/SNV_Tom_AP_AdLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata_SNV_tom_AP variety
fullData_tom_AP <- fullData_tom 
fullData_tom_AP<- fullData_tom_AP[-outs,]
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_variety

pdf("results/fulldata_SNV_tom_AP_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomleft", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_SNV_tom_AP_variety.xlsx", colNames=TRUE)
dev.off()

#LDA SNV_tom_AP_lowestLev-ByVariety
fullData_tom_AP <-ssc(fullData_tom_AP, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL
dadata <- fullData_tom_AP$header
dadata$NIR <- fullData_tom_AP$NIR
dadata$col <- fullData_tom_AP$colRep
col = dadata$col$C_variety

pdf("results/fulldata_SNV_tom_AP_LowesLev_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomleft", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:14,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_SNV_tom_AP_Lowestlev_variety.xlsx", colNames=TRUE)
dev.off()


