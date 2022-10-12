library(aquap2)
library(openxlsx)
library(mda)
source("R-code/FDAnir.R")



#fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = "custom@importTRH_VoltcraftFizAut_Lew.R", slType = "xls", multiplyRows = F, dol = F)
fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = F, slType = "xls", multiplyRows = F, dol = F)
plot_spectra(fullData, pg.fns="test_paprika")
#############

fullData<-reColor(fullData)
cb<-c("C_type", "C_mixture")
plot_spectra(fullData, colorBy= cb, pg.fns="tes_paprika2")
selectWls(fullData, 1000, 1200)
plot_spectra(fullData, colorBy= cb, pg.fns="tes_paprika3")


###############

fullData <- reColor(fullData)
#write.xlsx("fullData.xlsx", file="prueba")
#write.xlsx("fullData.xlsx", file = "fullData.xlsx", colNames = TRUE, borders = "columns")
cb<- c("C_samplename","C_type" ,"C_mixture", "C_variety", "C_adultlev", "C_adultNr","C_conSNr","C_Repl")

plot_spectra(fullData, colorBy = cb, pg.fns= "rawSpec")

#paprika

fullData_paprika<-fullData
fullData_paprika <- ssc(fullData, C_type %in% "paprika", include = TRUE) 
#fullData_paprika <- ssc(fullData, "paprika" %in% C_type, include = TRUE) 
fullData_paprika <- reColor(fullData_paprika)
fullData_paprika <- selectWls(fullData_paprika, 950,1650)
plot_spectra(fullData_paprika, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paprika_rawSpec")
#fullData_paprika <- selectWls(fullData_paprika, 800,1200)
#plot_spectra(fullData_paprika, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paprika_rawSpectrita")


fullData$header=as.factor(fullData$header)
fullData_paprika$header$C_all=as.factor(fullData_paprika$header$C_all)
cu <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T, 
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado" ,"Y_sample","Y_adulterants")))
plot_pls(cu, pg.fns= "PLS_paprika1")

cu2 <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                   pls.regOn =c("Y_corn","Y_color","Y_avocado" ,"Y_sample","Y_adulterants")))
plot_pls(cu2, pg.fns= "PLS_paprika_MixtureValidation")


#cu3 <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T, pls.valid="C_variety",
#                                    pls.regOn =c("Y_corn","Y_color","Y_avocado" ,"Y_sample","Y_adulterants")))
#plot_pls(cu3, pg.fns= "PLS_paprika_VarietyValidation")

#cu <- gdmm(fullData_paprika, getap(do.pca= T, pca.colorBy=c("C_variety","C_adultlev","C_adultNr"), do.pls= T, 
#                                   pls.regOn =c("Y_corn","Y_color","Y_avocado" ,"Y_sample","Y_adulterants")))

#plot_pls(cu, pg.fns= "PLS_paprika2")
#plot_pca(cu,pg.fns="PCA_paprika2")

#fullData_paprika_SH <- ssc(fullData_paprika, C_variety %in% "LE", include = TRUE) 
#fullData_paprika_SH <- reColor(fullData_paprika_SH)
#fullData_paprika_SH <- selectWls(fullData_paprika_SH, 950,1650)
#plot_spectra(fullData_paprika_SH, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paprika_SH_rawSpec")

#cu4<- gdmm(fullData_paprika_SH, getap(do.pca= T, pca.colorBy=c("C_adultlev","C_adultNr"), do.pls= T, 
#                                   pls.regOn =c("Y_corn","Y_color","Y_avocado" ,"Y_sample","Y_adulterants")))

#plot_pls(cu4, pg.fns= "PLS_paprika_SH2")
#plot_pca(cu4,pg.fns="PCA_paprika_SH2")




#LDA paprika mixture
dadata <- fullData_paprika$header
dadata$NIR <- fullData_paprika$NIR
dadata$col <- fullData_paprika$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprika_mixture.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:95,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures_D1")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)


TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:95,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:95,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprika_mixture.xlsx", colNames=TRUE)
dev.off()

#mixture but with no controls
fullData_paprika_NoControl <- ssc(fullData_paprika, C_mixture %in% 
                                   c("Control_SH","Control_TU", "Control_M1",
                                     "Control_NA","Control_TY","Control_LE"), include = FALSE) 

dadata <- fullData_paprika_NoControl$header
dadata$NIR <- fullData_paprika_NoControl$NIR
dadata$col <- fullData_paprika_NoControl$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprika_mixture_No_control.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures_D1")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprika_mixture_no_control.xlsx", colNames=TRUE)
dev.off()

#LDA paprika variety
dadata <- fullData_paprika$header
dadata$NIR <- fullData_paprika$NIR
dadata$col <- fullData_paprika$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paprika_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_paprika_variety.xlsx", colNames=TRUE)
dev.off()


#Lowest adulteration level (0.5%)

fullData_paprika_lowestLev <- ssc(fullData_paprika, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL

cu <- gdmm(fullData_paprika_lowestLev, getap(do.pca= F, do.pls= T, 
                                     pls.regOn =c("Y_corn","Y_color","Y_avocado")))
plot_pls(cu, pg.fns= "PLS_paprika_LowestLev")

cu2 <- gdmm(fullData_paprika_lowestLev, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                      pls.regOn =c("Y_corn","Y_color","Y_avocado")))
plot_pls(cu2, pg.fns= "PLS_paprika_MixtureValidation_LowestLev")

#LDA lowest adulteration level
dadata <- fullData_paprika_lowestLev$header
dadata$NIR <- fullData_paprika_lowestLev$NIR
dadata$col <- fullData_paprika_lowestLev$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprika_mixture_LowestLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:50,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:50,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:50,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprika_mixture_LowestLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paprika_lowestLev$header
dadata$NIR <- fullData_paprika_lowestLev$NIR
dadata$col <- fullData_paprika_lowestLev$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paprika_fullData_paprika_lowestLev_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topleft", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomleft", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprika_fullData_paprika_lowestLev_variety.xlsx", colNames=TRUE)
dev.off()


############################################################
############################################################
############################################################
########### MSC PRETREATED #######################
############################################################
############################################################
############################################################

fullData_paprika <- ssc(fullData, C_type %in% "paprika", include = TRUE) 
fullData_paprika <- reColor(fullData_paprika)
fullData_paprika <- do_msc(fullData_paprika)
fullData_paprika <- selectWls(fullData_paprika, 950,1650)
plot_spectra(fullData_paprika, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paprika_MSCSpec")

cu <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T,
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_sample","Y_adulterants")))
plot_pls(cu, pg.fns= "PLS_MSCpaprika")

cu2 <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                   pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_sample","Y_adulterants")))
plot_pls(cu2, pg.fns= "PLS_MSCpaprika_MixtureValidation")

#LDA fulldata mixture
dadata <- fullData_paprika$header
dadata$NIR <- fullData_paprika$NIR
dadata$col <- fullData_paprika$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaMSC_mixture.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:90,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:90,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprika_mixture_MSC.xlsx", colNames=TRUE)
dev.off()

#mixture but with no controls
fullData_paprika_NoControl <- ssc(fullData_paprika, C_mixture %in% 
                                   c("Control_SH","Control_TU", "Control_M1",
                                     "Control_NA","Control_TY","Control_LE"), include = FALSE) 

dadata <- fullData_paprika_NoControl$header
dadata$NIR <- fullData_paprika_NoControl$NIR
dadata$col <- fullData_paprika_NoControl$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaMSC_mixture_No_control.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:100,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:100,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:100,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprikaMSC_mixture_no_control.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paprika$header
dadata$NIR <- fullData_paprika$NIR
dadata$col <- fullData_paprika$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paprikaMSC_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topleft", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topleft", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprikaMSC_variety.xlsx", colNames=TRUE)
dev.off()


#Lowest adulteration level (0.5%)

fullData_paprika_lowestLev <- ssc(fullData_paprika, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL

cu <- gdmm(fullData_paprika_lowestLev, getap(do.pca= F, do.pls= T, 
                                            pls.regOn =c("Y_corn","Y_color","Y_avocado")))
plot_pls(cu, pg.fns= "PLS_paprikaMSC_LowestLev")

cu2 <- gdmm(fullData_paprika_lowestLev, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                      pls.regOn =c("Y_corn","Y_color","Y_avocado")))
plot_pls(cu2, pg.fns= "PLS_paprikaMSC_MixtureValidation_LowestLev")

#LDA fulldata mixture
dadata <- fullData_paprika_lowestLev$header
dadata$NIR <- fullData_paprika_lowestLev$NIR
dadata$col <- fullData_paprika_lowestLev$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaMSC_mixture_LowestLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:40,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:40,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:40,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprikaMSC_mixture_LowestLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paprika_lowestLev$header
dadata$NIR <- fullData_paprika_lowestLev$NIR
dadata$col <- fullData_paprika_lowestLev$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paprikaMSC_fullData_paprika_lowestLev_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_paprikaMSC_fullData_paprika_lowestLev_variety.xlsx", colNames=TRUE)
dev.off()




############################################################
############################################################
############################################################
########### Sgolay pretreated #######################
############################################################
############################################################
############################################################


fullData_paprika <- ssc(fullData, C_type %in% "paprika", include = TRUE) 
fullData_paprika <- reColor(fullData_paprika)
fullData_paprika <- do_sgolay(fullData_paprika,p = 2, n = 51)
fullData_paprika <- selectWls(fullData_paprika, 950,1650)
plot_spectra(fullData_paprika, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paprika_SGOLAYSpec")

fullData_paprika$header$C_all=as.factor(fullData_paprika$header$C_all)
cu <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T,
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_sample","Y_adulterants")))
plot_pls(cu, pg.fns= "PLS_SGOLAYpaprika")

cu2 <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                   pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_sample","Y_adulterants")))
plot_pls(cu2, pg.fns= "PLS_SGOLAYpaprika_MixtureValidation")

#LDA fulldata mixture
dadata <- fullData_paprika$header
dadata$NIR <- fullData_paprika$NIR
dadata$col <- fullData_paprika$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaSGOLAY_mixture.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:90,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprika_mixture_SGOLAY.xlsx", colNames=TRUE)
dev.off()

#mixture but with no controls
fullData_paprika_NoControl <- ssc(fullData_paprika, C_mixture %in% 
                                   c("Control_SH","Control_TU", "Control_M1",
                                     "Control_NA","Control_TY","Control_LE"), include = FALSE) 

dadata <- fullData_paprika_NoControl$header
dadata$NIR <- fullData_paprika_NoControl$NIR
dadata$col <- fullData_paprika_NoControl$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaSGOLAY_mixture_No_control.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:100,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:90,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprikaSGOLAY_mixture_no_control.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paprika$header
dadata$NIR <- fullData_paprika$NIR
dadata$col <- fullData_paprika$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paprikaSGOLAY_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_paprikaSGOLAY_variety.xlsx", colNames=TRUE)
dev.off()


#Lowest adulteration level (0.5%)

fullData_paprika_lowestLev <- ssc(fullData_paprika, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL

cu <- gdmm(fullData_paprika_lowestLev, getap(do.pca= F, do.pls= T, 
                                            pls.regOn =c("Y_corn","Y_color","Y_avocado")))
plot_pls(cu, pg.fns= "PLS_paprikaSGOLAY_LowestLev")

cu2 <- gdmm(fullData_paprika_lowestLev, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                      pls.regOn =c("Y_corn","Y_color","Y_avocado")))
plot_pls(cu2, pg.fns= "PLS_paprikaSGOLAY_MixtureValidation_LowestLev")

#LDA fulldata mixture
dadata <- fullData_paprika_lowestLev$header
dadata$NIR <- fullData_paprika_lowestLev$NIR
dadata$col <- fullData_paprika_lowestLev$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaSGOLAY_mixture_LowestLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:20,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:20,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:20,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprikaSGOLAY_mixture_LowestLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paprika_lowestLev$header
dadata$NIR <- fullData_paprika_lowestLev$NIR
dadata$col <- fullData_paprika_lowestLev$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paprikaSGOLAY_fullData_paprika_lowestLev_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_paprikaSGOLAY_fullData_paprika_lowestLev_variety.xlsx", colNames=TRUE)
dev.off()

############################################################
############################################################
############################################################
########### Sgolay and MSC pretreated #######################
############################################################
############################################################
############################################################

fullData_paprika <- ssc(fullData, C_type %in% "paprika", include = TRUE) 
fullData_paprika <- reColor(fullData_paprika)
fullData_paprika <- do_sgolay(fullData_paprika,p = 2, n = 21)
fullData_paprika <- do_msc(fullData_paprika)
fullData_paprika <- selectWls(fullData_paprika, 950,1650)
plot_spectra(fullData_paprika, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paprika_M&GOLSpec")

cu <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T,
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_sample","Y_adulterants")))
plot_pls(cu, pg.fns= "PLS_M&GOLpaprika")

cu2 <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                   pls.regOn =c("Y_corn","Y_color","Y_avocado" ,"Y_sample","Y_adulterants")))
plot_pls(cu2, pg.fns= "PLS_M&GOLpaprika_MixtureValidation")

#LDA fulldata mixture
dadata <- fullData_paprika$header
dadata$NIR <- fullData_paprika$NIR
dadata$col <- fullData_paprika$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaM&GOL_mixture.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprika_mixture_M&GOL.xlsx", colNames=TRUE)
dev.off()

#mixture but with no controls
fullData_paprika_NoControl <- ssc(fullData_paprika, C_mixture %in% 
                                   c("Control_SH","Control_TU", "Control_M1",
                                     "Control_NA","Control_TY","Control_LE"), include = FALSE) 

dadata <- fullData_paprika_NoControl$header
dadata$NIR <- fullData_paprika_NoControl$NIR
dadata$col <- fullData_paprika_NoControl$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaM&GOL_mixture_No_control.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprikaM&GOL_mixture_no_control.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paprika$header
dadata$NIR <- fullData_paprika$NIR
dadata$col <- fullData_paprika$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paprikaM&GOL_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

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
write.xlsx(conftable , file = "results/fulldata_paprikaM&GOL_variety.xlsx", colNames=TRUE)
dev.off()


#Lowest adulteration level (0.5%)

fullData_paprika_lowestLev <- ssc(fullData_paprika, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL

cu <- gdmm(fullData_paprika_lowestLev, getap(do.pca= F, do.pls= T, 
                                            pls.regOn =c("Y_corn","Y_color","Y_avocado")))
plot_pls(cu, pg.fns= "PLS_paprikaM&GOL_LowestLev")

cu2 <- gdmm(fullData_paprika_lowestLev, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                      pls.regOn =c("Y_corn","Y_color","Y_avocado")))
plot_pls(cu2, pg.fns= "PLS_paprikaM&GOL_MixtureValidation_LowestLev")

#LDA fulldata mixture
dadata <- fullData_paprika_lowestLev$header
dadata$NIR <- fullData_paprika_lowestLev$NIR
dadata$col <- fullData_paprika_lowestLev$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaM&GOL_mixture_LowestLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:40,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:40,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:40,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprikaM&GOL_mixture_LowestLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paprika_lowestLev$header
dadata$NIR <- fullData_paprika_lowestLev$NIR
dadata$col <- fullData_paprika_lowestLev$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paprikaM&GOL_fullData_paprika_lowestLev_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_paprikaM&GOL_fullData_paprika_lowestLev_variety.xlsx", colNames=TRUE)
dev.off()



############################################################
############################################################
############################################################
########### SNV pretreated #######################
############################################################
############################################################
############################################################

fullData_paprika <- ssc(fullData, C_type %in% "paprika", include = TRUE) 
fullData_paprika <- reColor(fullData_paprika)
fullData_paprika <- do_snv(fullData_paprika)
fullData_paprika <- selectWls(fullData_paprika, 950,1650)
plot_spectra(fullData_paprika, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paprika_SNVSpec")

cu <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T,
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_sample","Y_adulterants")))
plot_pls(cu, pg.fns= "PLS_SNVpaprika")

cu2 <- gdmm(fullData_paprika, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                   pls.regOn =c("Y_corn","Y_color","Y_avocado" ,"Y_sample","Y_adulterants")))
plot_pls(cu2, pg.fns= "PLS_SNVpaprika_MixtureValidation")

#LDA fulldata mixture
dadata <- fullData_paprika$header
dadata$NIR <- fullData_paprika$NIR
dadata$col <- fullData_paprika$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaSNV_mixture.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprika_mixture_SNV.xlsx", colNames=TRUE)
dev.off()

#mixture but with no controls
fullData_paprika_NoControl <- ssc(fullData_paprika, C_mixture %in% 
                                   c("Control_SH","Control_TU", "Control_M1",
                                     "Control_NA","Control_TY","Control_LE"), include = FALSE) 

dadata <- fullData_paprika_NoControl$header
dadata$NIR <- fullData_paprika_NoControl$NIR
dadata$col <- fullData_paprika_NoControl$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaSNV_mixture_No_control.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprikaSNV_mixture_no_control.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paprika$header
dadata$NIR <- fullData_paprika$NIR
dadata$col <- fullData_paprika$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paprikaSNV_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_paprikaSNV_variety.xlsx", colNames=TRUE)
dev.off()


#Lowest adulteration level (0.5%)

fullData_paprika_lowestLev <- ssc(fullData_paprika, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL

cu <- gdmm(fullData_paprika_lowestLev, getap(do.pca= F, do.pls= T, 
                                            pls.regOn =c("Y_corn","Y_color","Y_avocado")))
plot_pls(cu, pg.fns= "PLS_paprikaSNV_LowestLev")

cu2 <- gdmm(fullData_paprika_lowestLev, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                      pls.regOn =c("Y_corn","Y_color","Y_avocado"\)))
plot_pls(cu2, pg.fns= "PLS_paprikaSNV_MixtureValidation_LowestLev")

#LDA fulldata mixture
dadata <- fullData_paprika_lowestLev$header
dadata$NIR <- fullData_paprika_lowestLev$NIR
dadata$col <- fullData_paprika_lowestLev$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paprikaSNV_mixture_LowestLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:40,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:40,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:40,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprikaSNV_mixture_LowestLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paprika_lowestLev$header
dadata$NIR <- fullData_paprika_lowestLev$NIR
dadata$col <- fullData_paprika_lowestLev$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paprikaSNV_fullData_paprika_lowestLev_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomleft", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("topright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paprikaSNV_fullData_paprika_lowestLev_variety.xlsx", colNames=TRUE)
dev.off()