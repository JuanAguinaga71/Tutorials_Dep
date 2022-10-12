library(aquap2)
library(openxlsx)
library(mda)
source("R-code/FDAnir.R")



#fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = "custom@importTRH_VoltcraftFizAut_Lew.R", slType = "xls", multiplyRows = F, dol = F)
fullData <- gfd(filetype = "custom@MetriExcelTimeStImp.R", trhLog = F, slType = "xls", multiplyRows = F, dol = F)
fullData <- reColor(fullData)

cb<- c("C_samplename","C_type" ,"C_mixture", "C_variety", "C_adultlev", "C_adultNr","C_conSNr","C_Repl")

#plot_spectra(fullData, colorBy = cb, pg.fns= "rawSpec")

#Fulldata
fullData_paptom <- ssc(fullData, C_type %in% c("paprika","tomato"), include = TRUE) 
fullData_paptom <- reColor(fullData_paptom)
fullData_paptom <- selectWls(fullData_paptom, 950,1650)
plot_spectra(fullData_paptom, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paptom_rawSpec")

cu <- gdmm(fullData_paptom, getap(do.pca= F, do.pls= T, 
                          pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pls(cu, pg.fns= "PLS_paptom")

cu2 <- gdmm(fullData_paptom, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pls(cu2, pg.fns= "PLS_paptom_MixtureValidation")

#LDA fulldata mixture
dadata <- fullData_paptom$header
dadata$NIR <- fullData_paptom$NIR
dadata$col <- fullData_paptom$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_tompap_mixture.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_tompap_mixture.xlsx", colNames=TRUE)
dev.off()

#mixture but with no controls
fullData_paptom_NoControl <- ssc(fullData_paptom, C_mixture %in% 
                                   c("Control_SH","Control_TU", "Control_M1",
                                     "Control_NA","Control_TY","Control_LE"), include = FALSE) 

dadata <- fullData_paptom_NoControl$header
dadata$NIR <- fullData_paptom_NoControl$NIR
dadata$col <- fullData_paptom_NoControl$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_tompap_mixture_No_control.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:70,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:70,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:70,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_tompap_mixture_no_control.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paptom$header
dadata$NIR <- fullData_paptom$NIR
dadata$col <- fullData_paptom$colRep
col = dadata$col$C_variety

pdf("results/fulldata_tompap_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_tompap_variety.xlsx", colNames=TRUE)
dev.off()


#Lowest adulteration level (0.5%)

fullData_paptom_lowestLev <- ssc(fullData_paptom, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL

cu <- gdmm(fullData_lowestLev, getap(do.pca= F, do.pls= T, 
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto")))
plot_pls(cu, pg.fns= "PLS_paptom_LowestLev")

cu2 <- gdmm(fullData_lowestLev, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                   pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto")))
plot_pls(cu2, pg.fns= "PLS_paptom_MixtureValidation_LowestLev")

#LDA fulldata mixture
dadata <- fullData_paptom_lowestLev$header
dadata$NIR <- fullData_paptom_lowestLev$NIR
dadata$col <- fullData_paptom_lowestLev$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_tompap_mixture_LowestLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:70,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:70,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_tompap_mixture_LowestLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paptom_lowestLev$header
dadata$NIR <- fullData_paptom_lowestLev$NIR
dadata$col <- fullData_paptom_lowestLev$colRep
col = dadata$col$C_variety

pdf("results/fulldata_tompap_fullData_lowestLev_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_tompap_fullData_lowestLev_variety.xlsx", colNames=TRUE)
dev.off()


############################################################
############################################################
############################################################
########### MSC PRETREATED #######################
############################################################
############################################################
############################################################

fullData_paptom <- ssc(fullData, C_type %in% c("paprika","tomato"), include = TRUE) 
fullData_paptom <- reColor(fullData_paptom)
fullData_paptom <- do_msc(fullData_paptom)
fullData_paptom <- selectWls(fullData_paptom, 950,1650)
plot_spectra(fullData_paptom, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paptom_MSCSpec")

cu <- gdmm(fullData_paptom, getap(do.pca= F, do.pls= T,
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pls(cu, pg.fns= "PLS_MSCpaptom")

cu2 <- gdmm(fullData_paptom, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                   pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pls(cu2, pg.fns= "PLS_MSCpaptom_MixtureValidation")

#LDA fulldata mixture
dadata <- fullData_paptom$header
dadata$NIR <- fullData_paptom$NIR
dadata$col <- fullData_paptom$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomMSC_mixture.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:90,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptom_mixture_MSC.xlsx", colNames=TRUE)
dev.off()

#mixture but with no controls
fullData_paptom_NoControl <- ssc(fullData_paptom, C_mixture %in% 
                                   c("Control_SH","Control_TU", "Control_M1",
                                     "Control_NA","Control_TY","Control_LE"), include = FALSE) 

dadata <- fullData_paptom_NoControl$header
dadata$NIR <- fullData_paptom_NoControl$NIR
dadata$col <- fullData_paptom_NoControl$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomMSC_mixture_No_control.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:100,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:90,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomMSC_mixture_no_control.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paptom$header
dadata$NIR <- fullData_paptom$NIR
dadata$col <- fullData_paptom$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paptomMSC_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomMSC_variety.xlsx", colNames=TRUE)
dev.off()


#Lowest adulteration level (0.5%)

fullData_paptom_lowestLev <- ssc(fullData_paptom, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL

cu <- gdmm(fullData_paptom_lowestLev, getap(do.pca= F, do.pls= T, 
                                     pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto")))
plot_pls(cu, pg.fns= "PLS_paptomMSC_LowestLev")

cu2 <- gdmm(fullData_lowestLev, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                      pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto")))
plot_pls(cu2, pg.fns= "PLS_paptomMSC_MixtureValidation_LowestLev")

#LDA fulldata mixture
dadata <- fullData_paptom_lowestLev$header
dadata$NIR <- fullData_paptom_lowestLev$NIR
dadata$col <- fullData_paptom_lowestLev$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomMSC_mixture_LowestLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomMSC_mixture_LowestLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_lowestLev$header
dadata$NIR <- fullData_lowestLev$NIR
dadata$col <- fullData_lowestLev$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paptomMSC_fullData_lowestLev_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_paptomMSC_fullData_lowestLev_variety.xlsx", colNames=TRUE)
dev.off()




############################################################
############################################################
############################################################
########### Sgolay pretreated #######################
############################################################
############################################################
############################################################


fullData_paptom <- ssc(fullData, C_type %in% c("paprika","tomato"), include = TRUE) 
fullData_paptom <- reColor(fullData_paptom)
fullData_paptom <- do_sgolay(fullData_paptom,p = 2, n = 21)
fullData_paptom <- selectWls(fullData_paptom, 950,1650)
plot_spectra(fullData_paptom, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paptom_SGOLAYSpec")

cu <- gdmm(fullData_paptom, getap(do.pca= F, do.pls= T,
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pls(cu, pg.fns= "PLS_SGOLAYpaptom")

cu2 <- gdmm(fullData_paptom, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                   pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pls(cu2, pg.fns= "PLS_SGOLAYpaptom_MixtureValidation")

#LDA fulldata mixture
dadata <- fullData_paptom$header
dadata$NIR <- fullData_paptom$NIR
dadata$col <- fullData_paptom$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomSGOLAY_mixture.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:90,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptom_mixture_SGOLAY.xlsx", colNames=TRUE)
dev.off()

#mixture but with no controls
fullData_paptom_NoControl <- ssc(fullData_paptom, C_mixture %in% 
                                   c("Control_SH","Control_TU", "Control_M1",
                                     "Control_NA","Control_TY","Control_LE"), include = FALSE) 

dadata <- fullData_paptom_NoControl$header
dadata$NIR <- fullData_paptom_NoControl$NIR
dadata$col <- fullData_paptom_NoControl$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomSGOLAY_mixture_No_control.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:100,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:90,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomSGOLAY_mixture_no_control.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paptom$header
dadata$NIR <- fullData_paptom$NIR
dadata$col <- fullData_paptom$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paptomSGOLAY_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomSGOLAY_variety.xlsx", colNames=TRUE)
dev.off()


#Lowest adulteration level (0.5%)

fullData_paptom_lowestLev <- ssc(fullData_paptom, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL

cu <- gdmm(fullData_paptom_lowestLev, getap(do.pca= F, do.pls= T, 
                                            pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto")))
plot_pls(cu, pg.fns= "PLS_paptomSGOLAY_LowestLev")

cu2 <- gdmm(fullData_lowestLev, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                      pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto")))
plot_pls(cu2, pg.fns= "PLS_paptomSGOLAY_MixtureValidation_LowestLev")

#LDA fulldata mixture
dadata <- fullData_paptom_lowestLev$header
dadata$NIR <- fullData_paptom_lowestLev$NIR
dadata$col <- fullData_paptom_lowestLev$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomSGOLAY_mixture_LowestLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomSGOLAY_mixture_LowestLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_lowestLev$header
dadata$NIR <- fullData_lowestLev$NIR
dadata$col <- fullData_lowestLev$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paptomSGOLAY_fullData_lowestLev_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_paptomSGOLAY_fullData_lowestLev_variety.xlsx", colNames=TRUE)
dev.off()

############################################################
############################################################
############################################################
########### Sgolay and MSC pretreated #######################
############################################################
############################################################
############################################################

fullData_paptom <- ssc(fullData, C_type %in% c("paprika","tomato"), include = TRUE) 
fullData_paptom <- reColor(fullData_paptom)
fullData_paptom <- do_sgolay(fullData_paptom,p = 2, n = 21)
fullData_paptom <- do_msc(fullData_paptom)
fullData_paptom <- selectWls(fullData_paptom, 950,1650)
plot_spectra(fullData_paptom, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paptom_M&GOLSpec")

cu <- gdmm(fullData_paptom, getap(do.pca= F, do.pls= T,
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pls(cu, pg.fns= "PLS_M&GOLpaptom")

cu2 <- gdmm(fullData_paptom, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                   pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pls(cu2, pg.fns= "PLS_M&GOLpaptom_MixtureValidation")

#LDA fulldata mixture
dadata <- fullData_paptom$header
dadata$NIR <- fullData_paptom$NIR
dadata$col <- fullData_paptom$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomM&GOL_mixture.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptom_mixture_M&GOL.xlsx", colNames=TRUE)
dev.off()

#mixture but with no controls
fullData_paptom_NoControl <- ssc(fullData_paptom, C_mixture %in% 
                                   c("Control_SH","Control_TU", "Control_M1",
                                     "Control_NA","Control_TY","Control_LE"), include = FALSE) 

dadata <- fullData_paptom_NoControl$header
dadata$NIR <- fullData_paptom_NoControl$NIR
dadata$col <- fullData_paptom_NoControl$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomM&GOL_mixture_No_control.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomM&GOL_mixture_no_control.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paptom$header
dadata$NIR <- fullData_paptom$NIR
dadata$col <- fullData_paptom$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paptomM&GOL_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomM&GOL_variety.xlsx", colNames=TRUE)
dev.off()


#Lowest adulteration level (0.5%)

fullData_paptom_lowestLev <- ssc(fullData_paptom, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL

cu <- gdmm(fullData_paptom_lowestLev, getap(do.pca= F, do.pls= T, 
                                            pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto")))
plot_pls(cu, pg.fns= "PLS_paptomM&GOL_LowestLev")

cu2 <- gdmm(fullData_lowestLev, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                      pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto")))
plot_pls(cu2, pg.fns= "PLS_paptomM&GOL_MixtureValidation_LowestLev")

#LDA fulldata mixture
dadata <- fullData_paptom_lowestLev$header
dadata$NIR <- fullData_paptom_lowestLev$NIR
dadata$col <- fullData_paptom_lowestLev$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomM&GOL_mixture_LowestLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomM&GOL_mixture_LowestLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_lowestLev$header
dadata$NIR <- fullData_lowestLev$NIR
dadata$col <- fullData_lowestLev$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paptomM&GOL_fullData_lowestLev_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_paptomM&GOL_fullData_lowestLev_variety.xlsx", colNames=TRUE)
dev.off()



############################################################
############################################################
############################################################
########### SNV pretreated #######################
############################################################
############################################################
############################################################

fullData_paptom <- ssc(fullData, C_type %in% c("paprika","tomato"), include = TRUE) 
fullData_paptom <- reColor(fullData_paptom)
fullData_paptom <- do_snv(fullData_paptom)
fullData_paptom <- selectWls(fullData_paptom, 950,1650)
plot_spectra(fullData_paptom, colorBy = c("C_mixture","C_variety", "C_adultlev", "C_adultNr"), pg.fns= "paptom_SNVSpec")

cu <- gdmm(fullData_paptom, getap(do.pca= F, do.pls= T,
                                  pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pls(cu, pg.fns= "PLS_SNVpaptom")

cu2 <- gdmm(fullData_paptom, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                   pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto" ,"Y_sample","Y_adulterants")))
plot_pls(cu2, pg.fns= "PLS_SNVpaptom_MixtureValidation")

#LDA fulldata mixture
dadata <- fullData_paptom$header
dadata$NIR <- fullData_paptom$NIR
dadata$col <- fullData_paptom$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomSNV_mixture.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptom_mixture_SNV.xlsx", colNames=TRUE)
dev.off()

#mixture but with no controls
fullData_paptom_NoControl <- ssc(fullData_paptom, C_mixture %in% 
                                   c("Control_SH","Control_TU", "Control_M1",
                                     "Control_NA","Control_TY","Control_LE"), include = FALSE) 

dadata <- fullData_paptom_NoControl$header
dadata$NIR <- fullData_paptom_NoControl$NIR
dadata$col <- fullData_paptom_NoControl$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomSNV_mixture_No_control.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.8)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomSNV_mixture_no_control.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_paptom$header
dadata$NIR <- fullData_paptom$NIR
dadata$col <- fullData_paptom$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paptomSNV_variety.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_variety", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of variety")
legend("bottomright", legend= unique(dadata$C_variety), col = unique(col), pch = 16, cex = 0.8)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomSNV_variety.xlsx", colNames=TRUE)
dev.off()


#Lowest adulteration level (0.5%)

fullData_paptom_lowestLev <- ssc(fullData_paptom, C_adultlev %in% 1, include = TRUE) #LOWEST ADULTERATION LEVEL

cu <- gdmm(fullData_paptom_lowestLev, getap(do.pca= F, do.pls= T, 
                                            pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto")))
plot_pls(cu, pg.fns= "PLS_paptomSNV_LowestLev")

cu2 <- gdmm(fullData_lowestLev, getap(do.pca= F, do.pls= T, pls.valid= c("C_mixture"),
                                      pls.regOn =c("Y_corn","Y_color","Y_avocado","Y_Anatto")))
plot_pls(cu2, pg.fns= "PLS_paptomSNV_MixtureValidation_LowestLev")

#LDA fulldata mixture
dadata <- fullData_paptom_lowestLev$header
dadata$NIR <- fullData_paptom_lowestLev$NIR
dadata$col <- fullData_paptom_lowestLev$colRep
col = dadata$col$C_mixture

pdf("results/fulldata_paptomSNV_mixture_LowestLev.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:80,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topright", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)


conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/fulldata_paptomSNV_mixture_LowestLev.xlsx", colNames=TRUE)
dev.off()

#LDA fulldata variety
dadata <- fullData_lowestLev$header
dadata$NIR <- fullData_lowestLev$NIR
dadata$col <- fullData_lowestLev$colRep
col = dadata$col$C_variety

pdf("results/fulldata_paptomSNV_fullData_lowestLev_variety.pdf")
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
write.xlsx(conftable , file = "results/fulldata_paptomSNV_fullData_lowestLev_variety.xlsx", colNames=TRUE)
dev.off()