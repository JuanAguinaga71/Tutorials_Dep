library(aquap2)
library(openxlsx)
library(mda)
source("R-code/FDAnir.R")

load("R-data/whey_metri.R")
whey_metri <- ssc(whey_metri, C_adultLev %in% c("0.5% Adulteration","Pure Whey"), include = TRUE)
whey_metri <- ssc(whey_metri, C_adultLev %in% c("1% Adulteration","Pure Whey"), include = TRUE)
whey_metri <- ssc(whey_metri, C_adultLev %in% c("1.5% Adulteration","Pure Whey"), include = TRUE)
whey_metri <- ssc(whey_metri, C_adultLev %in% c("2% Adulteration","Pure Whey"), include = TRUE)
whey_metri <- ssc(whey_metri, C_adultLev %in% c("2.5% Adulteration","Pure Whey"), include = TRUE)
whey_metri <- ssc(whey_metri, C_adultLev %in% c("3% Adulteration","Pure Whey"), include = TRUE)

#Mixtures
whey_metri <- ssc(whey_metri, C_mixture %in% c("U","G","T","M","Pure Whey"), include = TRUE)
whey_metri <- ssc(whey_metri, C_mixture %in% c("UG","UM","UT","GT","TM","GM"), include = TRUE)
whey_metri <- ssc(whey_metri, C_mixture %in% c("UGT","UGM","UTM","GTM","UGTM"), include = TRUE)


#Classification
#Mixtures
dadata <- whey_metri$header
dadata$NIR <- whey_metri$NIR
dadata$col <- whey_metri$colRep

col = dadata$col$C_mixture
pdf("results/AllmixLev1.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:30,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.6)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:30,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("bottomleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.9)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_mixture", NrPCs = 1:30,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_mixture), col = unique(col), pch = 16, cex = 0.9)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/AllmixLeve1.xlsx", colNames=TRUE)
dev.off()


# Adultlev
dadata <- whey_metri$header
dadata$NIR <- whey_metri$NIR
dadata$col <- whey_metri$colRep
col = dadata$col$C_adultLev

pdf("results/AllLev_Metri.pdf")
TranInd = -seq(1, nrow(dadata), 3)
da1 <- FDA_PCAmodels(dadata,"C_adultLev", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da1, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_adultLev), col = unique(col), pch = 16, cex = 0.8)

TranInd = -seq(2, nrow(dadata), 3)
da2 <- FDA_PCAmodels(dadata,"C_adultLev", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da2, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_adultLev), col = unique(col), pch = 16, cex = 0.5)

TranInd = -seq(3, nrow(dadata), 3)
da3 <- FDA_PCAmodels(dadata,"C_adultLev", NrPCs = 1:60,  TranInd = TranInd)
plotFDAmodel(da3, projCVres = TRUE, col = col[TranInd], col2 = col[-TranInd], main = "Discriminant analysis of mixtures")
legend("topleft", legend= unique(dadata$C_adultLev), col = unique(col), pch = 16, cex = 0.5)

conftable <- calcAveConfMatrix(da1, da2, da3)
conftable
write.xlsx(conftable , file = "results/All_Adultlevmixtures_metri.xlsx", colNames=TRUE)
dev.off()



#PLSR for Whey metri
cu <- gdmm(whey_metri, getap(do.pca= F, do.pls= T, do.sim= F, sim.vars= c("C_urea","C_taurine","C_glycine","C_melamine"),
                             do.svm= F, svm.classOn= c("C_adultLev","C_mixture"), do.nnet= F, nnet.classOn= c("C_adultLev","C_mixture"),
                             pls.regOn = c("Y_urea","Y_taurine","Y_glycine","Y_melamine","Y_Protein")))

plot_pls(cu, pg.fns= "Metri_TripQuadAdult_Alllevels")
plot_simca(cu, pg.fns= "Metri")
plot_svm(cu, pg.fns= "Metri")


#PLSR  independent test for metri
whey_metri_calib <- ssc(whey_metri, C_Repl == "R1"|C_Repl == "R2", include = TRUE) # create the calibration set
whey_metri_test <- ssc(whey_metri, C_Repl == "R1"|C_Repl == "R2", include = FALSE) # create the test set

cu_indepPred <- gdmm(whey_metria_calib, getap(do.pca = F, pca.colorBy = cb, pls.regOn = c("Y_urea","Y_taurine","Y_glycine","Y_melamine"))) # here you build your plsr model
plot_pls_indepPred(whey_metri_test, cu_indepPred,pg.fns = "PLSR_IndPredic_fulldata")


#Caliberation transfer
cu <- gdmm(whey_metri, getap(do.pca= F, do.pls= T, 
                             pls.regOn = c("Y_urea","Y_taurine","Y_glycine","Y_melamine","Y_Protein")))
cu1 <- gdmm(Nano_no_bag, getap(do.pca= F, do.pls= T, 
                             pls.regOn = c("Y_urea","Y_taurine","Y_glycine","Y_melamine","Y_Protein")))



metri_wavelengths <- aquap2::getWavelengths(whey_metri)
metri_wavelengths <- seq(952,1648,2)
resampleNIR_inner <- function(NIR, targetWls=metri_wavelengths, method="linear") {
  x <- aquap2::getWavelengths(NIR)
  if (is.null(targetWls)) {
    xNew <- seq(ceiling(x[1]/2) * 2, floor(x[length(x)]/2) * 2, 2)
  } else {
    xNew <- targetWls
  }
  NIRnew <- t(apply(NIR$NIR, 1, pracma::interp1, x = x, xi = xNew, method = method))
  colnames(NIRnew) <- paste0("X", xNew)
  return(NIRnew)
} # EOF
Nano_no_bagX <- resampleNIR_inner(Nano_no_bag)

Nano_no_bag_resampled <- Nano_no_bag

Nano_no_bag_resampled$NIR <- Nano_no_bagX

plot_spectra(Nano_no_bag_resampled)
