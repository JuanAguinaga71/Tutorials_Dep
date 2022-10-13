# extract some error values ----------------------------------------


#' @title Extract RMSEC of the PLSR model
#' @description Extracts the Root Mean Square Error of Calibration of the given PLSR model
#' @details See examples
#' @param plsModel The PLSR model calculated by gdmm when do.pls = TRUE.
#' @return The RMSEC of the PLSR model.
#' @examples 
#' \dontrun{
#' fd <- gfd()
#' cu <- gdmm(fd)
#' plsModel <- getcm(cu, 1, what = "plsr")
#' RMSEC <- getRMSEC(plsModel)
#' }
#' @seealso \code{\link{plot_pls-method}}
#' @family PLSR documentation
#' @export
getRMSEC <- function(plsModel) {
  a <- c(pls::RMSEP(plsModel, estimate="train")$val)
  b <- a[length(a)]
  out <- round(b, .ap2$stn$plsr_nrDigitsRMSEx)
} #EOF


#' @title Extract RMSECV of the PLSR model
#' @description Extracts the Root Mean Square Error of Cross-validation of the given PLSR model
#' @details See examples
#' @param plsModel The PLSR model calculated by gdmm when do.pls = TRUE.
#' @return The RMSECV of the PLSR model.
#' @examples 
#' \dontrun{
#' fd <- gfd()
#' cu <- gdmm(fd)
#' plsModel <- getcm(cu, 1, what = "plsr")
#' RMSECV <- getRMSECV(plsModel)
#' }
#' @seealso \code{\link{plot_pls-method}}
#' @family PLSR documentation
#' @export
getRMSECV <- function(plsModel) {
  a <- c(pls::RMSEP(plsModel, estimate="CV")$val)
  b <- a[length(a)]
  out <- round(b, .ap2$stn$plsr_nrDigitsRMSEx)
} #EOF


#' @title Extract RMSEP of the PLSR model
#' @description Extracts the Root Mean Square Error of Prediction of the given PLSR model
#' @details See examples
#' @param plsModel The PLSR model calculated by gdmm when do.pls = TRUE.
#' @param newdata a dataframe containing the spectral data used for independent prediction
#' @return The RMSEP of the PLSR model.
#' @examples 
#' \dontrun{
#' fd <- gfd()
#' cu <- gdmm(fd)
#' plsModel <- getcm(cu, 1, what = "plsr")
#' RMSEP <- getRMSEP(plsModel)
#' }
#' @seealso \code{\link{plot_pls-method}}
#' @family PLSR documentation
#' @export
getRMSEP <- function(plsModel, newdata) {
  a <- c(pls::RMSEP(plsModel, estimate="test", newdata=newdata)$val)
  b <- a[length(a)]
  out <- round(b, .ap2$stn$plsr_nrDigitsRMSEx)
} # EOF


#' @title Extract R-square of calibration of the PLSR model
#' @description Extracts the determination coefficient of calibration of the given PLSR model
#' @details See examples
#' @param plsModel The PLSR model calculated by gdmm when do.pls = TRUE.
#' @return The R^2 of calibration of the PLSR model.
#' @examples 
#' \dontrun{
#' fd <- gfd()
#' cu <- gdmm(fd)
#' plsModel <- getcm(cu, 1, what = "plsr")
#' R2C <- getR2C(plsModel)
#' }
#' @seealso \code{\link{plot_pls-method}}
#' @family PLSR documentation
#' @export
getR2C <- function(plsModel) {
  a <- c(pls::R2(plsModel, estimate="train")$val)
  b <- a[length(a)]
  out <- round(b, .ap2$stn$plsr_nrDigitsRMSEx)
} #EOF


#' @title Extract R-square of cross-validation of the PLSR model
#' @description Extracts the determination coefficient of cross-validation of the given PLSR model
#' @details See examples
#' @param plsModel The PLSR model calculated by gdmm when do.pls = TRUE.
#' @return The R^2 of cross-validation of the PLSR model.
#' @examples 
#' \dontrun{
#' fd <- gfd()
#' cu <- gdmm(fd)
#' plsModel <- getcm(cu, 1, what = "plsr")
#' R2CV <- getR2CV(plsModel)
#' }
#' @seealso \code{\link{plot_pls-method}}
#' @family PLSR documentation
#' @export
getR2CV <- function(plsModel) {
  a <- c(pls::R2(plsModel, estimate="CV")$val)
  b <- a[length(a)]
  out <- round(b, .ap2$stn$plsr_nrDigitsRMSEx)
} #EOF


#' @title Extract R-square of prediction of the PLSR model
#' @description Extracts the determination coefficient of prediction of the given PLSR model
#' @details See examples
#' @param plsModel The PLSR model calculated by gdmm when do.pls = TRUE.
#' @param newdata a dataframe containing the spectral data used for independent prediction
#' @return The R^2 of prediction of the PLSR model.
#' @examples 
#' \dontrun{
#' fd <- gfd()
#' cu <- gdmm(fd)
#' plsModel <- getcm(cu, 1, what = "plsr")
#' R2P <- getR2P(plsModel)
#' }
#' @seealso \code{\link{plot_pls-method}}
#' @family PLSR documentation
#' @export
getR2P <- function(plsModel, newdata) {
  a <- c(pls::R2(plsModel, estimate="test", newdata=newdata)$val)
  b <- a[length(a)]
  out <- round(b, .ap2$stn$plsr_nrDigitsRMSEx)
} # EOF
####
getVecRMSEC <- function(plsModel) {			## for error plot
  out <- (c(pls::RMSEP(plsModel, estimate="train", intercept=TRUE)$val[,1,]))
} # EOF

getVecRMSECV <- function(plsModel) {		## for error plot
  out <- (c(pls::RMSEP(plsModel, estimate="CV", intercept=TRUE)$val[,1,]))
} # EOF

getVecRMSECV_adj <- function(plsModel) {	## for error plot
  out <- (c(pls::RMSEP(plsModel, estimate="adjCV", intercept=TRUE)$val[,1,]))
} # EOF
####
convertToRDP <- function(errorValue, ClassVar, header) {	# the dataset with all the observations that form the range	
  sdY <- sd(header[,ClassVar], na.rm=TRUE)
  out <- round(sdY/errorValue, .ap2$stn$plsr_nrDigitsRMSEx)
} # EOF
