#This are instructions for outlier detection (do not run here but in main code)
#################################################################################
################################################################################
#OUTLIER DETECTION for LDA models
#1) Run the model and visually detect if the outliers are in the training or validation set 
#2) Create an identifier for the data to be analized. e.g. C_ID 
fullDataLDA$header$C_ID <- seq(1, nrow(fullDataLDA), 1)
fullDataLDA <- reColor(fullDataLDA)
#3)Identify the outliers(clicking)
#run j=(1, 2, 3, etc) hasta plotFDA, hasta que se cargue el da1 correcto 
plotFDAmodel(da1, projCVres = TRUE, col = sampleColor[tri], col2 = sampleColor[-tri], pch = 16 , pch2 = 1)
identify(da1$scoorsVal)
#4)Identify the outlayers in the data (depending if it is training or validation)
fullData$header[-tri,][c(1:3,101:10)] ##to identify the outlayers in the LDA plot 
#5)kick out the outlayers
indexout <- which(fullDataLDA$header$C_ID %in% c(430:432, 538:540)) ##to kick the outlayers based in the identifier C_ID
fullDataLDA1 <- fullDataLDA[-indexout,]
############################################################################
###########################################################################