## 
##
##

library(ERP)

dataIn = read.table("/PEPS-protocol-phase2/Segmented_Data/myTextFile_allsuj.txt", header = FALSE, sep = "", dec = ".")
time_pt<-seq(from=0, to= 1050, by=1.9531)
dim(dataIn)

covariates<-dataIn[,1:4]
erpdata<-dataIn[,-(1:4)]

##****************************************************************************************************
# Extract the data for Human congruent-incong for all midline electrodes
dataIn.congincongFCz <- subset(dataIn, 
                           (dataIn$V2 == "FCz") & (dataIn$V4 == "Congruent-Incong"))
dataIn.congincongFz <- subset(dataIn, 
                            (dataIn$V2 == "Fz") & (dataIn$V4 == "Congruent-Incong"))
dataIn.congincongCz <- subset(dataIn, 
                               (dataIn$V2 == "Cz") & (dataIn$V4 == "Congruent-Incong"))
dataIn.congincongPz <- subset(dataIn, 
                              (dataIn$V2 == "Pz") & (dataIn$V4 == "Congruent-Incong"))

# Extract the data for Agent-Human congruent for all midline electrodes
dataIn.congFCz <- subset(dataIn, 
                               (dataIn$V2 == "FCz") & (dataIn$V4 == "Congruent"))
dataIn.congFz <- subset(dataIn, 
                              (dataIn$V2 == "Fz") & (dataIn$V4 == "Congruent"))
dataIn.congCz <- subset(dataIn, 
                              (dataIn$V2 == "Cz") & (dataIn$V4 == "Congruent"))
dataIn.congPz <- subset(dataIn, 
                              (dataIn$V2 == "Pz") & (dataIn$V4 == "Congruent"))

# Extract the data for Agent-Human Incongruent for all midline electrodes
dataIn.incongFCz <- subset(dataIn, 
                         (dataIn$V2 == "FCz") & (dataIn$V4 == "Incongruent"))
dataIn.incongFz <- subset(dataIn, 
                        (dataIn$V2 == "Fz") & (dataIn$V4 == "Incongruent"))
dataIn.incongCz <- subset(dataIn, 
                        (dataIn$V2 == "Cz") & (dataIn$V4 == "Incongruent"))
dataIn.incongPz <- subset(dataIn, 
                        (dataIn$V2 == "Pz") & (dataIn$V4 == "Incongruent"))

# Define the erp-data and the covariates for this effect of interest - congruent-incong
covariates.congincongFCz <- dataIn.congincongFCz[, 1:4] 
erpdata.congincongFCz <- dataIn.congincongFCz[, -(1:4)]  

covariates.congincongFz <- dataIn.congincongFz[, 1:4] 
erpdata.congincongFz <- dataIn.congincongFz[, -(1:4)] 

covariates.congincongCz <- dataIn.congincongCz[, 1:4] 
erpdata.congincongCz <- dataIn.congincongCz[, -(1:4)] 

covariates.congincongPz <- dataIn.congincongPz[, 1:4] 
erpdata.congincongPz <- dataIn.congincongPz[, -(1:4)] 

# Define the erp-data and the covariates for this effect of interest - congruent
covariates.congFCz <- dataIn.congFCz[, 1:4] 
erpdata.congFCz <- dataIn.congFCz[, -(1:4)]  

covariates.congFz <- dataIn.congFz[, 1:4] 
erpdata.congFz <- dataIn.congFz[, -(1:4)] 

covariates.congCz <- dataIn.congCz[, 1:4] 
erpdata.congCz <- dataIn.congCz[, -(1:4)] 

covariates.congPz <- dataIn.congPz[, 1:4] 
erpdata.congPz <- dataIn.congPz[, -(1:4)] 

# Define the erp-data and the covariates for this effect of interest - incongruent
covariates.incongFCz <- dataIn.incongFCz[, 1:4] 
erpdata.incongFCz <- dataIn.incongFCz[, -(1:4)]  

covariates.incongFz <- dataIn.incongFz[, 1:4] 
erpdata.incongFz <- dataIn.incongFz[, -(1:4)] 

covariates.incongCz <- dataIn.incongCz[, 1:4] 
erpdata.incongCz <- dataIn.incongCz[, -(1:4)] 

covariates.incongPz <- dataIn.incongPz[, 1:4] 
erpdata.incongPz <- dataIn.incongPz[, -(1:4)] 



par(mar = rep(2,2))
par(mfrow=c(2,2))

# Set up design matrix to test for Human vs. Agent effect for congruent-incong condition
design <- model.matrix(~ C(V1, sum) + relevel(V3, ref = "Agent"), data = covariates.congincongPz)
design0 <-model.matrix(~ C(V1, sum), data = covariates.congincongPz)

# We will average the curves over a pre-defined number of time points (by default 10)
# The Benjamini-Hochberg (BH) procedure is applied by default to ensure that the fdr
# is controlled at a preset level alpha (by default alpha = .05).
avetest <- erpavetest(erpdata.congincongPz, design, design0, nintervals = 20, method = "none", alpha=.05)

# Plot the difference curve and show the significant time points.
erpplot(erpdata.congincongPz, design, effect = ncol(design), lwd = 2, interval = "simultaneous",
        frames = time_pt, y = rev(c(-5, 5)), xlab = "Time", ylab = "Condition Effect")
title(main = "Agent - Human diff curve \n Congruent video incong, Pz", font.main = 2)
abline(v=250, col="red")
points(time_pt[avetest$significant], rep(0, length(avetest$significant)), pch = 20, col = "goldenrod")
abline(v = time_pt[avetest$breaks], lty = 2, col = "darkgray")


# Applying the Adaptive Factor-Adjustment Method --------------------------

fabh_Fz <- erpfatest(erpdata.incongFz, design, design0, nbf = NULL, wantplot = TRUE)
fabh2_Fz <- erpfatest(erpdata.incongFz, design, design0, nbf = fabh$nbf)

fabh_FCz <- erpfatest(erpdata.incongFCz, design, design0, nbf = NULL, wantplot = TRUE)
fabh2_FCz <- erpfatest(erpdata.incongFCz, design, design0, nbf = fabh$nbf)

fabh_Cz <- erpfatest(erpdata.incongCz, design, design0, nbf = NULL, wantplot = TRUE)
fabh2_Cz <- erpfatest(erpdata.incongCz, design, design0, nbf = fabh$nbf)

fabh_Pz <- erpfatest(erpdata.incongFCz, design, design0, nbf = NULL, wantplot = TRUE)
fabh2_Pz <- erpfatest(erpdata.incongFCz, design, design0, nbf = fabh$nbf)

par(mar = rep(2,2))
par(mfrow=c(2,2))

erpplot(erpdata.incongPz, design, effect = ncol(design), lwd = 2,
        interval = "simultaneous", frames = time_pt, ylim = rev(c(-6, 6)),
        xlab = "Time (ms)", ylab = "Condition effect curve")
title("Agent - Human Diff Curve \n Incongruent Pz", font = 2)
abline(v=250, col="red")
points(time_pt[fabh2_Pz$significant], rep(0, length(fabh2_Pz$significant)),
       pch = 20, col = "goldenrod")