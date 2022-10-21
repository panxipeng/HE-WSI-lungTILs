##########################################
# TimeROC TILs    
# TTLs projects 

##########################################

#lib
library(timeROC)
library(survival)
library(risksetROC)
library(openxlsx)


# 0. load the dataset
shengyi.data<-read.xlsx('shengyi_twice.xlsx')


# 1. median cut-off value
temp = shengyi.data$inflam_per_tumor > median(na.omit(shengyi.data$inflam_per_tumor))
inflam_per_tumor_24=as.vector(ifelse(temp,"high","low"))
shengyi.data$inflam_per_tumor_24 <- inflam_per_tumor_24
temp = shengyi.data$inflam_per_stromal > median(na.omit(shengyi.data$inflam_per_stromal))
inflam_per_stromal_24=as.vector(ifelse(temp,"high","low"))
shengyi.data$inflam_per_stromal_24 <- inflam_per_stromal_24


# 2. change stage to binary
tmp1 = shengyi.data$AJCC_stage
tmp1 = ifelse((tmp1 == 3),1,0)
shengyi.data$AJCC_stage_bin <- factor(tmp1)


# 3. factors
for (i in 3:13)
{
  fvars = colnames(shengyi.data)[i]
  shengyi.data[fvars]<-lapply(shengyi.data[fvars],factor)
}
for (i in 32:length(shengyi.data))
{
  fvars = colnames(shengyi.data)[i]
  shengyi.data[fvars]<-lapply(shengyi.data[fvars],factor)
}


# 4.only our features(oof)
# shengyi
tmp1 = shengyi.data$inflam_per_tumor
shengyi.data$inflam_per_tumor_log <- log(tmp1 + 1)
tmp1 = shengyi.data$inflam_per_stromal
shengyi.data$inflam_per_stromal_log <- log(tmp1 + 1)
f24_ooftwo_sy <- coxph(Surv(OS_month, OS_yn) ~ inflam_per_tumor_log + inflam_per_stromal_log, data = shengyi.data)
coxsummary=summary(f24_ooftwo_sy) 
coef6 = as.numeric(coxsummary$coefficients[,"coef"][1])
coef7 = as.numeric(coxsummary$coefficients[,"coef"][2])
score_shengyi_2factors = coef6*shengyi.data$inflam_per_tumor_log + coef7*shengyi.data$inflam_per_stromal_log
temp = score_shengyi_2factors > median(na.omit(score_shengyi_2factors))
# score_24_2factors=as.vector(ifelse(temp,"high","low"))
score_24_2factors=as.vector(ifelse(temp,"low","high"))
shengyi.data$score_24_2factors <- factor(score_24_2factors)


# 5.modeling
# training
f1 <- coxph(Surv(OS_month,OS_yn) ~ age_bin65 + AJCC_stage_bin, data = shengyi.data)
f2 <- coxph(Surv(OS_month,OS_yn) ~ age_bin65 + AJCC_stage_bin + score_24_2factors, data = shengyi.data)
predict_value1 <- predict(f1, shengyi.data)
predict_value1 <- 1/(1+exp(-predict_value1))
predict_value2 <- predict(f2, shengyi.data)
predict_value2 <- 1/(1+exp(-predict_value2))


#6. ROC
roc_clin <- timeROC(T=shengyi.data$OS_month,
                   delta=shengyi.data$OS_yn,
                   marker=predict_value1,
                   # weighting="cox",
                   cause=1,
                   times=seq(12,60,1),
                   iid = TRUE
                   )

roc_full <- timeROC(T=shengyi.data$OS_month,
                    delta=shengyi.data$OS_yn,
                    marker=predict_value2,
                    # weighting="cox",
                    cause=1,
                    times=seq(12,60,1),
                    iid = TRUE
                    )


#AUC and 95%CI
auc_clin <- round(roc_clin$AUC[c('t=60')],4)
auc_clin
confint(roc_clin,level = 0.95)$CI_AUC

auc_full <- round(roc_full$AUC[c('t=60')],4)
auc_full
confint(roc_full,level = 0.95)$CI_AUC


# ref_model VS full_model
# Plot the ROC curve
# png(filename="D:/projects/R_project/R_lung_TILs/OS_ROCAUC_figs/ROC-gdph_m60.png",width = 400,height = 280)
pdf("D:/projects/R_project/R_lung_TILs/OS_ROCAUC_figs/ROC-gdph_m60.pdf", width = 6, height = 4.5, family = "Times", onefile = FALSE)
# plot(dt, xlab="Risk scores for OS", ylab="Risk scores for DFS", col="blue", cex.axis=1.5,cex.lab=1.5,cex.main=2)
plot(roc_clin,time=60,col="darkcyan",lwd=2,title = FALSE)
plot(roc_full,time=60,col="tomato",lwd=2,title = FALSE, add=TRUE)
title(main=list("ROC(D1)",cex=1.25,col="black",font=1))
legend("bottomright", c("Reference model","Full model"), col=c("darkcyan","tomato"),lty=1,lwd=2,bty ="n")
dev.off()


# Plot the AUC curve
# png(filename="D:/projects/R_project/R_lung_TILs/OS_ROCAUC_figs/AUC-gdph_m60.png",width = 400,height = 280)
pdf("D:/projects/R_project/R_lung_TILs/OS_ROCAUC_figs/AUC-gdph_m60.pdf", width = 6, height = 4.5, family = "Times", onefile = FALSE)
plotAUCcurve(roc_clin,col = "darkcyan")
plotAUCcurve(roc_full,col ="tomato",add=TRUE)
title(main=list("AUC(D1)",cex=1.25,col="black",font=1))
legend("topright",c("Reference model","Full model"), col=c("darkcyan","tomato"),lty=1,lwd=2,bty="n")
dev.off()


