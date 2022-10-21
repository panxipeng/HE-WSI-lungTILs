##########################################
# Survival analysis    shengyi
# TTLs projects 

##########################################

library(openxlsx)
library("survival")
library("survminer")
library(gbm)


# 0. load the dataset
shengyi.data<-read.xlsx('shengyi_twice.xlsx')
head(shengyi.data)


# 1. binary through 1/4,median,3/4 cut-off value
# median cut-off value
temp = shengyi.data$inflam_per_tumor > median(na.omit(shengyi.data$inflam_per_tumor))
inflam_per_tumor_24=as.vector(ifelse(temp,"high","low"))
shengyi.data$inflam_per_tumor_24 <- inflam_per_tumor_24

temp = shengyi.data$inflam_per_stromal > median(na.omit(shengyi.data$inflam_per_stromal))
inflam_per_stromal_24=as.vector(ifelse(temp,"high","low"))
shengyi.data$inflam_per_stromal_24 <- inflam_per_stromal_24


# stage1+stage2=0, stage3=1
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



# 4. multivariate Cox regression analysis
# 1)TNM stage
f24_stage_sy <- coxph(Surv(OS_month, OS_yn) ~ AJCC_stage_bin, data = shengyi.data)
coxsummary=summary(f24_stage_sy) 
output = cbind(HR=coxsummary$coefficients[,"exp(coef)"],
               wald.test=coxsummary$coefficients[,"z"],
               pvalue=coxsummary$coefficients[,"Pr(>|z|)"],
               concordance=coxsummary$concordance)
write.csv(output, file = "D:/projects/R_project/R_lung_TILs/OS_results_table/shengyi_stage24_model.csv")

#lp
lp <- predict(object=f24_stage_sy, newdata =shengyi.data, type="lp")
#baseline hazard
basehaz.gbm(shengyi.data$OS_month, shengyi.data$OS_yn, lp, t.eval = 60, smooth = TRUE, cumulative = TRUE)


