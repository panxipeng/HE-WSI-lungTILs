##########################################
# Survival analysis   
# TTLs projects 

##########################################

library(openxlsx)
library("survival")
library("survminer")



# 0. load the dataset
shengyi.data<-read.xlsx('shengyi_twice.xlsx')



# 1. factors
# stage1+stage2=0, stage3=1
tmp1 = shengyi.data$AJCC_stage
tmp1 = ifelse((tmp1 == 3),1,0)
shengyi.data$AJCC_stage_bin <- factor(tmp1)

shengyi.data$age_bin65 <- factor(shengyi.data$age_bin65)



# 2. multivariate Cox regression analysis
# 1)only stage
f24_stage_sy <- coxph(Surv(OS_month, OS_yn) ~ AJCC_stage_bin, data = shengyi.data)
ci <- summary(f24_stage_sy)$concordance[1] 
se <- summary(f24_stage_sy)$concordance[2] 
cilow <- ci - 1.96 * se #95%CI
ciup <- ci + 1.96 * se
sy_stage_CI <- paste0(round(ci,3)," (",round(cilow,3),"-",round(ciup,3),")")
sy_stage_AIC <- round(AIC(f24_stage_sy),1)  


# 2)only our features(oof), WELL
tmp1 = shengyi.data$inflam_per_tumor
shengyi.data$inflam_per_tumor_log <- log(tmp1 + 1)
tmp1 = shengyi.data$inflam_per_stromal
shengyi.data$inflam_per_stromal_log <- log(tmp1 + 1)

f24_ooftwo_sy <- coxph(Surv(OS_month, OS_yn) ~ inflam_per_tumor_log + inflam_per_stromal_log, data = shengyi.data)
coxsummary=summary(f24_ooftwo_sy)
coef6 = as.numeric(coxsummary$coefficients[,"coef"][1])
coef7 = as.numeric(coxsummary$coefficients[,"coef"][2])
score_shengyi_2factors = coef6*shengyi.data$inflam_per_tumor_log + coef7*shengyi.data$inflam_per_stromal_log

# median cut-off value
temp = score_shengyi_2factors > median(na.omit(score_shengyi_2factors))
# score_24_2factors=as.vector(ifelse(temp,"high","low"))
score_24_2factors=as.vector(ifelse(temp,"low","high"))
shengyi.data$score_24_2factors <- score_24_2factors

f24_ooftwo_sy <- coxph(Surv(OS_month, OS_yn) ~ score_24_2factors, data = shengyi.data)
ci <- summary(f24_ooftwo_sy)$concordance[1] 
se <- summary(f24_ooftwo_sy)$concordance[2] 
cilow <- ci - 1.96 * se #95%CI
ciup <- ci + 1.96 * se
sy_score_CI <- paste0(round(ci,3)," (",round(cilow,3),"-",round(ciup,3),")")
sy_score_AIC <- round(AIC(f24_ooftwo_sy),1)  


# 3)stage and WELL
f24_stageand2factors_sy <- coxph(Surv(OS_month, OS_yn) ~ AJCC_stage_bin + score_24_2factors, data = shengyi.data)
ci <- summary(f24_stageand2factors_sy)$concordance[1] 
se <- summary(f24_stageand2factors_sy)$concordance[2] 
cilow <- ci - 1.96 * se #95%CI
ciup <- ci + 1.96 * se
sy_stageandscore_CI <- paste0(round(ci,3)," (",round(cilow,3),"-",round(ciup,3),")")
sy_stageandscore_AIC <- round(AIC(f24_stageand2factors_sy),1)  


# 4)clinical factors
fclinical_step_sy <- coxph(Surv(OS_month, OS_yn) ~ age_bin65 + AJCC_stage_bin, data = shengyi.data)
ci <- summary(fclinical_step_sy)$concordance[1] 
se <- summary(fclinical_step_sy)$concordance[2] 
cilow <- ci - 1.96 * se #95%CI
ciup <- ci + 1.96 * se
sy_clinic_CI <- paste0(round(ci,3)," (",round(cilow,3),"-",round(ciup,3),")")
sy_clinic_AIC <- round(AIC(fclinical_step_sy),1)  


# 5)full model
f24_stepcombinedV2_sy <- coxph(Surv(OS_month, OS_yn) ~ age_bin65 + AJCC_stage_bin + score_24_2factors, data = shengyi.data)
ci <- summary(f24_stepcombinedV2_sy)$concordance[1] 
se <- summary(f24_stepcombinedV2_sy)$concordance[2] 
cilow <- ci - 1.96 * se #95%CI
ciup <- ci + 1.96 * se
sy_full_CI <- paste0(round(ci,3)," (",round(cilow,3),"-",round(ciup,3),")")
sy_full_AIC <- round(AIC(f24_stepcombinedV2_sy),1) 

output = cbind(sy_stage_CI=sy_stage_CI, sy_stage_AIC=sy_stage_AIC,
               sy_score_CI=sy_score_CI, sy_score_AIC=sy_score_AIC,
               sy_stageandscore_CI=sy_stageandscore_CI, sy_stageandscore_AIC=sy_stageandscore_AIC,
               sy_clinic_CI=sy_clinic_CI, sy_clinic_AIC=sy_clinic_AIC,
               sy_full_CI=sy_full_CI, sy_full_AIC=sy_full_AIC)
write.csv(output, file = "D:/projects/R_project/R_lung_TILs/OS_results_table/shengyi_CIandAIC.csv")




