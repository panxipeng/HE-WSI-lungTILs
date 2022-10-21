##########################################
# Univariable and multivariable Cox regression analyses for OS in four cohorts
# TTLs projects 

##########################################
library("survival")
library("survminer")



# 0. load the dataset
shengyi.data<-read.xlsx('shengyi_twice.xlsx')
head(shengyi.data)



# 1. binary through 1/4,median,3/4 cut-off value
# 1/4 cut-off value
temp = shengyi.data$inflam_per_tumor > quantile(na.omit(shengyi.data$inflam_per_tumor),probs = 0.25)
inflam_per_tumor_14=as.vector(ifelse(temp,"high","low"))
shengyi.data$inflam_per_tumor_14 <- inflam_per_tumor_14

temp = shengyi.data$inflam_per_stromal > quantile(na.omit(shengyi.data$inflam_per_stromal),probs = 0.25)
inflam_per_stromal_14=as.vector(ifelse(temp,"high","low"))
shengyi.data$inflam_per_stromal_14 <- inflam_per_stromal_14


# median cut-off value
temp = shengyi.data$inflam_per_tumor > median(na.omit(shengyi.data$inflam_per_tumor))
inflam_per_tumor_24=as.vector(ifelse(temp,"high","low"))
shengyi.data$inflam_per_tumor_24 <- inflam_per_tumor_24

temp = shengyi.data$inflam_per_stromal > median(na.omit(shengyi.data$inflam_per_stromal))
inflam_per_stromal_24=as.vector(ifelse(temp,"high","low"))
shengyi.data$inflam_per_stromal_24 <- inflam_per_stromal_24


# 3/4 cut-off value
temp = shengyi.data$inflam_per_tumor > quantile(na.omit(shengyi.data$inflam_per_tumor),probs = 0.75)
inflam_per_tumor_34=as.vector(ifelse(temp,"high","low"))
shengyi.data$inflam_per_tumor_34 <- inflam_per_tumor_34

temp = shengyi.data$inflam_per_stromal > quantile(na.omit(shengyi.data$inflam_per_stromal),probs = 0.75)
inflam_per_stromal_34=as.vector(ifelse(temp,"high","low"))
shengyi.data$inflam_per_stromal_34 <- inflam_per_stromal_34



# 2. Determine the optimal cutpoint of variables via Overall Survival
res.cut <- surv_cutpoint(shengyi.data, time = "OS_month", event = "OS_yn",
                         variables = c("inflam_per_tumor", "inflam_per_stromal"))
summary(res.cut)
cp = summary(res.cut)$cutpoint

# Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)

# Plug the binary data into the shengyi.data 
shengyi.data$inflam_per_tumor_cutbin <- res.cat$inflam_per_tumor
shengyi.data$inflam_per_stromal_cutbin <- res.cat$inflam_per_stromal


# stage1+stage2=0, stage3=1
tmp1 = shengyi.data$AJCC_stage
tmp1 = ifelse((tmp1 == 3),1,0)
shengyi.data$AJCC_stage_bin <- factor(tmp1)


# tumor site
tmp1 = shengyi.data$tumor_site
tmp1 = ifelse((tmp1 == "2" | tmp1 == "5"), 0, ifelse((tmp1 == "1" | tmp1 == "3" | tmp1 == "4"), 1, 2))
shengyi.data$tumor_site_bin <- factor(tmp1)



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



# 4. Single factor analysis via overall survival with clinical factors
out=data.frame()
for (i in 3:13)
{
  res.cox <- coxph(Surv(OS_month, OS_yn) ~ shengyi.data[,i], data = shengyi.data)
  m=data.frame(ShowRegTable(res.cox))
  colnames(m)[1]="95%CI"
  coxsummary=summary(res.cox)
  out=rbind (out,cbind(factor=colnames(shengyi.data)[i],
                       HR=coxsummary$coefficients[,"exp(coef)"],
                       wald.test=coxsummary$coefficients[,"z"],
                       CI =m[1],
                       pvalue=coxsummary$coefficients[,"Pr(>|z|)"]))
}
head(out)


# 5. Single factor analysis via overall survival with our factors.
out1=data.frame()
for (i in 32:length(shengyi.data))
{
  res.cox <- coxph(Surv(OS_month, OS_yn) ~ shengyi.data[,i], data = shengyi.data)
  m=data.frame(ShowRegTable(res.cox))
  colnames(m)[1]="95%CI"
  coxsummary=summary(res.cox)
  out1=rbind (out1,cbind(factor=colnames(shengyi.data)[i],
                         HR=coxsummary$coefficients[,"exp(coef)"],
                         wald.test=coxsummary$coefficients[,"z"],
                         CI =m[1],
                         pvalue=coxsummary$coefficients[,"Pr(>|z|)"]))
}
head(out1)
out2=rbind (out,out1)
write.csv(out2, file = "D:/projects/R_project/R_lung_TILs/OS_results_table/shengyi_single_factor_analysis.csv")



# 6. multivariate Cox regression analysis
# only inflam_per_stromal_24
f24_stromal_sy <- coxph(Surv(OS_month, OS_yn) ~ inflam_per_stromal_24, data = shengyi.data)
coxsummary=summary(f24_stromal_sy) 
output = cbind(HR=coxsummary$coefficients[,"exp(coef)"],
               wald.test=coxsummary$coefficients[,"z"],
               pvalue=coxsummary$coefficients[,"Pr(>|z|)"],
               concordance=coxsummary$concordance)
write.csv(output, file = "D:/projects/R_project/R_lung_TILs/OS_results_table/shengyi_stromal24_model.csv")

HRandCI=coxsummary$conf.int
print("the HR is:")
print(as.numeric(HRandCI[1]))
print(as.numeric(HRandCI[3]))
print(as.numeric(HRandCI[4]))

# Plot the K-M curve and get the p-value
# png(filename="D:/projects/R_project/R_lung_TILs/OS_result_figs/km-gdph_inflamperstromal24.png",width = 600,height = 420)
pdf("D:/projects/R_project/R_lung_TILs/OS_result_figs/km-gdph_inflamperstromal24.pdf", width = 8, height = 6, family = "Times", onefile = FALSE)
fit12 <- survfit(Surv(OS_month, OS_yn) ~ inflam_per_stromal_24, data = shengyi.data) 
ggsurvplot(fit12,
           conf.int = FALSE,    #If TRUE, plots confidence interval.
           # linetype = "strata",        #Change line type by groups
           # surv.median.line = "hv",    #character vector for drawing a horizontal/vertical line at median survival.
           palette = "lancet",
           xlab = "Time(months)", # x
           ylab = "OS probability", # y
           break.x.by = 20, font.tickslab = 18, font.x =18, font.y =18,
           
           legend=c(0.12,0.3), font.legend = 16,
           legend.title = "Groups",
           legend.labs = c("Low risk", "High risk"),
           
           risk.table = TRUE,  fontsize = 6.5, tables.theme = theme_cleantable(),
           
           # pval = TRUE,        #If logical and TRUE, the p-value is added on the plot.
           pval.size = 6.5, pval.coord = c(0.1,0.05), 
           pval.method = TRUE, pval.method.coord = c(0.1,0.25), #the method that calculate the p-value.
           pval = paste0("p=0.012, ","HR=",round(HRandCI[1],2),"(",round(HRandCI[3],2),"-",round(HRandCI[4],2),")")
)
dev.off()



# only inflam_per_tumor_24
f24_tumor_sy <- coxph(Surv(OS_month, OS_yn) ~ inflam_per_tumor_24, data = shengyi.data)
coxsummary=summary(f24_tumor_sy)
output = cbind(HR=coxsummary$coefficients[,"exp(coef)"],
               wald.test=coxsummary$coefficients[,"z"],
               pvalue=coxsummary$coefficients[,"Pr(>|z|)"],
               concordance=coxsummary$concordance)
write.csv(output, file = "D:/projects/R_project/R_lung_TILs/OS_results_table/shengyi_tumor24_model.csv")

HRandCI=coxsummary$conf.int
print("the HR is:")
print(as.numeric(HRandCI[1]))
print(as.numeric(HRandCI[3]))
print(as.numeric(HRandCI[4]))

# Plot the K-M curve and get the p-value
# png(filename="D:/projects/R_project/R_lung_TILs/OS_result_figs/km-gdph_inflampertumor24.png",width = 600,height = 420)
pdf("D:/projects/R_project/R_lung_TILs/OS_result_figs/km-gdph_inflampertumor24.pdf", width = 8, height = 6, family = "Times", onefile = FALSE)
fit12 <- survfit(Surv(OS_month, OS_yn) ~ inflam_per_tumor_24, data = shengyi.data) 
ggsurvplot(fit12,
           conf.int = FALSE,    #If TRUE, plots confidence interval.
           # linetype = "strata",        #Change line type by groups
           # surv.median.line = "hv",    #character vector for drawing a horizontal/vertical line at median survival.
           palette = "lancet",
           xlab = "Time(months)", # x
           ylab = "OS probability", # y
           break.x.by = 20, font.tickslab = 18, font.x =18, font.y =18,
           
           legend=c(0.12,0.3), font.legend = 16,
           legend.title = "Groups",
           legend.labs = c("Low risk", "High risk"),
           
           risk.table = TRUE,  fontsize = 6.5, tables.theme = theme_cleantable(),
           
           # pval = TRUE,        #If logical and TRUE, the p-value is added on the plot.
           pval.size = 6.5, pval.coord = c(0.1,0.05), 
           pval.method = TRUE, pval.method.coord = c(0.1,0.25), #the method that calculate the p-value.
           # pval = "p=0.00014, HR=2.78(1.60-4.81)"                #need input manually.
           pval = paste0("p<0.001, ","HR=",round(HRandCI[1],2),"(",round(HRandCI[3],2),"-",round(HRandCI[4],2),")")
)
dev.off()



# only our features(oof),
tmp1 = shengyi.data$inflam_per_tumor
shengyi.data$inflam_per_tumor_log <- log(tmp1 + 1)
tmp1 = shengyi.data$inflam_per_stromal
shengyi.data$inflam_per_stromal_log <- log(tmp1 + 1)

f24_ooftwo_sy <- coxph(Surv(OS_month, OS_yn) ~ inflam_per_tumor_log + inflam_per_stromal_log, data = shengyi.data)
coxsummary=summary(f24_ooftwo_sy) 
output = cbind(HR=coxsummary$coefficients[,"exp(coef)"],
               wald.test=coxsummary$coefficients[,"z"],
               pvalue=coxsummary$coefficients[,"Pr(>|z|)"],
               concordance=coxsummary$concordance)
write.csv(output, file = "D:/projects/R_project/R_lung_TILs/OS_results_table/shengyi_ourfeatures24_model.csv")

coef6 = as.numeric(coxsummary$coefficients[,"coef"][1])
coef7 = as.numeric(coxsummary$coefficients[,"coef"][2])
score_shengyi_2factors = coef6*shengyi.data$inflam_per_tumor_log + coef7*shengyi.data$inflam_per_stromal_log

# median cut-off value
temp = score_shengyi_2factors > median(na.omit(score_shengyi_2factors))
# score_24_2factors=as.vector(ifelse(temp,"high","low"))
score_24_2factors=as.vector(ifelse(temp,"low","high"))
shengyi.data$score_24_2factors <- score_24_2factors

ftmp <- coxph(Surv(OS_month, OS_yn) ~ score_24_2factors, data = shengyi.data)
coxsummary=summary(ftmp)
HRandCI=coxsummary$conf.int
print("the HR is:")
print(as.numeric(HRandCI[1]))
print(as.numeric(HRandCI[3]))
print(as.numeric(HRandCI[4]))

# Plot the K-M curve and get the p-value
# png(filename="D:/projects/R_project/R_lung_TILs/OS_result_figs/km-gdph_2factors.png",width = 600,height = 420)
pdf("D:/projects/R_project/R_lung_TILs/OS_result_figs/km-gdph_2factors.pdf", width = 8, height = 6, family = "Times", onefile = FALSE)
fit12 <- survfit(Surv(OS_month, OS_yn) ~ score_24_2factors, data = shengyi.data) 
ggsurvplot(fit12,
           conf.int = FALSE,    #If TRUE, plots confidence interval.
           # linetype = "strata",        #Change line type by groups
           # surv.median.line = "hv",    #character vector for drawing a horizontal/vertical line at median survival.
           palette = "lancet",
           xlab = "Time(months)", # x
           ylab = "OS probability", # y
           break.x.by = 20, font.tickslab = 18, font.x =18, font.y =18,
           
           legend=c(0.12,0.3), font.legend = 16,
           legend.title = "Groups",
           legend.labs = c("Low risk", "High risk"),
           
           risk.table = TRUE,  fontsize = 6.5, tables.theme = theme_cleantable(),
           
           # pval = TRUE,        #If logical and TRUE, the p-value is added on the plot.
           pval.size = 6.5, pval.coord = c(0.1,0.05), 
           pval.method = TRUE, pval.method.coord = c(0.1,0.25), #the method that calculate the p-value.
           # pval = "p=2e-04, HR=2.68(1.56-4.62)"                #need input manually.
           pval = paste0("p<0.001, ","HR=",round(HRandCI[1],2),"(",round(HRandCI[3],2),"-",round(HRandCI[4],2),")")
)
dev.off()



# clinical and 2 factors combinedV2  
# add at 17/5/2022
f24_stepcombinedV2_sy <- coxph(Surv(OS_month, OS_yn) ~ age_bin65 + gender + smoking + AJCC_stage_bin
             + score_24_2factors, data = shengyi.data)
coxsummary=summary(f24_stepcombinedV2_sy)  

m=data.frame(ShowRegTable(f24_stepcombinedV2_sy))
output = cbind(HRinter=m$exp.coef...confint.,
               wald.test=coxsummary$coefficients[,"z"],
               pvalue=coxsummary$coefficients[,"Pr(>|z|)"],
               concordance=coxsummary$concordance)
write.csv(output, file = "D:/projects/R_project/R_lung_TILs/OS_results_table/shengyi_24_modelV2.csv")

# computing the score value and plot the K-M curve
coef21 = as.numeric(coxsummary$coefficients[,"coef"][1])
coef22 = as.numeric(coxsummary$coefficients[,"coef"][2])
coef23 = as.numeric(coxsummary$coefficients[,"coef"][3])
coef24 = as.numeric(coxsummary$coefficients[,"coef"][4])
coef25 = as.numeric(coxsummary$coefficients[,"coef"][5])

temp = as.numeric(shengyi.data$age_bin65) > 1
age_bin65vv=as.vector(ifelse(temp,1,0))

temp = as.numeric(shengyi.data$gender) > 1
gendervv=as.vector(ifelse(temp,1,0))

temp = as.numeric(shengyi.data$smoking) > 1
smokingvv=as.vector(ifelse(temp,1,0))

temp = as.numeric(shengyi.data$AJCC_stage_bin) > 1
AJCC_stage_binvv=as.vector(ifelse(temp,1,0))

temp = shengyi.data$score_24_2factors == "low"
score_24_2factorsvv=as.vector(ifelse(temp,1,0))

score_shengyiV2 = coef21*age_bin65vv + coef22*gendervv + coef23*smokingvv + coef24*AJCC_stage_binvv + coef25*score_24_2factorsvv

temp = score_shengyiV2 > median(na.omit(score_shengyiV2))
# score_24=as.vector(ifelse(temp,"high","low"))
score_24=as.vector(ifelse(temp,"low","high"))
shengyi.data$score_24V2 <- score_24

# 
ftmp <- coxph(Surv(OS_month, OS_yn) ~ score_24V2, data = shengyi.data)
coxsummary=summary(ftmp)
HRandCI=coxsummary$conf.int
print("the HR is:")
print(as.numeric(HRandCI[1]))
print(as.numeric(HRandCI[3]))
print(as.numeric(HRandCI[4]))

# Plot the K-M curve and get the p-value
# png(filename="D:/projects/R_project/R_lung_TILs/OS_result_figs/km-gdph_clinical_and_2factorsV2.png",width = 600,height = 420)
pdf("D:/projects/R_project/R_lung_TILs/OS_result_figs/km-gdph_clinical_and_2factorsV2.pdf", width = 8, height = 6, family = "Times", onefile = FALSE)
fit12 <- survfit(Surv(OS_month, OS_yn) ~ score_24V2, data = shengyi.data) 
ggsurvplot(fit12,
           conf.int = FALSE,    #If TRUE, plots confidence interval.
           # linetype = "strata",        #Change line type by groups
           # surv.median.line = "hv",    #character vector for drawing a horizontal/vertical line at median survival.
           palette = "lancet",
           xlab = "Time(months)", # x
           ylab = "OS probability", # y
           break.x.by = 20, font.tickslab = 18, font.x =18, font.y =18,
           
           legend=c(0.12,0.3), font.legend = 16,
           legend.title = "Groups",
           legend.labs = c("Low risk", "High risk"),
           
           risk.table = TRUE,  fontsize = 6.5, tables.theme = theme_cleantable(),
           
           # pval = TRUE,        #If logical and TRUE, the p-value is added on the plot.
           pval.size = 6.5, pval.coord = c(0.1,0.05), 
           pval.method = TRUE, pval.method.coord = c(0.1,0.25), #the method that calculate the p-value.
           # pval = "p<0.0001, HR=9.11(4.62-17.95)"                #need input manually.
           pval = paste0("p<0.001, ","HR=",round(HRandCI[1],2),"(",round(HRandCI[3],2),"-",round(HRandCI[4],2),")")
)
dev.off()


