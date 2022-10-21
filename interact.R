##########################################
# Interaction    
# TTLs projects 

##########################################
library(openxlsx)
library("survival")
library("survminer")


# 0. read the data
data <- read.xlsx('fourcenters.xlsx')
rownames(data) <- data[,2]


# 1. median cut-off value
temp = data$inflam_per_tumor > median(na.omit(data$inflam_per_tumor))
inflam_per_tumor_24=as.vector(ifelse(temp,"high","low"))
data$inflam_per_tumor_24 <- inflam_per_tumor_24

temp = data$inflam_per_stromal > median(na.omit(data$inflam_per_stromal))
inflam_per_stromal_24=as.vector(ifelse(temp,"high","low"))
data$inflam_per_stromal_24 <- inflam_per_stromal_24

# computation by cox model with inflam_per_tumor_log and inflam_per_stromal_log
coef6 = -0.440005
coef7 = -0.418071
tmp1 = data$inflam_per_tumor
data$inflam_per_tumor_log <- log(tmp1 + 1)
tmp1 = data$inflam_per_stromal
data$inflam_per_stromal_log <- log(tmp1 + 1)
score_tmp = coef6*data$inflam_per_tumor_log + coef7*data$inflam_per_stromal_log
temp = score_tmp > median(na.omit(score_tmp))
# score_24_2factors=as.vector(ifelse(temp,"high","low"))
score_24_2factors=as.vector(ifelse(temp,"low","high"))
data$score_24_2factors <- factor(score_24_2factors)


# 2. subgroup
# K-M curves
ID <- data$pathID[data$AJCC_stage_detail == "IB" ]
df_1b <- data[ID,]

ID <- data$pathID[data$AJCC_stage_detail == "IIA" ]
df_2a <- data[ID,]

ID <- data$pathID[data$AJCC_stage_detail == "IB" | data$AJCC_stage_detail == "IIA"]
df_1b2a <- data[ID,]


ID <- data$pathID[data$AJCC_stage == "1" | data$AJCC_stage == "2"]
df_1and2 <- data[ID,]

ID <- data$pathID[data$AJCC_stage == 1]
df_stage_1 <- data[ID,]
ID <- data$pathID[data$AJCC_stage == 2]
df_stage_2 <- data[ID,]
ID <- data$pathID[data$AJCC_stage == 3]
df_stage_3 <- data[ID,]



f1<- coxph(Surv(OS_month, OS_yn) ~ adjuvant_chemotherapy + score_24_2factors + adjuvant_chemotherapy*score_24_2factors, data = df_1b)
summary(f1)

f2<- coxph(Surv(OS_month, OS_yn) ~ adjuvant_chemotherapy + score_24_2factors + adjuvant_chemotherapy*score_24_2factors, data = df_2a)
summary(f2)

f3<- coxph(Surv(OS_month, OS_yn) ~ adjuvant_chemotherapy + score_24_2factors + adjuvant_chemotherapy*score_24_2factors, data = df_1b2a)
summary(f3)

f4<- coxph(Surv(DFS_month, DFS_yn) ~ adjuvant_chemotherapy + score_24_2factors + adjuvant_chemotherapy*score_24_2factors, data = df_1b)
summary(f4)

f5<- coxph(Surv(DFS_month, DFS_yn) ~ adjuvant_chemotherapy + score_24_2factors + adjuvant_chemotherapy*score_24_2factors, data = df_2a)
summary(f5)

f6<- coxph(Surv(DFS_month, DFS_yn) ~ adjuvant_chemotherapy + score_24_2factors + adjuvant_chemotherapy*score_24_2factors, data = df_1b2a)
summary(f6)




