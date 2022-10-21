
##########################################
# clinicopathologic features

##########################################

library(Rmisc)
library(openxlsx)
library(survival)

# 0. load the dataset
shengyi.data<-read.xlsx('shengyi_twice.xlsx')


############## age #############
# https://blog.csdn.net/Anne999/article/details/65627685
shapiro.test(shengyi.data$age)


############################ 
summary(shengyi.data$age)



###########################
tmp <- shengyi.data$smoking
all <- length(tmp)
s1 <- length(tmp[tmp==0])
s2 <- length(tmp[tmp==1])
R1 <- s1/all
R2 <- s2/all



############## OS month #############
# https://wirereed.com/t/f43f7ecdd277765caa414df68a81fd10
rkm_sy=survfit(Surv(OS_month, OS_yn==0)~1,data=shengyi.data)
rkm_sy


############## DFS month #############
rkm_sy=survfit(Surv(DFS_month, DFS_yn==0)~1,data=shengyi.data)
rkm_sy



# Pearson's Chi-squared test
age <- c(181,95,97,42,98,17,136,127)
age1 <- matrix(age,nrow = 2,ncol = 4)
age1
chisq.test(age1,correct = T)


gender <- c(141,135,86,53,58,57,116,147)
gender1 <- matrix(gender,nrow = 2,ncol = 4)
gender1
chisq.test(gender1)


smoke <- c(208,68,65,74,74,41,37,226)
smoke1 <- matrix(smoke, nrow = 2,ncol = 4)
smoke1
chisq.test(smoke1)

