library(caret)
library(survival)
library(survminer)
library(survMisc)
library(glmnet)
library(timeROC)

#ICC
library(irr)
library(magrittr)

df_reader1 <- read.csv("test_reader1.csv")
df_reader2 <- read.csv("test_reader2.csv")

can_ta <- rep(1:2163)

for (i in 1:2163){
  a <- icc(data.frame(df_reader1[i],df_reader2[i]), model = "oneway")
  can_ta[i] <- a$value
}

write.csv(can_ta, "icc.csv")

#===========Spearman correlation analysis of voxels and radiomics features
library(psych)
r_mat <- corr.test(df_tra[,4],df_tra[5:1725], method = 'Spearman')
r_mat_r <- r_mat$r
r_mat_p <- r_mat$p
r <- data.frame(r = t(r_mat_r), p.adjust = t(r_mat_p), row.names = colnames(r_mat_r))

df_voxel_3D <- df_tra[,c(5:694)][,-which(r_mat_p < 0.05 & abs(r_mat_r) > 0.5)]
df_voxel_PTV <- df_tra[,c(1096:1725)][,-which(r_mat_p < 0.05 & abs(r_mat_r) > 0.5)]

#=========univariate cox analysis
p_value <- rep(1, 1454)
for(j in 3:1454){
  dfInter2 <- df_voxel[c(1,2,j)]
  colnames(dfInter2) <- c("time","event", "radio")  # survival analysis
  fit2 <- coxph(Surv(time,event) ~ radio, data = dfInter2)
  
  p_value[j] <- summary(fit2)$coefficients[5]
}
names(p_value) <- colnames(df_voxel)
index <- (1:1454)[p_value <= 0.05]

df_uni_tra <- df_voxel[c(index)]
df_uni_tra <- cbind(df_tra[,c(2,3)],df_uni_tra)


#========varclus procedure

#==============lasso
y <- Surv(df_varclus$time, df_varclus$status)

x <- as.matrix(df_varclus[,-c(1,2)])

fit<-glmnet(x, y, family="cox")

plot(fit,xvar = "lambda", label = TRUE)

cv_fit <- cv.glmnet(x, y, type.measure = "deviance", family="cox", nfolds = 5)

plot(cv_fit)

cv_fit$lambda.1se

coefficients<-coef(cv_fit,s=cv_fit$lambda.1se) 
Active.Index<-which(coefficients!=0)     #coefficient was not zero
Active.coefficients<-coefficients[Active.Index]   
print(Active.coefficients)
print(Active.Index)


feature_lasso <- rownames(coefficients)[Active.Index]
df_lasso <- df_cor_tra[,feature_lasso[1:length(feature_lasso)]]
feature_num_lasso  <- dim(df_lasso)[2]


score <- df_lasso
for(i in 1: feature_num_lasso){ 
  score[i] <- Active.coefficients[i] * df_lasso[i]
}
score <- apply(score, 1, sum) # calculate the rad-scores
df_lasso$score <- score 
df_lasso <- data.frame(df_varclus[,c(1,2)],df_lasso)
fit <- coxph(Surv(time, status) ~ score, data = df_lasso) # rad-score could classify the patients into two groups according risk
c1 <- cutp(fit, defCont=3)
c1[]
class <- as.numeric(score >= 0.4343483)
df_lasso$class <- class
fit_s <- survfit(Surv(time, status) ~ class, data = df_lasso)
ggsurvplot(fit_s, data = df_lasso, size = 1,
           palette = c("#3193A9", "#FF9B37"),
           conf.int = T, legend = c(0.8, 0.2),
           pval = T,
           risk.table = T, tables.height = 0.3,
           legend.labs = c("Low risk", "High risk"),
           xlab = "Time(months)",xlim = c(0,72), break.x.by = 12, axes.offset = FALSE)

#======================survival comparison between GGO and low risk
df_survival_tra <- rbind(df_lasso[c(1,2,26)],df_GGO_1[2:4])

fit_s <- survfit(Surv(time, status) ~ class, data = df_survival_tra)

ggsurvplot(fit_s, data = df_survival_tra, size = 1,
           palette = c("#3193A9", "#FF9B37","#7AAB43"),
           conf.int = T, legend = c(0.8, 0.2),
           pval = T,
           risk.table = T, tables.height = 0.3,
           legend.labs = c("Low risk", "High risk","GGO"),
           xlab = "Time(months)",xlim = c(0,72), break.x.by = 12, axes.offset = FALSE)


#============validation 1
#===============================validation 1
df_val1_lasso <- df_val1[,feature_lasso[1:length(feature_lasso)]]
feature_num_lasso_val1  <- dim(df_val1_lasso)[2]
df_val1_lasso <- scale(df_val1_lasso)

score_val1 <- df_val1_lasso
for(i in 1: feature_num_lasso_val1){ 
  score_val1[i] <- Active.coefficients[i] * df_val1_lasso[i]
}
score_val1 <- apply(score_val1, 1, sum) # calculate the rad-scores
df_val1_lasso <- as.data.frame(df_val1_lasso)
df_val1_lasso$score_val1 <- score_val1
df_val1_lasso <- data.frame(df_val1[,c(2:3)],df_val1_lasso)

class_val1 <- as.numeric(score_val1 >= 0.4343483)
df_val1_lasso$class_val1 <- class_val1
fit_v <- survfit(Surv(time, status) ~ class_val1, data = df_val1_lasso)
ggsurvplot(fit_v, data = df_val1_lasso, size = 1,
           palette = c("#3193A9", "#FF9B37"),
           conf.int = T, legend = c(0.8, 0.2),
           pval = T,
           risk.table = T, tables.height = 0.3,
           legend.labs = c("Low risk", "High risk"),
           xlab = "Time(months)",xlim = c(0,72), break.x.by = 12, axes.offset = FALSE)

#============validation 2
df_val2_lasso <- df_val2[,feature_lasso[1:length(feature_lasso)]]
feature_num_lasso_val2  <- dim(df_val2_lasso)[2]
df_val2_lasso <- data.frame(time, status, df_val2_lasso) 
score_val2 <- df_val2_lasso
for(i in 1: feature_num_lasso_val2){ 
  score_val2[i] <- Active.coefficients[i] * df_val2_lasso[i]
}
score_val2 <- apply(score_val2, 1, sum) # calculate the rad-scores
df_val2_lasso$score_val2 <- score_val2
class_val2 <- as.numeric(score_val2 >= 0.4343483)
df_val2_lasso$class_val2 <- class_val2
fit_v <- survfit(Surv(time, status) ~ class_val2, data = df_val1_lasso)
ggsurvplot(fit_v, data = df_val2_lasso, size = 1,
           palette = c("#3193A9", "#FF9B37"),
           conf.int = T, legend = c(0.8, 0.2),
           pval = T,
           risk.table = T, tables.height = 0.3,
           legend.labs = c("Low risk", "High risk"),
           xlab = "Time(months)",xlim = c(0,72), break.x.by = 12, axes.offset = FALSE)

#=================Time ROC
library(timeROC)

time_roc_res <- timeROC(T = df_lasso$time, delta = df_lasso$status, 
                        marker = df_lasso$score, cause = 1, weighting="marginal",
                        times = c(3 * 12, 5 * 12),ROC = TRUE,iid = TRUE)

time_roc_res$AUC
confint(time_roc_res, level = 0.95)$CI_AUC

par(xaxt = 'n')
plotAUCcurve(time_roc_3D, add = F, conf.int = F, col = "#0072B5FF")
par(xaxt = 'l')
axis(side = 1,labels = c(12, 24, 36, 48, 60, 72), at = c(12, 24, 36, 48, 60, 72), las = 1, cex.axis = 1, tick = -0.02)

#===================C-index
library(survcomp)
f1 <- cph(Surv(time, status) ~ score, data = df_lasso)
cindex_tra <- concordance.index(predict(f1),
                                surv.time = df_lasso$time,
                                surv.event = df_lasso$status,method = "noether")
print(c(cindex_tra$c.index, cindex_tra$lower, cindex_tra$upper))



#===================Calibration curves
units(df_nom_val1$RFStime) <- "mouth"

f1 <- cph(Surv(df_nom_val1$RFStime,df_nom_val1$RFS) ~ age + overall_stage_clinical + score, data = df_nom_val1, x=T,y=T,surv=T,time.inc=36)
f2 <- cph(Surv(df_nom_val1$RFStime,df_nom_val1$RFS) ~ age + overall_stage_clinical + score, data = df_nom_val1, x=T,y=T,surv=T,time.inc=60)

cal1 <- calibrate(f1, cmethod='KM', method="boot", u=36, m=54, B=163)
cal2 <- calibrate(f2, cmethod='KM', method="boot", u=60, m=54, B=163)

## 绘制校正曲线，??rms::calibrate查看详细参数说明
par(mar=c(8,5,3,2),cex = 1.0)
plot(cal1,lwd=2,lty=1,
     errbar.col= "#3193A9",
     xlim=c(0.4,0.98),ylim=c(0.4,0.98),
     xlab="Predicted Probability of 3-Year RFS",
     ylab="Actual 3-Year RFS (proportion)",
     col="#3193A9")
par(new = T)
par(mar=c(8,5,3,2),cex = 1.0)
plot(cal2,lwd=2,lty=1,
     errbar.col= "#FF9B37",
     xlim=c(0.4,0.98),ylim=c(0.4,0.98),
     col= "#FF9B37")

#===DCA
library(rmda)

df_nom_val1$sur_binary_5 <- as.numeric(df_nom_val1$RFStime >= 60)
bb <-  decision_curve(sur_binary_5 ~ age + overall_stage_clinical, confidence.intervals = 0.95,
                      data = df_nom_val1)
plot_decision_curve(bb, confidence.intervals = F, col="#3193A9")
par(new = T)
bb1 <-  decision_curve(sur_binary_5 ~ score,
                       data = df_nom_val1)
plot_decision_curve(bb1, confidence.intervals = F, col="#FF9B37")
par(new = T)
bb2 <-  decision_curve(sur_binary_5 ~ age + overall_stage_clinical + score,
                       data = df_nom_val1)
plot_decision_curve(bb2, confidence.intervals = F, col="red")

#==================================计算NRI
library(nricens)

z.std <- as.matrix(df_nom_val2[,c( "age","overall_stage_clinical")])
z.new <- as.matrix(df_nom_val2[,c("score","age", "overall_stage_clinical")])

mstd = coxph(Surv(RFStime,RFS) ~  age  + overall_stage_clinical, data = df_nom_val2, x=TRUE)
mnew= coxph(Surv(RFStime,RFS) ~  age  + overall_stage_clinical + score, data = df_nom_val2, x=TRUE)

#predicted risk
p.std= get.risk.coxph(mstd, t0=60)
p.new= get.risk.coxph(mnew, t0=60)

nricens(time = df_nom_val2$RFStime, event= df_nom_val2$RFS, 
        mdl.std= mstd, mdl.new = mnew, 
        t0 = 60, cut = c(0.4, 0.6),
        niter = 1000, updown = 'category')

#
#==================================IDI
library(survIDINRI)
t0 = 60
indata0=as.matrix(subset(df_nom_val2,select=c(RFStime, RFS, age, overall_stage_clinical)))
indata1=as.matrix(subset(df_nom_val2,select=c(RFStime, RFS, age, overall_stage_clinical, score)))

covs0 <- as.matrix(indata0[,c(-1,-2)])
covs1<-as.matrix(indata1[,c(-1,-2)])

x<- IDI.INF(indata0[,1:2],covs0, covs1, t0, npert=1000)

IDI.INF.OUT(x)
