library(caret)
library(survival)
library(survminer)
library(survMisc)
library(glmnet)
library(randomForest)
library(timeROC)

#ICC
library(irr)
library(magrittr)

df_reader1 <- read.csv("test_reader1.csv")
df_reader2 <- read.csv("test_reader2.csv")

can_ta <- rep(1:2071)

for (i in 1:2071){
  a <- icc(data.frame(df_reader1[i],df_reader2[i]), model = "oneway")
  can_ta[i] <- a$value
}

write.csv(can_ta, "icc3.csv")

#=========univariate cox analysis
df <- read.csv("Radiomics_2D&3D&PTV_total.csv")
df_tra <- df[which(df$cohorts==1),]
df_val1 <- df[which(df$cohorts==2),]
df_val2 <- df[which(df$cohorts==3),]
df_GGO_1 <- read.csv("GGO_tra_val1.csv")
df_GGO_2 <- read.csv("GGO_val2.csv")

dfBinary_tra <- df_tra[c(3:2638)]

p_value <- rep(1, 2636)

for(j in 3:2636){
  dfInter2 <- dfBinary_tra[c(1,2,j)]
  colnames(dfInter2) <- c("time","event", "radio")  # survival analysis
  fit2 <- coxph(Surv(time,event) ~ radio, data = dfInter2)
  
  p_value[j] <- summary(fit2)$coefficients[5]
}

index <- (1:2636)[p_value <= 0.05]

df_uni_tra <- dfBinary_tra[c(index)]

#===========correlation
corMatrix_GTV =  cor(df_uni_tra, y = NULL, use = "ev",method = "pearson")
highly_correlated_columns_GTV = findCorrelation(
  corMatrix_GTV,
  cutoff = 0.8,
  verbose = FALSE,
  names = FALSE,
  exact = TRUE
)
df_cor_tra <- df_uni_tra[, -highly_correlated_columns_GTV]


library(corrplot)
corMatrix_cor =  cor(df_cor_tra, y = NULL, use = "ev",method = "pearson")

pdf(file = "pearson correlation.pdf", w = 7, h = 7 )
corrplot(corMatrix_cor, method = "color", tl.pos = 'n')
dev.off()

#===========LASSO -----To select the optimal feature set for signature development
#=============Normalization
scale_func <- function(x){
  return (scale(x, center=T, scale=T)) 
}
data_scale <- apply(df_cor_tra,2,scale_func)

#==============lasso
y <- Surv(df_tra$time, df_tra$status)

x <- as.matrix(data_scale)

fit<-glmnet(x, y, family="cox")

plot(fit,xvar = "lambda", label = TRUE)

cv_fit <- cv.glmnet(x, y, type.measure = "deviance", family="cox", nfolds = 10)

plot(cv_fit)

# coefficients of slected features in the lasso
# coef(cv_fit)
cv_fit$lambda
cv_fit$lambda.min
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
time <- df_tra$time
status <- df_tra$status
df_lasso <- data.frame(time,status,df_lasso)

fit_s <- survfit(Surv(time, status) ~ class, data = df_lasso)

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

df_val1_lasso <- df_val1[,feature_lasso[1:length(feature_lasso)]]
feature_num_lasso_val1  <- dim(df_val1_lasso)[2]

time <- df_val1$time
status <- df_val1$status

df_val1_lasso <- data.frame(time, status, df_val1_lasso) 

score_val1 <- df_val1_lasso
for(i in 1: 22){ 
  score_val1[i] <- Active.coefficients[i] * df_val1_lasso[i]
}

score_val1 <- apply(score_val1, 1, sum) # calculate the rad-scores

df_val1_lasso$score_val1 <- score_val1


#======================survival comparison between GGO and low risk
colnames(df_GGO_1)[colnames(df_GGO_1) == 'class'] <- 'class_val1'

df_survival_val1 <- rbind(df_val1_lasso[c(1,2,26)],df_GGO_1[2:4])

fit_s <- survfit(Surv(time, status) ~ class_val1, data = df_survival_val1)

ggsurvplot(fit_s, data = df_survival_val1, size = 1,
           palette = c("#3193A9", "#FF9B37","#7AAB43"),
           conf.int = T, legend = c(0.8, 0.2),
           pval = T,
           risk.table = T, tables.height = 0.3,
           legend.labs = c("Low risk", "High risk","GGO"),
           xlab = "Time(months)",xlim = c(0,72), break.x.by = 12, axes.offset = FALSE)

#============validation 2
df_val2_lasso <- df_val2[,feature_lasso[1:length(feature_lasso)]]

feature_num_lasso_val2  <- dim(df_val2_lasso)[2]

df_val2_lasso <- data.frame(time, status, df_val2_lasso) 

score_val2 <- df_val2_lasso
for(i in 1: 22){ 
  score_val2[i] <- Active.coefficients[i] * df_val2_lasso[i]
}

score_val2 <- apply(score_val2, 1, sum) # calculate the rad-scores

df_val2_lasso$score_val2 <- score_val2

#======================survival comparison between GGO and low risk
colnames(df_GGO_2)[colnames(df_GGO_2) == 'class'] <- 'class_val1'

df_survival_val2 <- rbind(df_val2_lasso[c(1,2,26)],df_GGO_2[2:4])

fit_s <- survfit(Surv(time, status) ~ class_val2, data = df_survival_val2)

ggsurvplot(fit_s, data = df_survival_val2, size = 1,
           palette = c("#3193A9", "#FF9B37","#7AAB43"),
           conf.int = T, legend = c(0.8, 0.2),
           pval = T,
           risk.table = T, tables.height = 0.3,
           legend.labs = c("Low risk", "High risk","GGO"),
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

units(df_nom_tra$RFStime) <- "mouth"

f1 <- cph(Surv(df_nom_tra$RFStime,df_nom_tra$RFS) ~ age_c + overall_stage_clinical + score, data = df_nom_tra, x=T,y=T,surv=T,time.inc=36)
f2 <- cph(Surv(df_nom_tra$RFStime,df_nom_tra$RFS) ~ age_c + overall_stage_clinical + score, data = df_nom_tra, x=T,y=T,surv=T,time.inc=60)

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
p.std= get.risk.coxph(mstd, t0=36)
p.new= get.risk.coxph(mnew, t0=36)

nricens(time = df_nom_val2$RFStime, event= df_nom_val2$RFS, 
        mdl.std= mstd, mdl.new = mnew, 
        t0 = 60, cut = c(0.4, 0.6),
        niter = 1000, updown = 'category')

#
#==================================IDI
library(survIDINRI)
t0= 60
indata0=as.matrix(subset(df_nom_val2,select=c(RFStime, RFS, age, overall_stage_clinical)))
indata1=as.matrix(subset(df_nom_val2,select=c(RFStime, RFS, age, overall_stage_clinical, score)))

covs0 <- as.matrix(indata0[,c(-1,-2)])
covs1<-as.matrix(indata1[,c(-1,-2)])

x<- IDI.INF(indata0[,1:2],covs0, covs1, t0, npert=1000)

IDI.INF.OUT(x)

IDI.INF.GRAPH(x)   



