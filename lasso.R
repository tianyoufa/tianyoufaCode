rm(list = ls())
setwd("E:/cancerPrognosis/MyLUAD/code/processing/lasso/")

library("glmnet")
library("survival")
library('stats')
library('ggplot2')
library('ggpubr')
library('survminer')

# 读取文件
cli_exp <- read.table("xgb_gbm_rxpress.txt",header=T,sep="\t",row.names=1,check.names=F)     #读取文件

# 获取模型需要的下x,y值
x <- as.matrix(cli_exp[,c(9:ncol(cli_exp))])
y <- data.matrix(Surv(cli_exp$time, cli_exp$status))
# 模型拟合
fit <- glmnet(x, y, family = "cox", maxit = 1000,alpha = 1)
# 开始画图
plot(fit,xvar="lambda",label=TRUE)
# 交叉验证
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
plot(cvfit)

abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")

# 变量筛选
# 获取拟合的系数
coef <- coef(fit, s = cvfit$lambda.min)
# 获取变量索引
index <- which(coef != 0)
# 根据索引查找拟合值
actCoef <- coef[index]
# 根据值查找基因名
lassoGene <- row.names(coef)[index]
kmeans_data <- cli_exp[,lassoGene]
colnames(kmeans_data)
# 整合临床信息 列名
lassoGene <- c('time', 'status', 'stage', lassoGene)
lassoExp <- cli_exp[,lassoGene]
colnames(lassoExp)
# 为数据添加标签

# 设置随机种子
set.seed(1)
kmeans_3 <- kmeans(kmeans_data,center = 3)
#分类贴标签
label <- kmeans_3$cluster
lassoExp$label <- label


# 存储数据
write.table(lassoExp, "lassoExp.txt", col.names = colnames(lassoExp), row.names = rownames(lassoExp), quote = F, sep = "\t")

rownames(lassoExp)



# 设置随机种子
set.seed(5)
kmeans_4 <- kmeans(kmeans_data,center = 4)
#分类贴标签
label <- kmeans_4$cluster
lassoExp$label <- label


# 画分组的图
colnames(lassoExp)[21] = 'group'
lassoExp$time <- lassoExp$time / 365
k_label <- survfit(Surv(time, status) ~ group, data = lassoExp)
ggsurvplot(k_label, 
           pval=TRUE, 
           data =lassoExp,
           font.main = 2,
           font.x = 12,
           font.y = 12,
           font.tickslab = 16,
           font.tickslab.size =4,
           size=1,
           censor.size=2,
           axis.line=36,
           xlab="Time (year)",
           ylab="Survival rate",
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
)












