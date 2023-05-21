rm(list = ls())
setwd("E:/cancerPrognosis/MyLUAD/code/processing/survival/")

library(survival)
library(survminer)

# 生存曲线绘制函数
Graph <- function(i) {
  lasso<-lasso[order(lasso[,i]),]#升序排列
  risk<-rep('Low',length(lasso[,1]))#先令所有等级为低等级
  lasso$risk<-risk
  lasso$risk[207:414]<-'High'
  fit <- survfit(Surv(time, status) ~ risk, data = lasso)
  ggsurvplot(fit, 
             pval=TRUE, 
             data =lasso,
             font.main = 2,
             font.x = 12,
             font.y = 12,
             font.tickslab = 16,
             font.tickslab.size =4,
             size=1,
             censor.size=2,
             axis.line=36,
             legend.title = i,
             xlab="Time (year)",
             ylab="Survival rate",
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             surv.median.line = "hv", # Specify median survival
             ggtheme = theme_bw(), # Change ggplot2 theme
             # legend.labs = 
             #   c("low risk", "high risk"),    # change legend labels.
             palette = c("#FF6103", "#3D9140")
  )
}


#读取文件
lasso <- read.table("lassoExp.txt",header=T,sep="\t",row.names=1,check.names=F)
upExp <- read.table("upExp.txt",header=T,sep="\t",row.names=1,check.names=F)
downExp <- read.table("downExp.txt",header=T,sep="\t",row.names=1,check.names=F)
lasso$time <- lasso$time / 365
upExp$time <- upExp$time / 365
downExp$time <- downExp$time / 365

# 画分组的图
colnames(lasso)[21] = 'group'
k_label <- survfit(Surv(time, status) ~ group, data = lasso)
ggsurvplot(k_label, 
           pval=TRUE, 
           data =lasso,
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


# 绘制筛选出的17个基因的生存曲线图
setwd("E:/cancerPrognosis/MyLUAD/code/processing/survival/graph/")

upGraph <- function(i) {
  upExp<-upExp[order(upExp[,i]),]#升序排列
  risk<-rep('low',length(upExp[,1]))#先令所有等级为低等级
  upExp$risk<-risk
  upExp$risk[207:414]<-'high'
  fit <- survfit(Surv(time, status) ~ risk, data = upExp)
  ggsurvplot(fit, 
             pval=TRUE, 
             data =upExp,
             font.main = 2,
             font.x = 12,
             font.y = 12,
             font.tickslab = 16,
             font.tickslab.size =4,
             size=1,
             censor.size=2,
             axis.line=36,
             legend.title = i,
             xlab="Time (year)",
             ylab="Survival rate",
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             surv.median.line = "hv", # Specify median survival
             ggtheme = theme_bw(), # Change ggplot2 theme
             # legend.labs = 
             #   c("low risk", "high risk"),    # change legend labels.
             palette = c("#FF6103", "#3D9140")
  )
}

downGraph <- function(i) {
  downExp<-downExp[order(downExp[,i]),]#升序排列
  risk<-rep('high',length(downExp[,1]))#先令所有等级为低等级
  downExp$risk<-risk
  downExp$risk[207:414]<-'low'
  fit <- survfit(Surv(time, status) ~ risk, data = downExp)
  ggsurvplot(fit, 
             pval=TRUE, 
             data =downExp,
             font.main = 2,
             font.x = 12,
             font.y = 12,
             font.tickslab = 16,
             font.tickslab.size =4,
             size=1,
             censor.size=2,
             axis.line=36,
             legend.title = i,
             xlab="Time (year)",
             ylab="Survival rate",
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             surv.median.line = "hv", # Specify median survival
             ggtheme = theme_bw(), # Change ggplot2 theme
             # legend.labs = 
             #   c("low risk", "high risk"),    # change legend labels.
             palette = c("#FF6103", "#3D9140")
  )
}

upName <- colnames(upExp)[4:14]

for(i in upName){
  g <- upGraph(i)
  print(g)
  ggsave(paste(i,".png"), plot = print(g),units = "in")
}

downName <- colnames(downExp)[4:9]

for(i in downName){
  g <- downGraph(i)
  print(g)
  ggsave(paste(i,".png"), plot = print(g),units = "in")
}
















k_gene <- survfit(Surv(time, status) ~ AGER, data = lasso)
ggsurvplot(k_gene, 
           pval=TRUE, 
           data =lasso,
           font.main = 2,
           font.x = 12,
           font.y = 12,
           font.tickslab = 16,
           font.tickslab.size =4,
           size=1,
           censor.size=2,
           axis.line=36,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           
)










