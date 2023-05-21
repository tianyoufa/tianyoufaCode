rm(list = ls())
setwd('E:/cancerPrognosis/MyLUAD/code/processing/')

# 安装CRAN来源的R包 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 安装包
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("factoextra")
BiocManager::install("ggplot2")
BiocManager::install("ggrepel")
BiocManager::install("EnhancedVolcano")

# 导入包
library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(factoextra)
library(ggrepel)
library(EnhancedVolcano)
library(ggpubr)

# 参数定义
normalNum <- 45 # ormal组样品数目
tumorNum <- 450 # tumor组样品数目
fdrFilter <- 0.05 # fdr临界值
logFCfilter <- 2 # logFC临界值

# 导入文件 此时的文件展示为行名表示基因名，列名表示样本信息名，共28585行，495列
exp <- read.table('normal_tumor.txt',sep="\t",header=T,check.names=F,row.names = 1)
# 样本注释 建立一个列表，共495个数据，其中前45个表示为normal，余下的表示450个为tumor
group <- c(rep("normal", normalNum), rep("tumor",tumorNum)) %>% factor(., levels = c("normal", "tumor"), ordered = F)
# 将group转换为一个带有顺序因子矩阵
group <- model.matrix(~factor(group)+0) 
# 给group矩阵赋列名，第一列为normal，第二列为tumor
colnames(group) <- c("normal", "tumor")
# 差异分析
# 设计对比矩阵
exp.matrix <- makeContrasts(normal - tumor, levels = group)
# 拟合模型
# lmFit()通过为每一个基因拟合线性模型来估计fold changes和standard errors。
exp.fit <- lmFit(exp, group)
# contrasts.fit()将原始模型重新定向到对比模型，得到新的系数和标准误差。
fit <- contrasts.fit(exp.fit, exp.matrix)
# eBayes()使用trend=TRUE对标准误差进行经验贝叶斯平滑，计算每个对比中每个基因的moderated t-statistic和log-odds。
fit <- eBayes(fit, trend = TRUE)
# topTable()给出一个最有可能在给定对比下差异表达的基因列表。
tempOutput <- topTable(fit, n = Inf, adjust = "fdr")
# 去掉数据中有NA的行或列
tempOutput <- na.omit(tempOutput)
# 更改列名
colnames(tempOutput)[5] = 'FDR'

# 跳转文件存储空间
setwd('E:/cancerPrognosis/MyLUAD/code/processing/limma/')
# 存储全部的差异基因分析
write.table(tempOutput, "tempOutput.txt", col.names = colnames(tempOutput), row.names = rownames(tempOutput), quote = F, sep = "\t")
# 筛选出差异基因并存储 共647个
diff_RNA <- tempOutput[(tempOutput$FDR < fdrFilter & (tempOutput$logFC > logFCfilter | tempOutput$logFC < (-logFCfilter))),]
dim(diff_RNA)
write.table(diff_RNA, "diff_RNA.txt", col.names = colnames(diff_RNA), row.names = rownames(diff_RNA), quote = F, sep = "\t")
# 获取并存储上调基因 共348个
diffup = diff_RNA[(diff_RNA$FDR < fdrFilter & diff_RNA$logFC<(-logFCfilter)),]
dim(diffup)
write.table(diffup, "diffup.txt", col.names = colnames(diffup), row.names = rownames(diffup), quote = F, sep = "\t")
# 获取并存储下调基因 共299个
diffdown = diff_RNA[(diff_RNA$FDR < fdrFilter & diff_RNA$logFC>logFCfilter),]
dim(diffdown)
write.table(diffdown, "diffdown.txt", col.names = colnames(diffdown), row.names = rownames(diffdown), quote = F, sep = "\t")

dim(tempOutput)


# 画火山图
# 读取画图需要的数据 FDR与logFC
volcano <- data.frame(gene = rownames(tempOutput),FDR = tempOutput$FDR,logFC = tempOutput$logFC,stringsAsFactors = F)
head(volcano)
# 添加标记信息 
volcano$type <- ifelse(volcano$FDR < 0.05 & volcano$logFC >logFCfilter,'down',
                       ifelse(volcano$FDR < 0.05 & volcano$logFC < (-logFCfilter) ,'up','stable')) 
head(volcano)
#设置标签——将p值小于0.03且差异倍数大于3的进行标注
volcano$label<-ifelse(volcano$FDR<0.03 & abs(volcano$logFC)>=3.5,"Y","N")
volcano$label<-ifelse(volcano$label == 'Y', as.character(volcano$gene), '') 
# 统计信息
sum(volcano$type == 'down')
sum(volcano$type == 'up')
sum(volcano$type == 'stable')
# 存储火山图数据
write.table(volcano, "volcano.txt", col.names = colnames(volcano), row.names = rownames(volcano), quote = F, sep = "\t")

# 开始绘图
graph <- ggplot(volcano,aes(x=logFC, y=-1*log10(FDR), label = rownames(volcano)))+ # 加载数据，定义横纵坐标
  geom_point(aes(color = type),alpha=0.5,size = 2.5)+ # 绘制散点图，分组依据是数据框的type列
  scale_color_manual(values=c("#0072B6","#BD0097","#BC3C28"))+ # 自定义颜色，将values更改成你想要的三个颜色
  geom_vline(xintercept=c(-logFCfilter,logFCfilter), colour="#0D2E05", lty = 4,lwd=0.6)+ # 在图上添加虚线
  geom_hline(yintercept = -log10(0.05),colour="#0D2E05", lty = 4,lwd=0.6)+ # 在图上添加虚线
  labs(x="logFC",y="-log10(FDR)", title = "volcanoplot")+ # 定义标题，x轴，y轴名称
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "top")+
  geom_text_repel(aes(x = logFC,
                      y = -log10(FDR),
                      label=label),
                  max.overlaps = 10000,
                  size=3,
                  box.padding=unit(0.8,'lines'),
                  point.padding=unit(0.8, 'lines'),
                  segment.color='#00BD19',
                  show.legend=FALSE) 
  
graph 


ggsave(plot = graph, filename = 'volcanoplot.pdf', width = 9,height = 9)





