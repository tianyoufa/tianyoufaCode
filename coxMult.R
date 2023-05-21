rm(list = ls())
setwd("E:/cancerPrognosis/MyLUAD/code/processing/coxMult/")

#加载要用到的包
library("survival")

# 导入文件 此时的文件展示为行名表示样本名，列名表示生存信息与miRNA信息，共463行，307列
cli_exp <- read.table('limma_clinical_express.txt',sep="\t",header=T,check.names=F,row.names = 1)

# 定义combines函数，整合P值小于0.1的数据
combines <- function(coef, Int){
  # 定义存储经过多变量cox分析后被筛选出来的分析数据的数据框
  mulgene <- data.frame('Characteristics'='', 'Hazard Ratio'='', 'CI95'='','P Value'='')[-1,]
  # 获取结束数据说表示的基因名
  rownames <- rownames(coef)
  # 循环遍历筛选P值小于0.1的基因数据
  for(i in seq(from = 1,to = length(rownames),by = 1)){
    if(coef[i,5]<=0.1){
      #将被筛选的数据整合成一个一行三列的数据框
      CI <- paste0(round(Int[i, 3:4], 2), collapse = "-")
      mulcox <- data.frame('Characteristics' = rownames[i],
                           'Hazard Ratio' = round(coef[i,2],3),
                           'CI95' = CI,
                           'P Value' = round(coef[i,5],2))
      # 将数据框进行拼接
      mulgene <- rbind.data.frame(mulgene,mulcox) 
    }
  }
  #返回数据框
  return (mulgene)
}

# 定义获取组合变量的函数
get_names <- function(number){
  j <- 1
  k <- 1
  # 定义连接符
  op <- '+'
  # 定义一个列表，存储多变量分析的基因组合名
  names <- list()
  # 获取基因id的字符串拼接，用于多变量分析
  name <- character()
  for(i in seq.default(from = 9, to = length(cli_exp), by = 1)){
    if(j == 1){
      name <- paste(name, colnames(cli_exp)[i])
      j <- j + 1
    }
    else{
      if(j != number){
        name <- paste(name, op, colnames(cli_exp)[i])
        j <- j + 1
      }
      else{
        if(j == number){
          name <- paste(name, op, colnames(cli_exp)[i])
          names[k] <- sub(' ', '', name)
          k <- k + 1
          name <- character()
          j <- 1
        }
      }
    }
  }
  return(names)
}

# 获取多变量分析筛选得到的基因
get_genes <- function(names){
  # 定义存储经过单变量cox分析后被筛选出来的分析数据的数据框
  genes <- data.frame('Characteristics'='', 'Hazard Ratio'='', 'CI95'='','P Value'='')[-1,]
  for(i in seq.default(from = 1, to = length(names), by = 1)){
    # 多变量分析
    mul <- as.formula(paste0('Surv(time, status)~',names[i]))
    # 进行多变量cox分析
    fit.quan <- coxph(mul,data = cli_exp)
    # 开始进行分析
    fit1 <- summary(fit.quan)
    # 获取c_index指数
    fit1$concordance
    # 获取分析的coefficients结果信息
    fit.coef.all <- fit1$coefficients
    # 获取分析的conf.int结果信息
    fit.int.all <- fit1$conf.int
    # 筛选符合条件的基因
    mulgene <- combines(fit.coef.all, fit.int.all)
    # print(dim(mulgene)[1])
    if(dim(mulgene)[1] != 0){
      # 将数据框进行合并
      genes <- rbind.data.frame(genes, mulgene)
    }
  }
  return(genes)
}

#封装定义回归分析函数
UniCox <- function(name){
  FML <- as.formula(paste0('Surv(time, status) ~ ', name))
  # 单变量cox回归分析
  res.cox <- coxph(FML, data = cli_exp)
  # 获取分析结果
  GSum <- summary(res.cox)
  #HR风险比
  HR <- round(GSum$coefficients[, 2], 2)
  #P值
  PValue <- round(GSum$coefficients[, 5], 3)
  #95%的置信区间，上下界
  CI <- paste0(round(GSum$conf.int[, 3:4], 2), collapse = "-")
  #将上述三个整合成一个一行三列的数据框
  unicox <- data.frame('Characteristics' = name,
                       'Hazard Ratio' = HR,
                       'CI95' = CI,
                       'P Value' = PValue)
  #返回数据框
  return (unicox)
}


#获取数据的列名(数据里包含的所有MiRNA的名称) 共647个
VarNames <- colnames(cli_exp)[9:length(cli_exp)]
length(VarNames)
#应用lapply()函数，缓解内存加快运算速度
univar_analysis <- lapply(VarNames, UniCox)
# 定义存储经过单变量cox分析后被筛选出来的分析数据的数据框
genes1 <- data.frame('Characteristics'='', 'Hazard Ratio'='', 'CI95'='','P Value'='')[-1,]
# 对每一个基因进行分析
for(i in VarNames){
  # 调用单变量回归分析函数
  a <- UniCox(i)
  if(a$P.Value < 0.1){
    # 将数据框进行合并
    genes1 <- rbind.data.frame(genes1,a)
  } 
}
# 观察有多少基因被选择 有308个
dim(genes1)
genes1

# 3个基因组合多变量分析
name3 <- get_names(3)
genes3 <- get_genes(name3)
# 观察有多少基因被选择 有214个
dim(genes3)
genes3

# 5个基因组合多变量分析
name5 <- get_names(5)
genes5 <- get_genes(name5)
# 观察有多少基因被选择 有174个
dim(genes5)
genes5

# 7个基因组合多变量分析
name7 <- get_names(7)
genes7 <- get_genes(name7)
# 观察有多少基因被选择 有138个
dim(genes7)
genes7

# 9个基因组合多变量分析
name9 <- get_names(9)
genes9 <- get_genes(name9)
# 观察有多少基因被选择 有136个
dim(genes9)
# 存储多变量分析结果
write.table(genes9,file="COXgenes9.txt",sep="\t",row.names=F,quote=F)
genes9


# 构建话韦恩图的数据
char1 <- c(genes1$Characteristics)
length(char1)

char3 <- c(genes3$Characteristics)
num <- length(char1) - length(char3)
gene1 <- rep(c('gene1'), each = num)
char3 <- c(char3,gene1)
length(char3)

char5 <- c(genes5$Characteristics)
num <- length(char1) - length(char5)
gene5 <- rep(c('gene5'), each = num)
char5 <- c(char5,gene5)
length(char5)

char7 <- c(genes7$Characteristics)
num <- length(char1) - length(char7)
gene7 <- rep(c('gene7'), each = num)
char7 <- c(char7,gene7)
length(char7)

char9 <- c(genes9$Characteristics)
num <- length(char1) - length(char9)
gene9 <- rep(c('gene9'), each = num)
char9 <- c(char9,gene9)
length(char9)

# 将数据转换为数据框
Characteristics <- data.frame('CoxPH1' = char1, 'CoxPH3' = char3, 
                              'CoxPH5' = char5,'CoxPH7' = char7,
                              'CoxPH9' = char9)

# 存储画韦恩图的数据
write.table(Characteristics,file="VennData.txt",sep="\t",row.names=F,quote=F)
























