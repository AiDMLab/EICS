setwd("D:\\guoman\\cibersort\\改进\\lasso\\20")
#a <- read.table("CIBERSORT-lasso-Results-9.txt",header=T,sep="\t",row.names=1,check.names=F)
#a <- read.table("D:\\guoman\\cibersort\\data\\CIBERSORT-Results.txt",header=T,sep="\t",row.names=1,check.names=F)
#b <- read.csv("D:\\guoman\\cibersort\\改进\\lasso\\100\\比例.csv",header=T,sep=",",row.names=1,check.names=F)
b <- read.table("D:\\guoman\\cibersort\\cibersort-改进\\原cibersort\\20个\\PBMCs-Fig3a-Flow-Cytometry.txt",header=T,sep="\t",row.names=1,check.names=F)
#a <- read.csv("D:\\guoman\\cibersort\\改进\\神经网络\\keras\\20\\result_yuan.csv",header=T,sep=",",row.names=1,check.names=T)
#a <- read.csv("D:\\guoman\\cibersort\\改进\\岭回归\\20\\result_yuan.csv",header=T,sep=",",row.names=1,check.names=T)
a <- read.table("D:\\guoman\\cibersort\\表达矩阵\\CIBERSORT-lasso-Results.txt",header=T,sep="\t",row.names=1,check.names=F)
b <- read.table("D:\\guoman\\cibersort\\EPIC_SVR\\data\\PBMCs-Fig3a-Flow.txt",header=T,sep="\t",row.names=1,check.names=F)
b <- read.csv("D:\\guoman\\cibersort\\EPIC_SVR\\作图\\mix100\\比例2.csv",header=T,sep=",",row.names=1,check.names=F)
a <- read.table("D:\\guoman\\cibersort\\EPIC_SVR\\作图\\mix100\\EPIC-SVR-65-mix-Results-1.txt",header=T,sep="\t",row.names=1,check.names=F)
b<-read.csv("D:\\guoman\\cibersort\\EPIC_SVR\\20实际数据\\HIC-10.csv",header=T,sep=",",row.names=1,check.names=F)


a <- data.matrix(a) #将数据框转换为数字矩阵
b <- data.matrix(b)

###################取出公共集#######################
agns <- row.names(a) #行名
bgns <- row.names(b)
binta <- bgns %in% bgns
b <- b[binta,]
aintb <- agns %in% row.names(b)
a <- a[aintb,]

######## a 标准化（化为百分之...）##################
a_stand <- NULL
num0 <- dim(a)[1] ##a的行数
tmpr <- rownames(a)  #a行名
for(i in 1:num0)
{
  s <- a[i,]
  if(sum(s)==0)
  {a_stand<-rbind(a_stand,s)}
  else
    {
      s <- s / sum(s)
      a_stand <- rbind(a_stand,s)
    }
}
rownames(a_stand) <- tmpr 

b_stand <- NULL
num0 <- dim(b)[1] ##a的行数
tmpr <- rownames(b)  #a行名
for(i in 1:num0)
{
  s <- b[i,]
  if(sum(s)==0)
  {b_stand<-rbind(b_stand,s)}
  else
  {
    s <- s / sum(s)
    b_stand <- rbind(b_stand,s)
  }
}
rownames(b_stand) <- tmpr 
###########初始化##############################
num1 <- dim(a_stand)[2]###########样本a的列数
num2 <- dim(b_stand)[2]###########样本b的列数
result <- NULL
i <- 1
rmse <- rep(0,num1) ###初始化。重复函数,第一个参数重复的内容，第二个参数重复的次数
corrv <- rep(0,num1)
corrv1 <- rep(0,num1)
###################循环求均方根误差和相关系数
while(i <= num1){
  a_standi <- a_stand[,i]
  tmpr <- colnames(a_stand)[i]#######取列名，colnames(b)b的全部列名，rownames(b)[i]第i个b的列名
  #############在b中寻找和a[i,]相同id的样本
  for(j in 1:num2){
    if(tmpr==colnames(b_stand)[j])
    {
      k <- j 
    }
  }
  bk=b_stand[,k]
  rmse[i] <- sqrt((mean((a_standi - bk)^2)))
  corrv[i] <- cor(a_standi,bk,method="spearman")
  corrv1[i] <- cor(a_standi,bk,method="pearson")
  print(paste("id:",tmpr,", RMSE:", rmse[i],", SPEARMAN:", corrv[i],", PEARSON:",corrv1[i], sep=""))
  result <- rbind(result, c(tmpr, rmse[i], corrv[i],corrv1[i]))
  i <- i+1
}
colnames(result) <- c("id", "RMSE", "SPEARMAN","PEARSON")
result_df <- as.data.frame(result)
result_df
write.csv(result_df, paste("./","相关系数.csv",sep=""))

###############作散点图
library(ggplot2)
file<-read.csv("相关系数.csv")
ggplot(data = file,aes(x=X,y=PEARSON))+geom_point(colour="blue")+ylim(0,1)+ labs(x="Sample",y="Pearson Correlation Coefficient")


cor.test(a_standi,bk,method="spearman")
cor.test(a_standi,bk,method="pearson")

