library(openxlsx)
library(actuar)
library(psych)
library(MASS)
library(mnormt)
library(mFilter)
library(tseries)
library(urca)
library(forecast)
library(tmvtnorm)
library(zoo)           
library(xts)            
library(plyr)
library(TSA)
library(WaveletComp)
library(gtools)

setwd("C:/Users/Administrator/Desktop/发表论文/2经济周期协同性-多方法综合评价")


##按指标读取

GDP=read.xlsx(xlsxFile = "GDP水平值.xlsx",sheet = 3,detectDates = TRUE, colNames = TRUE)

province_number<<-ncol(GDP)-1
Time<<-nrow(GDP)



#去除OLS趋势线
 GDP_deOLS=rep(list(0),province_number)
 GDP_ts=GDP
for(i in 2:(province_number+1))
{
  GDP_deOLS[[i-1]]=lm(GDP[,i]~GDP[,1])

  GDP_ts[,i]=GDP[,i]-GDP_deOLS[[i-1]]$fitted.values
}



#CF滤波
GDP_cf=GDP_ts
GDP_cf[,-1]=cffilter(GDP_ts[,-1],pl =6,pu =32, root = T,drift = TRUE)$cycle


#C-M指数
CM=rep(list(matrix(0,nrow=province_number,ncol=province_number)),Time)


for(i in 2:(province_number+1))
{
  for(j in 2:(province_number+1))
  {
    for(t in 1:(Time))
    {
      CM[[t]][i-1,j-1]=1-0.5*((GDP_cf[t,i]-mean(GDP_cf[,i]))/sd(GDP_cf[,i])-
                            (GDP_cf[t,j]-mean(GDP_cf[,j]))/sd(GDP_cf[,j]))^2
    }
  }
}


TCM=CM

for(i in 1:(province_number))
{
  for(j in 1:(province_number))
  {
    for(t in 1:(Time))
    {
      if(i!=j)
      {
        TCM[[t]][i,j]=0.5*log(1/(1-CM[[t]][i,j]))
      }
    }
  }
}

ACM=TCM
for(i in 1:(province_number))
{
  for(j in 1:(province_number))
  {
    for(t in 1:(Time))
    {
      if(i!=j)
      {
        ACM[[t]][i,j]=(exp(2*TCM[[t]][i,j])-1)/(exp(2*TCM[[t]][i,j])+1)
      }
    }
  }
}




trans_TCM=matrix(0,nrow=Time,ncol=4)

# TCM_test=matrix(0,nrow=Time,ncol=435)
# east_TCM_test=matrix(0,nrow=Time,ncol=55)
# mid_TCM_test=matrix(0,nrow=Time,ncol=28)
# west_TCM_test=matrix(0,nrow=Time,ncol=55)

#均值
for(t in 1:Time)
{
  trans_TCM[t,1]=mean(ACM[[t]][1:11,1:11][lower.tri(ACM[[t]][1:11,1:11])])
  trans_TCM[t,2]=mean(ACM[[t]][12:19,12:19][lower.tri(ACM[[t]][12:19,12:19])])
  trans_TCM[t,3]=mean(ACM[[t]][20:30,20:30][lower.tri(ACM[[t]][20:30,20:30])])
  trans_TCM[t,4]=mean(ACM[[t]][1:30,1:30][lower.tri(ACM[[t]][1:30,1:30])])
}


# matplot(trans_TCM,type='l')


#小波变换

epoch.seq = seq(from = as.POSIXct("2000-03-01"),
                 to = as.POSIXct("2018-08-01"), by = 3*30*24*3600)



#时间轴
ticks <- as.Date(c("2000-03-01","2001-09-01","2003-03-01","2004-09-01","2006-03-01",
             "2007-09-01","2009-03-01","2010-09-01","2012-03-01","2013-09-01",
             "2015-03-01","2016-09-01","2018-03-01"))

labels <- c("03-2000","09-2001","03-2003","09-2004","03-2006",
              "09-2007","03-2009","09-2010","03-2012","09-2013",
              "03-2015","09-2016","03-2018")


# wc.image(wavelet_coherence1, which.image = "wc", color.key = "quantile",
#          n.levels = 250,
#          siglvl.contour = 0.05, siglvl.arrow = 0.05,which.arrow.sig = "wc",
#          legend.params = list(lab = "小波动态相关系数",label.digits= 2),
#          periodlab = "季度频率",timelab="",
#          color.palette = "gray((n.levels):1/n.levels)",
#          col.ridge = "blue",label.time.axis = TRUE,show.date = TRUE,
#          date.format = "%Y-%m",
#          spec.time.axis = list(at = ticks, labels = labels, las = 1))


#江苏
# wavelet_jiangsu <- analyze.wavelet(GDP_cf,"江苏",
#                                         loess.span = 0, dt = 1, dj = 1/16,
#                                         # window.type.t = "ham", window.type.s = "ham",
#                                         # window.size.t = 3, window.size.s = 1,
#                                         make.pval = TRUE, n.sim = 100,lowerPeriod = 6,
#                                         upperPeriod = 32)


# wt.image(wavelet_jiangsu, color.key = "quantile",
#          n.levels = 250,
#          legend.params = list(lab = "江苏小波功率谱",label.digits= 2),
#          periodlab = "季度频率",timelab="",
#          color.palette = "gray((n.levels):1/n.levels)",
#          col.ridge = "blue",label.time.axis = TRUE,show.date = TRUE,
#          date.format = "%Y-%m",
#          spec.time.axis = list(at = ticks, labels = labels, las = 1))


compare_pair_list=matrix(0,nrow=75,ncol=435*2)

data_for_pair=GDP_cf[,-1]

#write.csv(x=compare_pair_list,file = "省区对.csv")

combine_number=combinations(n = 30,r = 2)

compare_pair_list=data_for_pair[,t(combine_number)]

wavelet_coherence=rep(list(0),435)

for(i in 1:435)
{
  wavelet_coherence[[i]]=analyze.coherency(compare_pair_list[,((i-1)*2+1):((i-1)*2+2)],
                                            my.pair = c(1,2),
                                            loess.span = 0, dt = 1, dj = 1/16,
                                            window.type.t = "ham", window.type.s = "ham",
                                            window.size.t = 3, window.size.s = 1,
                                            make.pval = TRUE, n.sim = 10,lowerPeriod = 6,
                                            upperPeriod = 32)
}


rou_wc=rep(list(0),435)



#实部
for(i in 1:435)
{
  rou_wc[[i]]=Re(wavelet_coherence[[i]]$Coherency)
  for(t in 1:75)
  {
    for(f in 1:39)
    {
      if((1-wavelet_coherence[[i]]$Coherence.pval[f,t])<0.9)
      {
        rou_wc[[i]][f,t]=0
      }
    }
  }
}


combination_names=combn(names(GDP_cf[,-1]),m = 2)



#均值

atanh_rou_twc=matrix(0,nrow = 75,ncol=435)

for(i in 1:435)
{
  for(t in 1:75)
  {
   for(j in 1:39)
   {
      atanh_rou_twc[t,i]=
        atanh_rou_twc[t,i]+
        log(1+rou_wc[[i]][j,t])-log(1-rou_wc[[i]][j,t])
      #rou_wc[[i]][j,t]

   }
  }
}

atanh_rou_twc=atanh_rou_twc/78
rou_twc=tanh(atanh_rou_twc)
#rou_twc=atanh_rou_twc


east_rou_twc=matrix(0,nrow=75,ncol=1)
mid_rou_twc=matrix(0,nrow=75,ncol=1)
west_rou_twc=matrix(0,nrow=75,ncol=1)
east=0
mid=0
west=0

for(t in 1:75)
{
  for(i in 1:435) 
  {
    if(combination_names[1,i]%in%c("北京","天津","辽宁","上海","江苏","浙江",
                                    "福建","山东","广东","海南")&
         combination_names[2,i]%in%c("北京","天津","辽宁","上海","江苏","浙江",
                                     "福建","山东","广东","海南"))
    {
      east_rou_twc[t,]=east_rou_twc[t,]+rou_twc[t,i]
      east=east+1
      
    }else if(combination_names[1,i]%in%c("山西","吉林","黑龙江","安徽","江西","河南",
                                         "湖北","湖南")&
             combination_names[2,i]%in%c("山西","吉林","黑龙江","安徽","江西","河南",
                                         "湖北","湖南"))
    {
      mid_rou_twc[t,]=mid_rou_twc[t,]+rou_twc[t,i]
      mid=mid+1
    }else if(combination_names[1,i]%in%c("内蒙古","四川","广西","重庆","贵州","云南",
                                         "陕西","甘肃","青海","宁夏","新疆")&
             combination_names[2,i]%in%c("内蒙古","四川","广西","重庆","贵州","云南",
                                         "陕西","甘肃","青海","宁夏","新疆"))
    {
      west_rou_twc[t,]=west_rou_twc[t,]+rou_twc[t,i]
      west=west+1
    }
  }
}



trans_rou_twc=cbind(east_rou_twc/(east/75),mid_rou_twc/(mid/75),
                    west_rou_twc/(west/75),rowMeans(rou_twc))

plot(rowMeans(rou_twc),type='l')
lines(east_rou_twc/(east/75),col='red')
lines(mid_rou_twc/(mid/75),col='blue')
lines(west_rou_twc/(west/75),col='green')




#阶段重合度


