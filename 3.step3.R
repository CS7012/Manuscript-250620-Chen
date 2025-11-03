#数据与代码声明
#如果没有购买SCI狂人团队或者生信狂人团队的正版会员
#没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
#如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿

####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
###作者邮箱：sxkrteam@shengxinkuangren.com
#install packages
#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

#install.packages("ggplot2")
rm(list = ls())
#library
library(TwoSampleMR)
library(ggplot2)
library(foreach)

setwd("F:\\133_双细胞孟德尔随机化\\double_cell\\2_imc_out")

iddf=read.table("imc731id.txt",header =T,sep = "\t")


#数据与代码声明
#如果没有购买SCI狂人团队或者生信狂人团队的正版会员
#没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
#如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿

dxw=as.vector(iddf$id)

result=data.frame()

foreach(i=dxw, .errorhandling = "pass") %do%{
  expo_rt<- read.table(file = paste0("localdata/",i,".txt"),header = T,sep = "\t")
  
#数据与代码声明
#如果没有购买SCI狂人团队或者生信狂人团队的正版会员
#没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
#如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿
  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = "GCST90018584.txt.gz",
    sep = "\t",
    snp_col = "rsids2",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value")
  
  
  

  harm_rt <- harmonise_data(
    exposure_dat =  expo_rt, 
    outcome_dat = outc_rt,action=2)
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure*harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  harm_rt$meanf<- mean( harm_rt$f)
  harm_rt<-harm_rt[harm_rt$f>10,]
  
  mr_result<- mr(harm_rt)
  result_or=generate_odds_ratios(mr_result) 
  if (mr_result$pval[3]<0.05){
  result=rbind(result,cbind(id=i,pvalue=result_or$pval[3])) 
  dir.create(i) 
  filename=i
  write.table(harm_rt, file =paste0(filename,"/harmonise.txt"),row.names = F,sep = "\t",quote = F)
  write.table(result_or[,5:ncol(result_or)],file =paste0(filename,"/OR.txt"),row.names = F,sep = "\t",quote = F)
  pleiotropy=mr_pleiotropy_test(harm_rt)
  write.table(pleiotropy,file = paste0(filename,"/pleiotropy.txt"),sep = "\t",quote = F)
  heterogeneity=mr_heterogeneity(harm_rt)
  write.table(heterogeneity,file = paste0(filename,"/heterogeneity.txt"),sep = "\t",quote = F)
  
  #####################################
  #数据与代码声明
  #如果没有购买SCI狂人团队或者生信狂人团队的正版会员
  #没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
  #如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿
  p1 <- mr_scatter_plot(mr_result, harm_rt)
  ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), width=8, height=8)
  
  #####################################
  #数据与代码声明
  #如果没有购买SCI狂人团队或者生信狂人团队的正版会员
  #没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
  #如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿
  
  singlesnp_res<- mr_singlesnp(harm_rt)
  singlesnpOR=generate_odds_ratios(singlesnp_res)
  write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.txt"),row.names = F,sep = "\t",quote = F)
  p2 <- mr_forest_plot(singlesnp_res)
  ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
  sen_res<- mr_leaveoneout(harm_rt)
  p3 <- mr_leaveoneout_plot(sen_res)
  ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), width=8, height=8)
  res_single <- mr_singlesnp(harm_rt)
  p4 <- mr_funnel_plot(singlesnp_res)
  ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
  presso=run_mr_presso(harm_rt,NbDistribution = 1000)
  capture.output(presso,file = paste0(filename,"/presso.txt"))
}
}
write.table(result,"immresult.txt",sep = "\t",quote = F,row.names = F)
#数据与代码声明
#如果没有购买SCI狂人团队或者生信狂人团队的正版会员
#没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
#如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿

