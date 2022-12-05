#!/usr/bin/Rscript

rm(list=ls())
library(stringr)
library(ggplot2)
library(MASS)

# one_drive_path = str_replace_all(Sys.getenv("OneDrive"),"\\\\", "/")
# setwd(paste0(one_drive_path,'/kernel_est/patients/code'))


patient = read.csv('../data/SARS-CoV-2+Patient.csv',header = F)
patient = as.matrix(patient[2:nrow(patient),3:ncol(patient)])
mode(patient) = 'numeric'
colnames(patient) = paste0('p',seq(1,ncol(patient)))
row.names(patient) = seq(1,nrow(patient))

health = read.csv('../data/SARS-CoV-2+Health.csv',header = F)
health = as.matrix(health[2:nrow(health),3:ncol(health)])
mode(health) = 'numeric'
colnames(health) = paste0('h',seq(1,ncol(health)))
row.names(health) = seq(1,nrow(health))



hs = stack(as.data.frame(t(health)))
ps = stack(as.data.frame(t(patient)))

threshold = 0.5
ps[ps[,1]>threshold,1] = threshold


hs$State = 'Health'
ps$State = 'Patient'

pdf(paste0(one_drive_path,'/kernel_est/patients/pics/sero/','ana_value.pdf'),width = 10,height = 5)
total_result = rbind(hs,ps)
p3<-ggplot(total_result, aes(x=ind, y=values,color = State)) +
  geom_boxplot(aes(fill = State))+
  theme_light()+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
  )+
  labs(x="Peptide Index",y="Immune Response")
p3
dev.off()




cal_threshold<-function(row_pep){
  return(mean(row_pep)+3*sd(row_pep))
}

thresholds = apply(health, 1, cal_threshold)
patient[patient>=thresholds] = 1
patient[patient<thresholds] = 0

cal_inconsis_ratio<-function(row_pep){
  # return(1-abs(sum(row_pep==1) - sum(row_pep==0))/length(row_pep))
  return(sum(row_pep==0)/(length(row_pep)))
  # return(sum(row_pep==0)*sum(row_pep==1)/(length(row_pep))**2)
}

alpha = 0.3

# plot a figure for count changes along with alpha
# scatter plot (specific alpha = 0.3)
pdf(paste0(one_drive_path,'/kernel_est/patients/pics/sero/','inconsis.pdf'),width = 10,height = 5)
ratios = apply(patient, 1, cal_inconsis_ratio)
inconsis = sum(ratios>alpha & ratios<=1-alpha)
result = data.frame(ratios)
colnames(result) = 'inconsis_ratio'
result['pep_i'] = 1:nrow(result)
ggplot(result,aes(x=pep_i,y=inconsis_ratio)) +
  geom_bar(stat="identity",fill="#f8766d", alpha=.9, width = 0.7) +
  # geom_ribbon(aes(x = seq(-1,43, length.out = 42), ymin=alpha, ymax=1-alpha), fill="#00bfc4", alpha = 0.3)+
  annotate("rect", xmin = 0, xmax = 43, ymin = alpha, ymax = 1-alpha,
           fill="#00bfc4", alpha = 0.3)+
  annotate("text", x = 45, y = 0.5, label = 'Inconsistent Zone',
           color="#00bfc4",alpha = 1, size = 7)+
  annotate("text", x = 44, y = alpha-0.03, label = 'alpha',parse = T,
           color="#00bfc4",alpha = 1, size = 7)+
  annotate("text", x = 44, y = 1 - alpha +0.04, label = '1 - alpha', parse = T,
           color="#00bfc4",alpha = 1, size = 7)+
  # geom_point(colour ="#f8766d", alpha=.9, size =3, shape = 17)+
  coord_flip() +
  xlim(0,45)+
  xlab("") +
  ylab('Non-response Ratio')+
  theme_light()+
  geom_hline(aes(yintercept=1-alpha), colour="#00bfc4", linetype="dashed", size=1.3)+
  geom_hline(aes(yintercept=alpha), colour="#00bfc4", linetype="dashed", size=1.3)

print(inconsis)
dev.off()


alphas = seq(0,0.5,by=0.05)
results = signif(sapply(alphas,function(alpha){sum(ratios>alpha & ratios<=1-alpha)/nrow(patient)}),
                 digit = 3)
results = t(data.frame(Alpha = signif(alphas,digit = 3),Inconsis_ratio = results))
results
write.table(results,file=paste0(one_drive_path,'/kernel_est/patients/pics/sero/','inconsis_alpha.csv'),
          quote=F, col.names = F, sep = ',')
# plot(results)





