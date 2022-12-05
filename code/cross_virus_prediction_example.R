
library(stringr)

source('kre.R')


set.seed(2021)

sars_data = read.csv('../data/SARS-CoV+-+B+cell+epitope.csv',sep = '',stringsAsFactors = F)
covid_data = read.csv('../data/SARS-CoV-2+-+B+cell+epitope.csv',sep = '',stringsAsFactors = F)
denv_data = read.csv('../data/DENV_B_cell_epitope.csv',sep = ',',header = F,stringsAsFactors = F)[,2:3]
zika_data = read.csv('../data/ZIKA_B_cell_epitope.csv',sep = ',',header = F,stringsAsFactors = F)[,2:3]

raw_data = c(list(sars_data),list(covid_data),list(denv_data),list(zika_data))
data_names = c('SARS-CoV','SARS-CoV-2','DENV','ZIKA')
for (j in 1:length(data_names)) {
  colnames(raw_data[[j]]) = c('pep','prob')
  raw_data[[j]]['virus'] = data_names[j]
}
totoal_data = Reduce(rbind, raw_data)

new_result = totoal_data[!duplicated(totoal_data$pep),]
raw_data = list()
for (i in 1:4) {
  temp_data = new_result[new_result['virus']==data_names[i],]
  raw_data = c(raw_data, list(temp_data[,c('pep','prob')]))
}


raw_data = list(sars_data,covid_data,denv_data,zika_data)
set.seed(2021)
k_fold = 5
result_coe_matrix = matrix(0,nrow = length(raw_data),ncol = length(raw_data))
result_p_value = matrix(0,nrow = length(raw_data),ncol = length(raw_data))
optimal_h = c()

# 'Logistic', 'Sigmoid', 'Gaussian'
K = kernel_sel('Gaussian')
pdf(paste0('','self_cv.pdf'),width = 12,height = 7)
par(mfrow=c(1,4))
for (i in 1:length(data_names)) {
  temp = training_fold(raw_data[[i]],k_fold = k_fold, threshold = 0,plot.cv = T,K = K)
  # temp = training(raw_data[[i]], threshold = 0,plot.cv = T,K = K)
  
  cv_result = temp$result
  
  result_coe_matrix[i,i] = cor(cv_result[,2],cv_result[,3],method = 'spearman')
  cc = cor.test(cv_result[,2],cv_result[,3],method = 'spearman',
                alternative = 'greater')
  result_p_value[i,i] = cc$p.value	
  
  optimal_h = c(optimal_h, temp$config$bw)
  for (j in 1:length(data_names)) {
    print(c(i,j))
    if(i!=j){
      tt = pred(raw_data[[i]],raw_data[[j]][,1], threshold = 0,config=temp$config,K = K)
      result_coe_matrix[i,j] = cor(raw_data[[j]][,2],tt[,2],method = 'spearman')
      cc = cor.test(raw_data[[j]][,2],tt[,2],method = 'spearman',
                    alternative = 'greater')
      result_p_value[i,j] = cc$p.value	
    }
  }
}
dev.off()
optimal_h

combine_result = sapply(1:length(data_names), function(x){
  paste0(trunc(result_coe_matrix[,x]),'(',trunc(result_p_value[,x]),')')
})

name_matrix<-function(A){
  row.names(A) = data_names
  colnames(A) = data_names
  return(A)
}
# combine_result = paste0(trunc(result_coe_matrix),'(',trunc(result_p_value),')')
combine_result = name_matrix(combine_result)

combine_result
optimal_h

# write.csv(optimal_h,paste0('optimal_h.csv'),row.names = F,quote=F)
# write.csv(combine_result,paste0('combine_result.csv'),quote=F)




