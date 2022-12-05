#!/usr/bin/Rscript

library('stringr')


source('kre_parallel.R')

options(warn=1)

cores = 10
data_names = read.csv('../Virus/Virus_Name.csv',stringsAsFactors = F)$Organism.Abbreviation
family_spe = ''


print(Sys.time())
cat(family_spe,'\n')

n = length(data_names)
raw_data = list()
del_index = c()
for (i in 1:n) {
  temp_virus_data = read.csv(paste0('../Virus/', data_names[i], '.csv'),stringsAsFactors = F)

  sel_criterion = (temp_virus_data[,'Upper.Bound.of.95..CI']-temp_virus_data[,'Lower.Bound.of.95..CI'])<0.2
  temp_virus_data = temp_virus_data[sel_criterion,][,c('Sequence', 'Response.Freq.')]

  # temp_virus_data = temp_virus_data[temp_virus_data[,'Subjects.Tested']>=10,][,c('Sequence', 'Response.Freq.')]

  # cat(paste(data_names[i],':',nrow(temp_virus_data)),'\n')
  if(nrow(temp_virus_data[temp_virus_data$Response.Freq.>0,])<=30){
    del_index = c(del_index, i)}
  else{
    raw_data = c(raw_data, list(temp_virus_data))
  }
}



if(!is.null(del_index)){
  data_names = data_names[-del_index]
  n = n - length(del_index)
}


for (j in 1:length(data_names)) {
  colnames(raw_data[[j]]) = c('pep','prob')
}



virus_dist_matrix = read.csv(paste0('../Virus/', 'Virus_Distance_Data', '.csv'),
  row.names= 1,stringsAsFactors = F)[data_names,data_names]
virus_dist_matrix = as.matrix(virus_dist_matrix)
virus_dist_matrix[upper.tri(virus_dist_matrix)] <- t(virus_dist_matrix)[upper.tri(virus_dist_matrix)]
diag(virus_dist_matrix) = 0
temp_vec = sort(unique(c(virus_dist_matrix)))
all_threshold = '3_threshold'
possible_thresholds = temp_vec[seq(2,length(temp_vec), by = 3)]

# all_threshold = ''
# possible_thresholds = temp_vec[seq(length(temp_vec),2, by = -5)]

cat('Possible Thresholds: ',possible_thresholds,'\n')

raw_l = 1:length(data_names)
seed = 2022
k_fold = 5

ave_trace = data.frame()
total_result = data.frame()
for (i in 1:length(data_names)) {
  # pdf(paste0(latex_path, 'Viruses.pdf'),width = 12,height = 7)
  temp_dist = virus_dist_matrix[i,]
  temp_measure = c()
  cat('\n')
  print(Sys.time())
  cat(paste0('Virus:', data_names[i]),'\n')

  cur_len_train = 0
  for (threshold_index in 1:length(possible_thresholds)) {
    threshold = possible_thresholds[threshold_index]
    # cat(paste0('i=',i))
    # pdf(paste0(latex_path, data_names[i], '.pdf'),width = 12,height = 7)
    
    add_list = which(temp_dist<=threshold&temp_dist>0)

    # if no viruses are in the range, select the closest one
    if(length(add_list)==0){
      add_list = order(temp_dist, decreasing=FALSE)[2]
    }

    if(length(add_list) > cur_len_train){
      ddata1 = Reduce(rbind, raw_data[add_list])
      
      temp = training_fold(ddata1,k_fold = k_fold, grid_len = 20,
        threshold = 0, plot.cv = F, only_config = T, seed = seed,
                     parallel=T, cores = cores)

      tt = pred(ddata1,raw_data[[i]][,1], threshold = 0,config=temp$config,
        parallel=T, cores = cores)
      cc = cor.test(raw_data[[i]][,2],tt[,2], method = 'spearman',
                    alternative = 'greater')
      temp_coe = as.numeric(cc$estimate)
      temp_p_value = as.numeric(cc$p.value)
      cur_len_train = length(add_list)
    }
    

    temp_measure = c(temp_measure, temp_coe)
    total_result = rbind(total_result,cbind(data_names[i], temp_coe, temp_p_value, threshold))
    cat(paste('Threshold:', threshold, ' ',
      'Number of training viruses:', length(add_list),' ',
      'Coe:', temp_coe, ' ',
      'p_value', temp_p_value),'\n')
  }
  ave_trace = rbind(ave_trace, cbind(threshold,
    trunc(median(temp_measure)),
    trunc(mean(temp_measure)), 
    trunc(sd(temp_measure))))

  cat(paste('Measure_coe (median, mean, sd):', trunc(median(temp_measure)),
    trunc(mean(temp_measure)), 
    trunc(sd(temp_measure))),'\n')

  write.csv(total_result,paste0('../result/',family_spe,all_threshold,'total_result.csv'),row.names = T,quote=F)
  write.csv(ave_trace,paste0('../result/',family_spe,all_threshold,'ave_trace.csv'),row.names = T,quote=F)
  
}

