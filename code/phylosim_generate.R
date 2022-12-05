#!/usr/bin/Rscript

rm(list=ls())
library(stringr)
library(ggplot2)
library(ggseqlogo)
library(MASS)
library(Biostrings)
library(phylosim)

set.seed(2021)

one_drive_path = str_replace_all(Sys.getenv("OneDrive"),"\\\\", "/")
latex_path = paste0(one_drive_path,'/kernel_est/new_simulation/code')
if(dir.exists(latex_path)==0){dir.create(latex_path)}
setwd(latex_path)

data(BLOSUM62)

##AA
dict = 'ACDEFGHIKLMNPQRSTVWY'
tmp_index = 1:str_length(dict)
dict = str_sub(dict,tmp_index,tmp_index)

digit_dict = as.character(1:length(dict))
names(digit_dict) = dict
names(dict) = digit_dict 

theta0 = rep(1/length(dict),length(dict))
names(theta0) = dict



generate_pattern<-function(n, len_seq, theta0){
  # a = sample(1:(len_seq-motif_len+1),size = 1,replace=T)
  temp_0 = sample(dict,size = len_seq*n, replace=T,prob = theta0)
  result = c()
  for (i in 1:n) {
    temp_index = 1:len_seq
    temp_ele = temp_0[temp_index]
    temp_0 = temp_0[-temp_index]
    temp_ele = paste0(temp_ele,collapse = '')
    result = c(result, temp_ele)
  }
  
  return(result)
}

evolution<-function(pattern,ntips){
  sim<-PhyloSim(
    phylo=rtree(ntips), # n tips + (n-1)internal nodes
    # we want only tips
    root.seq=AminoAcidSequence(string=pattern,processes=list(list(PAM())))
  );
  Simulate(sim);
  return(sapply(1:sim$ntips, function(x) as.character(sim$sequences[sim$tips][[x]])))
}

# randomly add some letters to make all lengths equal to len_seq

complete<-function(x, dict, len_seq, theta0){
  l = str_length(x)
  if(l<len_seq){
    x = paste0(c(x,sample(dict,size = len_seq-l, replace=T,prob = theta0)),collapse = '')
  }
  return(x)
}

# mean_diag = mean(diag(BLOSUM62[1:20,1:20]))
type = 'local'
print(type)
get_score<-function(pattern_tips,pattern){
  as.numeric(sapply(pattern_tips, function(x) pairwiseAlignment(
    AAString(x) ,AAString(pattern),
    substitutionMatrix = "BLOSUM62",
    scoreOnly=T,type = type)))
}


# length of each sequence 
len_seqs = c(20,40,60)

# the number of sequences
n = 500
print(n)
# the number of pattern
np = 5
len_p = c(10,12,14,16,18)
# pattern_score = runif(np)*10
repl = 1:20
group = c('g1','g2')

# add 5 more patterns; change two patterns
change_pattern_index = c(4,5)
cases = c('standard','more','change')

for (len_seq in len_seqs) {
  for (re in repl) {
    patterns = sapply(len_p, function(x) generate_pattern(1,x,theta0))
    
    for (temp_case in cases) {
      if(temp_case=='standard'){
        for (g in group) {
          ntips = floor(n/np)
          all_tips = lapply(patterns, function(x) evolution(x,ntips))
          
          score = lapply(1:np, function(x) get_score(all_tips[[x]],patterns[x]))
          result = lapply(1:np, function(x) data.frame(old_seq = all_tips[[x]],
                                                       id=x,
                                                       pattern = patterns[x],
                                                       score = score[[x]]))
          result = Reduce(rbind,result)
          
          new_seq = as.character(sapply(result$old_seq, function(x) complete(x, dict, len_seq, theta0)))
          result$new_seq = new_seq
          result = result[,c('new_seq','old_seq','pattern','id','score')]
          write.csv(result, paste0('../data/',n,'_',temp_case,'_',g,'_',
                                   len_seq,'_',re,'.csv'),row.names = F)
        }
      }
      
      if(temp_case=='more'){
        for (g in group) {
          if(g=='g2'){
            more_patterns = sapply(len_p, function(x) generate_pattern(1,x,theta0))
            new_patterns = c(patterns, more_patterns)
            new_np = 2*np
            ntips = floor(n/new_np)
            all_tips = lapply(new_patterns, function(x) evolution(x,ntips))
            
            score = lapply(1:new_np, function(x) get_score(all_tips[[x]],new_patterns[x]))
            result = lapply(1:new_np, function(x) data.frame(old_seq = all_tips[[x]],
                                                         id=x,
                                                         pattern = new_patterns[x],
                                                         score = score[[x]]))
            
            result = Reduce(rbind,result)
            
            new_seq = as.character(sapply(result$old_seq, function(x) complete(x, dict, len_seq, theta0)))
            result$new_seq = new_seq
            result = result[,c('new_seq','old_seq','pattern','id','score')]
            write.csv(result, paste0('../data/',n,'_',temp_case,'_',g,'_',
                                     len_seq,'_',re,'.csv'),row.names = F)
          }else{
            ntips = floor(n/np)
            all_tips = lapply(patterns, function(x) evolution(x,ntips))
            
            score = lapply(1:np, function(x) get_score(all_tips[[x]],patterns[x]))
            result = lapply(1:np, function(x) data.frame(old_seq = all_tips[[x]],
                                                         id=x,
                                                         pattern = patterns[x],
                                                         score = score[[x]]))
            
            result = Reduce(rbind,result)
            
            new_seq = as.character(sapply(result$old_seq, function(x) complete(x, dict, len_seq, theta0)))
            result$new_seq = new_seq
            result = result[,c('new_seq','old_seq','pattern','id','score')]
            write.csv(result, paste0('../data/',n,'_',temp_case,'_',g,'_',
                                     len_seq,'_',re,'.csv'),row.names = F)
          }
        }
      }
      
      if(temp_case=='change'){
        for (g in group) {
          if(g == 'g2'){
            temp_patterns = sapply(len_p[change_pattern_index], function(x) generate_pattern(1,x,theta0))
            new_patterns = patterns
            new_patterns[change_pattern_index] = new_patterns
            ntips = floor(n/np)
            all_tips = lapply(new_patterns, function(x) evolution(x,ntips))
            
            score = lapply(1:np, function(x) get_score(all_tips[[x]],new_patterns[x]))
            result = lapply(1:np, function(x) data.frame(old_seq = all_tips[[x]],
                                                         id=x,
                                                         pattern = new_patterns[x],
                                                         score = score[[x]]))
            
            result = Reduce(rbind,result)
            
            new_seq = as.character(sapply(result$old_seq, function(x) complete(x, dict, len_seq, theta0)))
            result$new_seq = new_seq
            result = result[,c('new_seq','old_seq','pattern','id','score')]
            write.csv(result, paste0('../data/',n,'_',temp_case,'_',g,'_',
                                     len_seq,'_',re,'.csv'),row.names = F)
          }else{
            ntips = floor(n/np)
            all_tips = lapply(patterns, function(x) evolution(x,ntips))
            
            score = lapply(1:np, function(x) get_score(all_tips[[x]],patterns[x]))
            result = lapply(1:np, function(x) data.frame(old_seq = all_tips[[x]],
                                                         id=x,
                                                         pattern = patterns[x],
                                                         score = score[[x]]))
            
            result = Reduce(rbind,result)
            
            new_seq = as.character(sapply(result$old_seq, function(x) complete(x, dict, len_seq, theta0)))
            result$new_seq = new_seq
            result = result[,c('new_seq','old_seq','pattern','id','score')]
            write.csv(result, paste0('../data/',n,'_',temp_case,'_',g,'_',
                                     len_seq,'_',re,'.csv'),row.names = F)
          }
          
        }
      }
      
      
      
    }
    
    
    
  }
}
















