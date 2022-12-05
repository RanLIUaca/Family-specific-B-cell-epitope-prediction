library(ggplot2)
library(ggforce)
library(ggrepel)
library(stringr)
library(tidyr)

one_drive_path = str_replace_all(Sys.getenv("OneDrive"),"\\\\", "/")

setwd(paste0(one_drive_path,'/kernel_est/new_virus/code'))
data_names = read.csv('../Virus/Virus_Name.csv',stringsAsFactors = F)$Organism.Abbreviation

latex_path = paste0(one_drive_path,'/kernel_est/new_virus/pics/')
if(dir.exists(latex_path)==0){dir.create(latex_path)}


sel_virus = 'SARS2'


n = length(data_names)
raw_data = list()
del_index = c()
for (i in 1:n) {
  temp_virus_data = read.csv(paste0('../Virus/', data_names[i], '.csv'),stringsAsFactors = F)

  sel_criterion = (temp_virus_data[,'Upper.Bound.of.95..CI']-temp_virus_data[,'Lower.Bound.of.95..CI'])<0.2
  temp_virus_data = temp_virus_data[sel_criterion,][,c('Sequence', 'Response.Freq.')]
  temp_virus_data = temp_virus_data[!duplicated(temp_virus_data$Sequence),]

  if(nrow(temp_virus_data[temp_virus_data$Response.Freq.>0,])<=30){
    del_index = c(del_index, i)}
}


if(!is.null(del_index)){
  data_names = data_names[-del_index]
  n = n - length(del_index)
}

const = 0

virus_dist_matrix = read.csv(paste0('../Virus/', 'Virus_Distance_Data', '.csv'),
  row.names= 1,stringsAsFactors = F)[data_names,data_names]
virus_dist_matrix = as.matrix(virus_dist_matrix)
virus_dist_matrix[upper.tri(virus_dist_matrix)] <- t(virus_dist_matrix)[upper.tri(virus_dist_matrix)]
virus_dist_matrix = virus_dist_matrix+const
diag(virus_dist_matrix) = 0


GenetDis<-as.dist(virus_dist_matrix,diag=FALSE)
mds.AntiDis<-cmdscale(GenetDis,k=2,add= TRUE)$points
sel_virus_index = which(rownames(virus_dist_matrix)==sel_virus)
mds.AntiDis[,1]=mds.AntiDis[,1]-mds.AntiDis[sel_virus_index,1]
mds.AntiDis[,2]=mds.AntiDis[,2]-mds.AntiDis[sel_virus_index,2]


far_index = which(data_names=='CVHOC')
best_t = sqrt(sum(mds.AntiDis[far_index,]**2))+1

AntiDisdata<-data.frame(name = rownames(virus_dist_matrix),
                        PC1=-mds.AntiDis[,1],
                        PC2=mds.AntiDis[,2],
                        distance = virus_dist_matrix[sel_virus_index,])

t <- seq(0, 2*pi, length.out = 500) # 50 data points from 0 up to 2*pi (angles)
x <- sin(t)
y <- cos(t)
circ1 <- data.frame(t, x1=x, y1=y)*best_t

fig3<- ggplot() + theme_bw() +
  geom_point(data=AntiDisdata,  mapping=aes(x=PC1,y=PC2,color = distance), alpha=0.8, size=3)  + 
  scale_color_gradient(low="blue", high="red")+
  geom_path(data=circ1, aes(x=x1,y=y1),color="darkred",size=1)+
  # theme(legend.position = "none")+
  # scale_x_continuous(breaks = seq(-4, 6, by = 2))+
  # scale_y_continuous(breaks = seq(-4, 6, by = 2))+ 
  labs(x = "", y = "") +
  coord_fixed() + 
  geom_text_repel(data=AntiDisdata,aes(x=PC1, y=PC2,label=name),max.overlaps=100)
fig3

pdf(paste0(latex_path,'new_viruses/',sel_virus,'_dist_show.pdf'),width = 10,height = 5)
fig3
dev.off()
