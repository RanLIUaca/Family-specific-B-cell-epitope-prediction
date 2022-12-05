#!/usr/bin/Rscript

rm(list=ls())
library(stringr)
library(ggplot2)
library(MASS)
library(patchwork)

# one_drive_path = str_replace_all(Sys.getenv("OneDrive"),"\\\\", "/")
# setwd(paste0(one_drive_path,'/kernel_est/patients/code'))



################## first pic
n = 10000

p0 = c(0, 1)
p1 = c(2, 0.8)
d =  0.8

set.seed(2021)

x0 = seq(-3, 5, length.out = n)
y1 = dnorm(x0, mean = p0[1], sd = p0[2])
y2 = dnorm(x0, mean = p1[1], sd = p1[2])

data1 = data.frame(x = x0, y1 = y1, y2 = y2)
ssize = 2
sssize = 10
pdf(paste0(one_drive_path,'/kernel_est/patients/pics/sero/','dist_1.pdf'),width = 10,height = 5)
f1<-ggplot(data = data1, aes(x = x))+
  geom_line(aes(y = y1), color = "#f8766d",size = ssize, linetype = 'solid')+
  geom_line(aes(y = y2), color = "#00bfc4",size = ssize, linetype = 'solid')+
  ylab("") + xlab("") + 
  scale_shape_manual(values=c(1,2))+
  geom_vline(xintercept = d, color = 'brown4',size = ssize-0.5, linetype = 'longdash')+
  geom_ribbon(aes(ymin = 0, ymax = ifelse(x>d &x<max(x0), y1, 0)), fill="#f8766d", alpha = 0.3)+
  geom_ribbon(aes(ymin = 0, ymax = ifelse(x<d &x>min(x0), y2, 0)), fill="#00bfc4", alpha = 0.3)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  theme(panel.border = element_rect(color = 'grey'))+
  annotate("text", x = -2, y = 0.35, label = "Non-epitope", color = "#f8766d", size = sssize)+
  annotate("text", x = 3.5+0.2, y = 0.35, label = "Epitope", color = "#00bfc4", size = sssize)
f1
# f1 + plot_annotation(tag_levels = list(c("A"))) 
dev.off()





################## second pic
n = 10000

p0 = c(0, 1)
p1 = c(0.4, 0.8)
d =  0.2

set.seed(2021)

x0 = seq(-3, p1[1]+3, length.out = n)
y1 = dnorm(x0, mean = p0[1], sd = p0[2])
y2 = dnorm(x0, mean = p1[1], sd = p1[2])

data1 = data.frame(x = x0, y1 = y1, y2 = y2)
ssize = 2
sssize = 10
pdf(paste0(one_drive_path,'/kernel_est/patients/pics/sero/','dist_2.pdf'),width = 10,height = 5)
f2<-ggplot(data = data1, aes(x = x))+
  geom_line(aes(y = y1), color = "#f8766d",size = ssize, linetype = 'solid')+
  geom_line(aes(y = y2), color = "#00bfc4",size = ssize, linetype = 'solid')+
  ylab("") + xlab("") + 
  scale_shape_manual(values=c(1,2))+
  geom_vline(xintercept = d, color = 'brown4',size = ssize-0.5, linetype = 'longdash')+
  geom_ribbon(aes(ymin = 0, ymax = ifelse(x>d &x<max(x0), y1, 0)), fill="#f8766d", alpha = 0.3)+
  geom_ribbon(aes(ymin = 0, ymax = ifelse(x<d &x>min(x0), y2, 0)), fill="#00bfc4", alpha = 0.3)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  theme(panel.border = element_rect(color = 'grey'))+
  annotate("text", x = -1.5-0.2, y = 0.35, label = "Non-epitope", color = "#f8766d", size = sssize)+
  annotate("text", x = p1[1]+1.5, y = 0.35, label = "Epitope", color = "#00bfc4", size = sssize)
f2
# f2+ plot_annotation(tag_levels = list(c("B"))) 
dev.off()





