#!/usr/bin/Rscript


library('stringr')
source('kre.R')

data_names = c('BCHK3','CVH22','CVHK1','CVHNL','CVHOC','DENV','MERS1','SARS1','SARS2','ZIKV')
family_pred = c(0.195798385,0.167147788,0.172366922,0.178119399,0.188044877,0.041531542,0.168292519,0.190990265,0.347098799,0.099344387)
family_pred_p = c(3.15E-16,1.76E-55,3.53E-61,5.80E-62,1.06E-78,0.002513632,3.83E-62,2.32E-92,3.08E-275,0.126716987)
all_pred = c(0.066677709,-0.022180776,-0.022582744,-0.02843512,-0.014810541,0.083706986,0.018424202,0.01262478,0.057795159,0.071868356)
all_pred_p = c(0.00317549,0.980392923,0.983850698,0.995661292,0.928306152,7.50E-09,0.035624203,0.091295488,5.34E-09,0.204626171)

sel_order = c('SARS1','SARS2','CVHOC','MERS1','CVHK1','CVH22','CVHNL','BCHK3','DENV','ZIKV')
sel_order_index = sapply(sel_order, function(x) which(data_names==x))
result = data.frame(Virus = data_names, family_pred = family_pred, all_pred = all_pred,
                    family_pred_p = family_pred_p, all_pred_p = all_pred_p)
result = result[sel_order_index,]
p1 = ggplot(data = result)+theme_bw()+
	geom_point(aes(x=factor(data_names),y=family_pred),color='red',shape=15,size=3)+geom_line(aes(x=factor(data_names),y=family_pred,group=1),color='red')+
	geom_point(aes(x=factor(data_names),y=all_pred),color='blue',shape=15,size=3)+geom_line(aes(x=factor(data_names),y=all_pred,group=1),color='blue')+
	# theme(legend.position="none")+
	xlab('Virus')+ylab('Coefficient')

types = c(rep('Family_specific',length(data_names)), rep('General',length(data_names)))
result2 = data.frame(Virus = c(data_names,data_names), predictions = c(family_pred,all_pred), 
                     significance= c(family_pred_p, all_pred_p), type_pred = types)
p2 = ggplot(data=result2, aes(x=factor(Virus), y=predictions, fill=factor(type_pred))) +theme_bw()+
  geom_bar(stat="identity", position=position_dodge(0.5), width=0.5, alpha=1)+
  # scale_fill_manual(values=c('red','blue'))+
  # theme(legend.position="none")+
  xlab('Virus')+ylab('Coefficient')+labs(fill = 'Prediction Type')+ylim(c(-0.05,0.4))+
  geom_text(aes(label = significance),position=position_dodge(width = 0.5),size = 3,vjust = -0.5)
p2





