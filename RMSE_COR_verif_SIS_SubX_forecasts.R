# TENGO QUE TERMINAR DE DEFINIR COMO HAGO PARA QUE ESTO NO ME SOBREESCRIBA LA DT DEL FCST DEL SIS PORQUE
# DESPUES NO PUEDO CALCULAR LAS ROC NEGATIVAS



# This scripts plots all the SIS index forecasts for the SubX models, the MME forecast and the
# MME spread (using all SubX models' ensemble mean)
#
# M. Alvarez - 2020
#_____________________________________________________________________________________________
rm(list=ls())
graphics.off()

setwd("/home/alvarez/")

# Call libraries to be used
library("ncdf4")
library("metR")
library('pracma')
library('lubridate')
library('reshape2')
library('data.table')
library("s2dverification") # No se si lo voy a usar al final
library("RColorBrewer")
library("maps")
library("ggplot2")
library("gridExtra")
library("grid")
library("verification")

#---------------------------------------------------------------------------------------
#  Main Program  
#---------------------------------------------------------------------------------------
load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_MME.rda")
dt.fcstSIS.MME = dt.fcstSIS
setnames(dt.fcstSIS.MME,old="SIS.index.fcst",new="SIS.fcst.MME")

rm(dt.fcstSIS)

load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_ECCC.rda")
dt.fcstSIS.ECCC = dt.fcstSIS
setnames(dt.fcstSIS.ECCC,old="SIS.index.fcst",new="SIS.fcst.ECCC")
rm(dt.fcstSIS)
load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_EMC.rda")
dt.fcstSIS.EMC = dt.fcstSIS
setnames(dt.fcstSIS.EMC,old="SIS.index.fcst",new="SIS.fcst.EMC")
rm(dt.fcstSIS)
load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_ESRL.rda")
dt.fcstSIS.ESRL = dt.fcstSIS
setnames(dt.fcstSIS.ESRL,old="SIS.index.fcst",new="SIS.fcst.ESRL")
rm(dt.fcstSIS)
load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_NRL.rda")
dt.fcstSIS.NRL = dt.fcstSIS
setnames(dt.fcstSIS.NRL,old="SIS.index.fcst",new="SIS.fcst.NRL")
rm(dt.fcstSIS)

load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_GMAO.rda")
dt.fcstSIS.GMAO = dt.fcstSIS
setnames(dt.fcstSIS.GMAO,old="SIS.index.fcst.run1",new="SIS.fcst.GMAO.run1")
setnames(dt.fcstSIS.GMAO,old="SIS.index.fcst.run2",new="SIS.fcst.GMAO.run2")
rm(dt.fcstSIS)
load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_RSMAS.rda")
dt.fcstSIS.RSMAS = dt.fcstSIS
setnames(dt.fcstSIS.RSMAS,old="SIS.index.fcst.run1",new="SIS.fcst.RSMAS.run1")
setnames(dt.fcstSIS.RSMAS,old="SIS.index.fcst.run2",new="SIS.fcst.RSMAS.run2")
rm(dt.fcstSIS)

# Cargo el RTSIS para la verificación
load("~/SubX_scripts/RTSIS_OA9899_1415.rda")


dt.fcstSIS=merge(dt.fcstSIS.ECCC,dt.fcstSIS.EMC,by=c("targetdate","startdate","L"),all=TRUE)
dt.fcstSIS=merge(dt.fcstSIS,dt.fcstSIS.ESRL,by=c("targetdate","startdate","L"),all=TRUE)
dt.fcstSIS=merge(dt.fcstSIS,dt.fcstSIS.NRL,by=c("targetdate","startdate","L"),all=TRUE)
dt.fcstSIS=merge(dt.fcstSIS,dt.fcstSIS.GMAO,by=c("targetdate","startdate","L"),all=TRUE)
dt.fcstSIS=merge(dt.fcstSIS,dt.fcstSIS.RSMAS,by=c("targetdate","startdate","L"),all=TRUE)
dt.fcstSIS=merge(dt.fcstSIS,dt.fcstSIS.MME,by=c("targetdate","startdate","L"),all=TRUE)

# 
# bin.obs.SIS = dt.RTSIS
# umb=1
# bin.obs.SIS = bin.obs.SIS[,bin.obs := SIS.index>=umb]
# bin.obs.SIS$bin.obs.num=0
# bin.obs.SIS[bin.obs==TRUE]$bin.obs.num = 1

# Join observed SIS values
dt.verif.SIS=merge(dt.fcstSIS,dt.RTSIS[targetdate>as.Date("1999-01-07")],by="targetdate",all.x=TRUE)

# Computation of RMSE for the ensemble mean and each member according to lead

rmse.mme=seq(1,28,1)
rmse.eccc=seq(1,28,1)
rmse.emc=seq(1,28,1)
rmse.esrl=seq(1,28,1)
rmse.nrl=seq(1,28,1)
rmse.gmao=seq(1,28,1)
rmse.rsmas=seq(1,28,1)

for(lead in 1:28){
  rmse.mme[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.MME)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.MME)))
  # Now for each member (without considering run2)
  rmse.eccc[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.ECCC)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.ECCC)))
  rmse.emc[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.EMC)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.EMC)))
  rmse.esrl[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.ESRL)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.ESRL)))
  rmse.nrl[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.NRL)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.NRL)))
  rmse.gmao[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.GMAO.run1)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.GMAO.run1)))
  rmse.rsmas[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.RSMAS.run1)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.RSMAS.run1)))
   
}

dt.rmse=data.table(
  L=seq(1,28,1),
  rmse.eccc=rmse.eccc,
  rmse.emc=rmse.emc,
  rmse.esrl=rmse.esrl,
  rmse.nrl=rmse.nrl,
  rmse.gmao=rmse.gmao,
  rmse.rsmas=rmse.rsmas,
  rmse.mme=rmse.mme
)

colorBlindBlack8  <- c("MME" = "#000000","ECCC" = "#E69F00","EMC" = "#56B4E9","ESRL" = "#009E73", 
                      "NRL" = "#CC79A7","GMAO" ="#0072B2","RSMAS" = "#D55E00")

g1 <- ggplot(data=dt.rmse)+
  geom_line(aes(x=L,y=rmse.mme,color="MME"),size=1.2)+
  geom_line(aes(x=L,y=rmse.eccc,color="ECCC"))+
  geom_line(aes(x=L,y=rmse.emc,color="EMC"))+
  geom_line(aes(x=L,y=rmse.esrl,color="ESRL"))+
  geom_line(aes(x=L,y=rmse.nrl,color="NRL"))+
  geom_line(aes(x=L,y=rmse.gmao,color="GMAO"))+
  geom_line(aes(x=L,y=rmse.rsmas,color="RSMAS"))+
  labs(y= "RMSE", x = "Lead time (days)")+
  scale_color_manual(values = colorBlindBlack8)+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0.5),breaks=seq(1,28,2)) + 
  scale_y_continuous(limits=c(0.5,1.3),breaks=seq(0.5,1.3,0.2),expand=expand_scale(mult = c(0.05, .1)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=15), axis.title = element_text(size=15))+
  theme(legend.position=c(0.87,0.25),legend.title=element_blank())

# Computation of CORR for the ensemble mean and each member according to lead

corr.mme=seq(1,28,1)
corr.eccc=seq(1,28,1)
corr.emc=seq(1,28,1)
corr.esrl=seq(1,28,1)
corr.nrl=seq(1,28,1)
corr.gmao=seq(1,28,1)
corr.rsmas=seq(1,28,1)

for(lead in 1:28){
  corr.mme[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.MME,use="pairwise.complete.obs")
  # Now for each member (without considering run2)
  corr.eccc[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.ECCC,use="pairwise.complete.obs")
  corr.emc[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.EMC,use="pairwise.complete.obs")
  corr.esrl[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.ESRL,use="pairwise.complete.obs")
  corr.nrl[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.NRL,use="pairwise.complete.obs")
  corr.gmao[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.GMAO.run1,use="pairwise.complete.obs")
  corr.rsmas[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.RSMAS.run1,use="pairwise.complete.obs")
  
}

dt.corr=data.table(
  L=seq(1,28,1),
  corr.eccc=corr.eccc,
  corr.emc=corr.emc,
  corr.esrl=corr.esrl,
  corr.nrl=corr.nrl,
  corr.gmao=corr.gmao,
  corr.rsmas=corr.rsmas,
  corr.mme=corr.mme
)

colorBlindBlack8  <- c("MME" = "#000000","ECCC" = "#E69F00","EMC" = "#56B4E9","ESRL" = "#009E73", 
                       "NRL" = "#CC79A7","GMAO" ="#0072B2","RSMAS" = "#D55E00")

g2 <- ggplot(data=dt.corr)+
  geom_line(aes(x=L,y=corr.mme,color="MME"),size=1.2)+
  geom_line(aes(x=L,y=corr.eccc,color="ECCC"))+
  geom_line(aes(x=L,y=corr.emc,color="EMC"))+
  geom_line(aes(x=L,y=corr.esrl,color="ESRL"))+
  geom_line(aes(x=L,y=corr.nrl,color="NRL"))+
  geom_line(aes(x=L,y=corr.gmao,color="GMAO"))+
  geom_line(aes(x=L,y=corr.rsmas,color="RSMAS"))+
  labs(y= "COR", x = "Lead time (days)")+
  scale_color_manual(values = colorBlindBlack8)+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0.5),breaks=seq(1,28,2)) + 
  scale_y_continuous(limits=c(0,0.9),breaks=seq(0,1,0.2),minor_breaks =seq(0,1,0.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=15), axis.title = element_text(size=15))+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))+
  theme(legend.position=c(0.87,0.75),legend.title=element_blank())

fig <- grid.arrange(g1,g2, ncol = 2,top = textGrob(paste0("SubX SIS index forecast, all inits"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/alvarez/SubX_processed_Rdata/verif_SIS_index_forecast_all.png"),plot=fig,width = 10, height = 5)

#---------------------------------------------------------------------------------------------------------------
### RMSE and CORR according to months of init

# ONLY DJF

# Computation of RMSE for the ensemble mean and each member according to lead

rmse.mme=seq(1,28,1)
rmse.eccc=seq(1,28,1)
rmse.emc=seq(1,28,1)
rmse.esrl=seq(1,28,1)
rmse.nrl=seq(1,28,1)
rmse.gmao=seq(1,28,1)
rmse.rsmas=seq(1,28,1)

for(lead in 1:28){
  rmse.mme[lead]=sqrt(sum((dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index-dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.MME)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.MME)))
  # Now for each member (without considering run2)
  rmse.eccc[lead]=sqrt(sum((dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index-dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.ECCC)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.ECCC)))
  rmse.emc[lead]=sqrt(sum((dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index-dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.EMC)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.EMC)))
  rmse.esrl[lead]=sqrt(sum((dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index-dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.ESRL)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.ESRL)))
  rmse.nrl[lead]=sqrt(sum((dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index-dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.NRL)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.NRL)))
  rmse.gmao[lead]=sqrt(sum((dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index-dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.GMAO.run1)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.GMAO.run1)))
  rmse.rsmas[lead]=sqrt(sum((dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index-dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.RSMAS.run1)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.RSMAS.run1)))
  
}

dt.rmse=data.table(
  L=seq(1,28,1),
  rmse.eccc=rmse.eccc,
  rmse.emc=rmse.emc,
  rmse.esrl=rmse.esrl,
  rmse.nrl=rmse.nrl,
  rmse.gmao=rmse.gmao,
  rmse.rsmas=rmse.rsmas,
  rmse.mme=rmse.mme
)

colorBlindBlack8  <- c("MME" = "#000000","ECCC" = "#E69F00","EMC" = "#56B4E9","ESRL" = "#009E73", 
                       "NRL" = "#CC79A7","GMAO" ="#0072B2","RSMAS" = "#D55E00")

g1 <- ggplot(data=dt.rmse)+
  geom_line(aes(x=L,y=rmse.mme,color="MME"),size=1.2)+
  geom_line(aes(x=L,y=rmse.eccc,color="ECCC"))+
  geom_line(aes(x=L,y=rmse.emc,color="EMC"))+
  geom_line(aes(x=L,y=rmse.esrl,color="ESRL"))+
  geom_line(aes(x=L,y=rmse.nrl,color="NRL"))+
  geom_line(aes(x=L,y=rmse.gmao,color="GMAO"))+
  geom_line(aes(x=L,y=rmse.rsmas,color="RSMAS"))+
  labs(y= "RMSE", x = "Lead time (days)")+
  scale_color_manual(values = colorBlindBlack8)+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0.5),breaks=seq(1,28,2)) + 
  scale_y_continuous(limits=c(0.5,1.3), breaks=seq(0.5,1.3,0.2),expand=expand_scale(mult = c(0.05, .1)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=15), axis.title = element_text(size=15))+
  theme(legend.position=c(0.87,0.25),legend.title=element_blank())

# Computation of CORR for the ensemble mean and each member according to lead

corr.mme=seq(1,28,1)
corr.eccc=seq(1,28,1)
corr.emc=seq(1,28,1)
corr.esrl=seq(1,28,1)
corr.nrl=seq(1,28,1)
corr.gmao=seq(1,28,1)
corr.rsmas=seq(1,28,1)

for(lead in 1:28){
  corr.mme[lead]=cor(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index,dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.MME,use="pairwise.complete.obs")
  # Now for each member (without considering run2)
  corr.eccc[lead]=cor(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index,dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.ECCC,use="pairwise.complete.obs")
  corr.emc[lead]=cor(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index,dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.EMC,use="pairwise.complete.obs")
  corr.esrl[lead]=cor(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index,dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.ESRL,use="pairwise.complete.obs")
  corr.nrl[lead]=cor(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index,dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.NRL,use="pairwise.complete.obs")
  corr.gmao[lead]=cor(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index,dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.GMAO.run1,use="pairwise.complete.obs")
  corr.rsmas[lead]=cor(dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.index,dt.verif.SIS[L==lead & (month(startdate)==1 | month(startdate)==2 | month(startdate)==12),]$SIS.fcst.RSMAS.run1,use="pairwise.complete.obs")
  
}

dt.corr=data.table(
  L=seq(1,28,1),
  corr.eccc=corr.eccc,
  corr.emc=corr.emc,
  corr.esrl=corr.esrl,
  corr.nrl=corr.nrl,
  corr.gmao=corr.gmao,
  corr.rsmas=corr.rsmas,
  corr.mme=corr.mme
)

colorBlindBlack8  <- c("MME" = "#000000","ECCC" = "#E69F00","EMC" = "#56B4E9","ESRL" = "#009E73", 
                       "NRL" = "#CC79A7","GMAO" ="#0072B2","RSMAS" = "#D55E00")

g2 <- ggplot(data=dt.corr)+
  geom_line(aes(x=L,y=corr.mme,color="MME"),size=1.2)+
  geom_line(aes(x=L,y=corr.eccc,color="ECCC"))+
  geom_line(aes(x=L,y=corr.emc,color="EMC"))+
  geom_line(aes(x=L,y=corr.esrl,color="ESRL"))+
  geom_line(aes(x=L,y=corr.nrl,color="NRL"))+
  geom_line(aes(x=L,y=corr.gmao,color="GMAO"))+
  geom_line(aes(x=L,y=corr.rsmas,color="RSMAS"))+
  labs(y= "COR", x = "Lead time (days)")+
  scale_color_manual(values = colorBlindBlack8)+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0.5),breaks=seq(1,28,2)) + 
  scale_y_continuous(limits=c(0,0.9),breaks=seq(0,1,0.2),minor_breaks =seq(0,1,0.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=15), axis.title = element_text(size=15))+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))+
  theme(legend.position=c(0.87,0.75),legend.title=element_blank())

fig <- grid.arrange(g1,g2, ncol = 2,top = textGrob(paste0("SubX SIS index forecast, DJF inits"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/alvarez/SubX_processed_Rdata/verif_SIS_index_forecast_DJFinits.png"),plot=fig,width = 10, height = 5)

#---------------------------------------------
# Computation of RMSE for the ensemble mean and each member according to lead

sis.in.inits=dt.verif.SIS[L==0]$SIS.index
# 486 start dates
# 252 SIS positive
# 71 SIS>1
# 71 SIS<-1

hist(sis.in.inits)

umb=-1 

sel.startdates = dt.verif.SIS[L==0 & SIS.index<=umb]$startdate

dt.verif.SIS=dt.verif.SIS[startdate %in% sel.startdates,]

rmse.mme=seq(1,28,1)
rmse.eccc=seq(1,28,1)
rmse.emc=seq(1,28,1)
rmse.esrl=seq(1,28,1)
rmse.nrl=seq(1,28,1)
rmse.gmao=seq(1,28,1)
rmse.rsmas=seq(1,28,1)

for(lead in 1:28){
  rmse.mme[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.MME)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.MME)))
  # Now for each member (without considering run2)
  rmse.eccc[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.ECCC)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.ECCC)))
  rmse.emc[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.EMC)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.EMC)))
  rmse.esrl[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.ESRL)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.ESRL)))
  rmse.nrl[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.NRL)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.NRL)))
  rmse.gmao[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.GMAO.run1)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.GMAO.run1)))
  rmse.rsmas[lead]=sqrt(sum((dt.verif.SIS[L==lead,]$SIS.index-dt.verif.SIS[L==lead,]$SIS.fcst.RSMAS.run1)^2,na.rm=TRUE)/sum(!is.na(dt.verif.SIS[L==lead,]$SIS.fcst.RSMAS.run1)))
  
}

dt.rmse=data.table(
  L=seq(1,28,1),
  rmse.eccc=rmse.eccc,
  rmse.emc=rmse.emc,
  rmse.esrl=rmse.esrl,
  rmse.nrl=rmse.nrl,
  rmse.gmao=rmse.gmao,
  rmse.rsmas=rmse.rsmas,
  rmse.mme=rmse.mme
)

colorBlindBlack8  <- c("MME" = "#000000","ECCC" = "#E69F00","EMC" = "#56B4E9","ESRL" = "#009E73", 
                       "NRL" = "#CC79A7","GMAO" ="#0072B2","RSMAS" = "#D55E00")

g1 <- ggplot(data=dt.rmse)+
  geom_line(aes(x=L,y=rmse.mme,color="MME"),size=1.2)+
  geom_line(aes(x=L,y=rmse.eccc,color="ECCC"))+
  geom_line(aes(x=L,y=rmse.emc,color="EMC"))+
  geom_line(aes(x=L,y=rmse.esrl,color="ESRL"))+
  geom_line(aes(x=L,y=rmse.nrl,color="NRL"))+
  geom_line(aes(x=L,y=rmse.gmao,color="GMAO"))+
  geom_line(aes(x=L,y=rmse.rsmas,color="RSMAS"))+
  labs(y= "RMSE", x = "Lead time (days)")+
  scale_color_manual(values = colorBlindBlack8)+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0.5),breaks=seq(1,28,2)) + 
  scale_y_continuous(limits=c(0.5,1.3),breaks=seq(0.5,1.3,0.2),expand=expand_scale(mult = c(0.05, .1)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=15), axis.title = element_text(size=15))+
  theme(legend.position=c(0.6,0.1),legend.title=element_blank(),legend.direction = "horizontal")

# Computation of CORR for the ensemble mean and each member according to lead

corr.mme=seq(1,28,1)
corr.eccc=seq(1,28,1)
corr.emc=seq(1,28,1)
corr.esrl=seq(1,28,1)
corr.nrl=seq(1,28,1)
corr.gmao=seq(1,28,1)
corr.rsmas=seq(1,28,1)

for(lead in 1:28){
  corr.mme[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.MME,use="pairwise.complete.obs")
  # Now for each member (without considering run2)
  corr.eccc[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.ECCC,use="pairwise.complete.obs")
  corr.emc[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.EMC,use="pairwise.complete.obs")
  corr.esrl[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.ESRL,use="pairwise.complete.obs")
  corr.nrl[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.NRL,use="pairwise.complete.obs")
  corr.gmao[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.GMAO.run1,use="pairwise.complete.obs")
  corr.rsmas[lead]=cor(dt.verif.SIS[L==lead,]$SIS.index,dt.verif.SIS[L==lead,]$SIS.fcst.RSMAS.run1,use="pairwise.complete.obs")
  
}

dt.corr=data.table(
  L=seq(1,28,1),
  corr.eccc=corr.eccc,
  corr.emc=corr.emc,
  corr.esrl=corr.esrl,
  corr.nrl=corr.nrl,
  corr.gmao=corr.gmao,
  corr.rsmas=corr.rsmas,
  corr.mme=corr.mme
)

colorBlindBlack8  <- c("MME" = "#000000","ECCC" = "#E69F00","EMC" = "#56B4E9","ESRL" = "#009E73", 
                       "NRL" = "#CC79A7","GMAO" ="#0072B2","RSMAS" = "#D55E00")

g2 <- ggplot(data=dt.corr)+
  geom_line(aes(x=L,y=corr.mme,color="MME"),size=1.2)+
  geom_line(aes(x=L,y=corr.eccc,color="ECCC"))+
  geom_line(aes(x=L,y=corr.emc,color="EMC"))+
  geom_line(aes(x=L,y=corr.esrl,color="ESRL"))+
  geom_line(aes(x=L,y=corr.nrl,color="NRL"))+
  geom_line(aes(x=L,y=corr.gmao,color="GMAO"))+
  geom_line(aes(x=L,y=corr.rsmas,color="RSMAS"))+
  labs(y= "COR", x = "Lead time (days)")+
  scale_color_manual(values = colorBlindBlack8)+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0.5),breaks=seq(1,28,2)) + 
  scale_y_continuous(limits=c(0,0.9),breaks=seq(0,1,0.2),minor_breaks =seq(0,1,0.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size=15), axis.title = element_text(size=15))+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15))+
  theme(legend.position=c(0.6,0.90),legend.title=element_blank(),legend.direction = "horizontal")

fig <- grid.arrange(g1,g2, ncol = 2,top = textGrob(paste0("SubX SIS index forecast, inits with obs. SIS<-1 (71 inits)"),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/alvarez/SubX_processed_Rdata/verif_SIS_index_forecast_inits_obsSIS_ltn1.png"),plot=fig,width = 10, height = 5)


