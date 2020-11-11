# SubX data - olr anomaly verification for SubX MME (SubX database)
#
# M. Alvarez - 2020
#-----------------------------------------------------------------------------------------------------------------------
rm(list=ls())

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

#---------------------------------------------------------------------------------------
#  functions of verification metrics 
#---------------------------------------------------------------------------------------

domae <- function(obs, fcst){
  mae=mean(abs(obs-fcst))
  return(mae)
}

dome <- function(obs, fcst){
  me=mean(obs-fcst)
  return(me)
}

dormse <- function(obs, fcst){
  rmse=sqrt(mean((obs-fcst)^2))
  return(rmse)
}

doacc <- function(obs, fcst){
  acc=cor(obs, fcst, use = "pairwise.complete.obs",method = "pearson")
  return(acc)
}

#---------------------------------------------------------------------------------------
#  Main Program  
#---------------------------------------------------------------------------------------
dt.verif = readRDS("/home/alvarez/SubX_processed_Rdata/toverif_MME_ONDEFMA.rds")
dt.verif=as.data.table(dt.verif)

# Comienzo eliminando las filas con week mayor a 4 (sólo había renombrado L en función de la week hasta 35)
dt.verif=dt.verif[week<5,]

# Reduzco hasta el 2014 como startdates para que estén siempre los mismos 6 modelos
dt.verif=dt.verif[year(startdate)<=2014,]

test=dt.verif[,rlutaem.w:=mean(MME,na.rm=TRUE),by=.(lat,lon,startdate,week)] # Funcionó!
test=test[,anom.w:=mean(anom,na.rm=TRUE),by=.(lat,lon,startdate,week)]

test$targetdate=NULL
test$L=NULL
test$rlutaem=NULL
test$anom=NULL

# Ahora elimino las filas duplicadas (porque los 7 días de cada semana quedaron con el mismo valor)
dt.verif.w=unique(test, incomparables=FALSE, fromLast=FALSE)
rm(test)

# Calculo las distintas métricas por cada lon/lat/targetweek
dt.verif.w=dt.verif.w[,me.w:=dome(anom.w,rlutaem.w),by=.(lat,lon,week)]
dt.verif.w=dt.verif.w[,mae.w:=domae(anom.w,rlutaem.w),by=.(lat,lon,week)]
dt.verif.w=dt.verif.w[,rmse.w:=dormse(anom.w,rlutaem.w),by=.(lat,lon,week)]
dt.verif.w=dt.verif.w[,acc.w:=doacc(anom.w,rlutaem.w),by=.(lat,lon,week)]

# Armo dt sólo de scores, eliminando startdates
dt.scores.w=unique(dt.verif.w, incomparables=FALSE, fromLast=FALSE,by=c("lat","lon","week"))
dt.scores.w$startdate=NULL # No me sirve más
dt.scores.w$rlutaem.w=NULL
dt.scores.w$anom.w=NULL

#---------------------------------------------------------------------------------------
#  Gráficos  
#---------------------------------------------------------------------------------------
# Settings
map.world <- map_data("world2")

# Cargo el patrón SIS ONDEFMA para marcar el área sobre cada gráfico


lons=seq(290,327.5,2.5)
lats=seq(5,-40,-2.5)

week.labs <- c("Week 1", "Week 2","Week 3","Week 4")
names(week.labs) <- c("1","2","3","4")

# Cargo el patrón SIS de octubre-abril
sis.oa=metR::ReadNetCDF("weteof.nc",out="array")
sis.oa=sis.oa[[1]]
dimnames(sis.oa)=list(lon=seq(290,327.5,2.5),lat=seq(-40,5,2.5))
melt.sis.oa=melt(sis.oa)
dt.sis.oa=as.data.table(melt.sis.oa)


# Plots
# ACC
g1 <- ggplot() +
  geom_contour_fill(data=dt.scores.w,aes(lon, lat, z = acc.w),breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)) +
  scale_x_longitude(breaks = c(300, 320)) +
  scale_y_latitude(breaks = c(-30,-15,0)) +
  scale_fill_distiller(name="ACC",palette="RdBu",direction=-1,
                       breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),
                       limits = c(-1, 1),
                       guide = guide_colorstrip(),
                       oob  = scales::squish) +
  geom_contour2(data=dt.sis.oa,aes(lon, lat, z = value),alpha=0.5,breaks=seq(3,18,5))+
  geom_contour2(data=dt.sis.oa,aes(lon, lat, z = value),alpha=0.75,breaks=seq(-18,-3,5),linetype="dashed")+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), xlim=range(290,330), ylim=range(-40,5), fill="NA", color="black", inherit.aes = F)+
  theme(axis.text=element_text(size=12))+
  facet_grid(. ~ week,labeller = labeller(week = week.labs))+
  theme(strip.text.x = element_text(size = 12, colour = "black"))+
  theme(strip.background = element_rect(color="black", fill="white", size=1.2, linetype="blank"))+
  coord_cartesian()

# RMSE
g2 <- ggplot(dt.scores.w,aes(lon, lat, z = rmse.w)) +
  geom_contour_fill(breaks = c(10,15,20,25,30,35)) +
  scale_x_longitude(breaks = c(300, 320)) +
  scale_y_latitude(breaks = c(-30,-15,0)) +
  scale_fill_distiller(name="RMSE",palette="YlOrRd",direction=1,
                       breaks = c(10,15,20,25,30,35),
                       limits = c(10,35),
                       guide = guide_colorstrip(),
                       oob  = scales::squish) +
  geom_contour2(data=dt.sis.oa,aes(lon, lat, z = value),alpha=0.5,breaks=seq(3,18,5))+
  geom_contour2(data=dt.sis.oa,aes(lon, lat, z = value),alpha=0.75,breaks=seq(-18,-3,5),linetype="dashed")+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), xlim=range(290,330), ylim=range(-40,5), fill="NA", color="black", inherit.aes = F)+
  theme(axis.text=element_text(size=12))+
  facet_grid(. ~ week,labeller = labeller(week = week.labs))+
  theme(strip.text.x = element_text(size = 12, colour = "black"))+
  theme(strip.background = element_rect(color="black", fill="white", size=1.2, linetype="blank"))+
  coord_cartesian()


# MAE
g3 <- ggplot(dt.scores.w,aes(lon, lat, z = mae.w)) +
  geom_contour_fill(breaks = c(5,10,15,20,25,30)) +
  scale_x_longitude(breaks = c(300, 320)) +
  scale_y_latitude(breaks = c(-30,-15,0)) +
  scale_fill_distiller(name="MAE",palette="YlOrRd",direction=1,
                       breaks = c(5,10,15,20,25,30),
                       limits = c(5,30),
                       guide = guide_colorstrip(),
                       oob  = scales::squish) +
  geom_contour2(data=dt.sis.oa,aes(lon, lat, z = value),alpha=0.5,breaks=seq(3,18,5))+
  geom_contour2(data=dt.sis.oa,aes(lon, lat, z = value),alpha=0.75,breaks=seq(-18,-3,5),linetype="dashed")+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), xlim=range(290,330), ylim=range(-40,5), fill="NA", color="black", inherit.aes = F)+
  theme(axis.text=element_text(size=12))+
  facet_grid(. ~ week,labeller = labeller(week = week.labs))+
  theme(strip.text.x = element_text(size = 12, colour = "black"))+
  theme(strip.background = element_rect(color="black", fill="white", size=1.2, linetype="blank"))+
  coord_cartesian()

# ME
g4 <- ggplot(dt.scores.w,aes(lon, lat, z = me.w)) +
  geom_contour_fill(breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)*1) +
  scale_x_longitude(breaks = c(300, 320)) +
  scale_y_latitude(breaks = c(-30,-15,0)) +
  scale_fill_distiller(name="ME",palette="RdBu",direction=-1,
                       breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)*1,
                       limits = c(-1, 1)*1,
                       guide = guide_colorstrip(),
                       oob  = scales::squish) +
  geom_contour2(data=dt.sis.oa,aes(lon, lat, z = value),alpha=0.5,breaks=seq(3,18,5))+
  geom_contour2(data=dt.sis.oa,aes(lon, lat, z = value),alpha=0.75,breaks=seq(-18,-3,5),linetype="dashed")+
  geom_map(dat=map.world, map = map.world, aes(map_id=region), xlim=range(290,330), ylim=range(-40,5), fill="NA", color="black", inherit.aes = F)+
  theme(axis.text=element_text(size=12))+
  facet_grid(. ~ week,labeller = labeller(week = week.labs))+
  theme(strip.text.x = element_text(size = 12, colour = "black"))+
  theme(strip.background = element_rect(color="black", fill="white", size=1.2, linetype="blank"))+
  coord_cartesian()

fig <- grid.arrange(g1,g2,g4, ncol = 1,top = textGrob("SubX MME rluta (99-14, Oct-Apr)",gp=gpar(fontsize=13,font=3)))
ggsave(filename="/home/alvarez/SubX_processed_Rdata/scores_map_MME.png",plot=fig,width = 10, height = 11)


##########################################################
# Rasters
##########################################################

# Voy a tener que separar en lat/lon de SACZ y lat/lon de SESA. Luego, por cada targetdate-startdate hacer el score para cada par lat/lon.
# Finalmente, promediar el score en el grupo de lat/lon.

dt.sis.oa$center="out"
dt.sis.oa$center[dt.sis.oa$value<=-8]="sesa"
dt.sis.oa$center[dt.sis.oa$value>=13]="sacz"
# ok, funciona bien! Necesito mergear  el dt.sis.oa con dt.verif para luego armar esa columna "sesa"/"sacz"

dt.verif2=merge(dt.verif,dt.sis.oa,by=c("lat","lon"))
dt.verif2=dt.verif2[,!c("week","rlutaem.w","anom.w","value")]

# Para cada startdate, targetdate y center hacer los scores
dt.verif2=dt.verif2[,me:=dome(anom,MME),by=.(startdate,targetdate,center)]
dt.verif2=dt.verif2[,mae:=domae(anom,MME),by=.(startdate,targetdate,center)]
dt.verif2=dt.verif2[,rmse:=dormse(anom,MME),by=.(startdate,targetdate,center)]
dt.verif2=dt.verif2[,acc:=doacc(anom,MME),by=.(startdate,targetdate,center)]

# Para los gráficos, armo la lista de L=1 (sábados)
satsL1=seq(as.Date("1999-01-09"),as.Date("2015-12-26"),14)


# Creo un par de data frames para guardar datos
#df.acc=data.frame(yi=seq(2001,2013,1),ye=seq(2002,2014,1),SACZ_w1=rep(NA_real_,13),SESA_w1=rep(NA_real_,13),SACZ_w2=rep(NA_real_,13),SESA_w2=rep(NA_real_,13),SACZ_w3=rep(NA_real_,13),SESA_w3=rep(NA_real_,13),SACZ_w4=rep(NA_real_,13),SESA_w4=rep(NA_real_,13))
#df.me=data.frame(yi=seq(2001,2013,1),ye=seq(2002,2014,1),SACZ_w1=rep(NA_real_,13),SESA_w1=rep(NA_real_,13),SACZ_w2=rep(NA_real_,13),SESA_w2=rep(NA_real_,13),SACZ_w3=rep(NA_real_,13),SESA_w3=rep(NA_real_,13),SACZ_w4=rep(NA_real_,13),SESA_w4=rep(NA_real_,13))

dt.acc=data.table(yi=seq(2001,2013,1))
dt.acc = dt.acc[, .(week = 1:4), by = c(colnames(dt.acc))]
dt.acc = dt.acc[, .(region = c("SACZ","SESA")), by = c(colnames(dt.acc))]
dt.acc$ACC=NA_real_
dt.acc$ME=NA_real_


# Acá debería empezar el for en los años
for(yyi in (1999:2013)){

  yye=yyi+1

  # A los gráficos les voy a agregar el promedio de los scores por semana.
  scores = dt.verif2[(((year(startdate)==yyi & month(startdate)>=10) | (year(startdate)==yye & month(startdate)<=4))),]

  accw1_sacz = mean(scores[center=="sacz" & L>=1 & L<=7]$acc)
  accw2_sacz = mean(scores[center=="sacz" & L>=8 & L<=14]$acc)
  accw3_sacz = mean(scores[center=="sacz" & L>=15 & L<=21]$acc)
  accw4_sacz = mean(scores[center=="sacz" & L>=23 & L<=28]$acc)

  accw1_sesa = mean(scores[center=="sesa" & L>=1 & L<=7]$acc)
  accw2_sesa = mean(scores[center=="sesa" & L>=8 & L<=14]$acc)
  accw3_sesa = mean(scores[center=="sesa" & L>=15 & L<=21]$acc)
  accw4_sesa = mean(scores[center=="sesa" & L>=23 & L<=28]$acc)

  # Quiero ir guardando estos datos
  # df.acc$SACZ_w1[df.acc$yi==yyi]=accw1_sacz
  # df.acc$SACZ_w2[df.acc$yi==yyi]=accw2_sacz
  # df.acc$SACZ_w3[df.acc$yi==yyi]=accw3_sacz
  # df.acc$SACZ_w4[df.acc$yi==yyi]=accw4_sacz
  #
  # df.acc$SESA_w1[df.acc$yi==yyi]=accw1_sesa
  # df.acc$SESA_w2[df.acc$yi==yyi]=accw2_sesa
  # df.acc$SESA_w3[df.acc$yi==yyi]=accw3_sesa
  # df.acc$SESA_w4[df.acc$yi==yyi]=accw4_sesa

  dt.acc$ACC[dt.acc$yi==yyi & dt.acc$week==1 & dt.acc$region=="SACZ"]=accw1_sacz
  dt.acc$ACC[dt.acc$yi==yyi & dt.acc$week==2 & dt.acc$region=="SACZ"]=accw2_sacz
  dt.acc$ACC[dt.acc$yi==yyi & dt.acc$week==3 & dt.acc$region=="SACZ"]=accw3_sacz
  dt.acc$ACC[dt.acc$yi==yyi & dt.acc$week==4 & dt.acc$region=="SACZ"]=accw4_sacz

  dt.acc$ACC[dt.acc$yi==yyi & dt.acc$week==1 & dt.acc$region=="SESA"]=accw1_sesa
  dt.acc$ACC[dt.acc$yi==yyi & dt.acc$week==2 & dt.acc$region=="SESA"]=accw2_sesa
  dt.acc$ACC[dt.acc$yi==yyi & dt.acc$week==3 & dt.acc$region=="SESA"]=accw3_sesa
  dt.acc$ACC[dt.acc$yi==yyi & dt.acc$week==4 & dt.acc$region=="SESA"]=accw4_sesa



  g5 <- ggplot(dt.verif2[(((year(startdate)==yyi & month(startdate)>=10) | (year(startdate)==yye & month(startdate)<=4)) & center=="sacz"),], aes(targetdate, L)) +
  geom_raster(aes(fill = cut(acc,c(-1,-0.75,-0.6,-0.45,-0.3,-0.15,0,0.15,0.3,0.45,0.6,0.75,1)))) +
  scale_fill_brewer("ACC",palette = "RdBu",direction=-1) +
  scale_x_date(name="Target date",breaks = satsL1[((year(satsL1)==yyi & month(satsL1)>=10) | (year(satsL1)==yye & month(satsL1)<=4))],minor_breaks = waiver(), date_labels = "%d%b", expand = expand_scale(mult = 0,add = 0))+
  scale_y_continuous(name="Lead (days)",breaks=(seq(1,28,2)),limits=c(0.5,28.5),expand=c(0,0))+
  labs(title=paste0("SACZ region - Oct-Apr ",as.character(yyi),"-",as.character(yye)),tag="a)")+
  theme_bw()+  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  geom_hline(yintercept = 7.5,size=0.25)+
  geom_hline(yintercept = 14.5,size=0.25)+
  geom_hline(yintercept = 21.5,size=0.25)+
  annotate("text",x=as.Date(paste0(as.character(yyi),"-10-13")), y=18, label=as.character(round(accw3_sacz,digits=2)),size=4)+
  annotate("text",x=as.Date(paste0(as.character(yyi),"-10-13")), y=25, label=as.character(round(accw4_sacz,digits=2)),size=4)+
  annotate("text",x=as.Date(paste0(as.character(yye),"-05-15")), y=4, label=as.character(round(accw1_sacz,digits=2)),size=4)+
  annotate("text",x=as.Date(paste0(as.character(yye),"-05-15")), y=11, label=as.character(round(accw2_sacz,digits=2)),size=4)


  g6 <-ggplot(dt.verif2[(((year(startdate)==yyi & month(startdate)>=10) | (year(startdate)==yye & month(startdate)<=4)) & center=="sesa"),], aes(targetdate, L)) +
  geom_raster(aes(fill = cut(acc,c(-1,-0.75,-0.6,-0.45,-0.3,-0.15,0,0.15,0.3,0.45,0.6,0.75,1)))) +
  scale_fill_brewer("ACC",palette = "RdBu",direction=-1) +
  scale_x_date(name="Target date",breaks = satsL1[((year(satsL1)==yyi & month(satsL1)>=10) | (year(satsL1)==yye & month(satsL1)<=4))],minor_breaks = waiver(), date_labels = "%d%b", expand = expand_scale(mult = 0,add = 0))+
  scale_y_continuous(name="Lead (days)",breaks=(seq(1,28,2)),limits=c(0.5,28.5),expand=c(0,0))+
  labs(title=paste0("SESA region - Oct-Apr ",as.character(yyi),"-",as.character(yye)),tag="b)")+
  theme_bw()+  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  geom_hline(yintercept = 7.5,size=0.25)+
  geom_hline(yintercept = 14.5,size=0.25)+
  geom_hline(yintercept = 21.5,size=0.25)+
  annotate("text",x=as.Date(paste0(as.character(yyi),"-10-13")), y=18, label=as.character(round(accw3_sesa,digits=2)),size=4)+
  annotate("text",x=as.Date(paste0(as.character(yyi),"-10-13")), y=25, label=as.character(round(accw4_sesa,digits=2)),size=4)+
  annotate("text",x=as.Date(paste0(as.character(yye),"-05-15")), y=4, label=as.character(round(accw1_sesa,digits=2)),size=4)+
  annotate("text",x=as.Date(paste0(as.character(yye),"-05-15")), y=11, label=as.character(round(accw2_sesa,digits=2)),size=4)


fig <- grid.arrange(g5,g6, ncol = 1,top = textGrob("SubX MME rluta",gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/alvarez/SubX_processed_Rdata/ACC_regionxtarget_MME_",as.character(yyi),"-",as.character(yye),".png"),plot=fig,width = 14, height = 10)


}

# Vamos con el Bias

# Acá debería empezar el for en los años
for(yyi in (1999:2013)){

  # A los gráficos les voy a agregar el promedio de los scores por semana.
  scores = dt.verif2[(((year(startdate)==yyi & month(startdate)>=10) | (year(startdate)==yye & month(startdate)<=4))),]

  mew1_sacz = mean(scores[center=="sacz" & L>=1 & L<=7]$me)
  mew2_sacz = mean(scores[center=="sacz" & L>=8 & L<=14]$me)
  mew3_sacz = mean(scores[center=="sacz" & L>=15 & L<=21]$me)
  mew4_sacz = mean(scores[center=="sacz" & L>=23 & L<=28]$me)

  mew1_sesa = mean(scores[center=="sesa" & L>=1 & L<=7]$me)
  mew2_sesa = mean(scores[center=="sesa" & L>=8 & L<=14]$me)
  mew3_sesa = mean(scores[center=="sesa" & L>=15 & L<=21]$me)
  mew4_sesa = mean(scores[center=="sesa" & L>=23 & L<=28]$me)



  # Quiero ir guardando estos datos
  # df.me$SACZ_w1[df.me$yi==yyi]=mew1_sacz
  # df.me$SACZ_w2[df.me$yi==yyi]=mew2_sacz
  # df.me$SACZ_w3[df.me$yi==yyi]=mew3_sacz
  # df.me$SACZ_w4[df.me$yi==yyi]=mew4_sacz
  #
  # df.me$SESA_w1[df.me$yi==yyi]=mew1_sesa
  # df.me$SESA_w2[df.me$yi==yyi]=mew2_sesa
  # df.me$SESA_w3[df.me$yi==yyi]=mew3_sesa
  # df.me$SESA_w4[df.me$yi==yyi]=mew4_sesa


  dt.acc$ME[dt.acc$yi==yyi & dt.acc$week==1 & dt.acc$region=="SACZ"]=mew1_sacz
  dt.acc$ME[dt.acc$yi==yyi & dt.acc$week==2 & dt.acc$region=="SACZ"]=mew2_sacz
  dt.acc$ME[dt.acc$yi==yyi & dt.acc$week==3 & dt.acc$region=="SACZ"]=mew3_sacz
  dt.acc$ME[dt.acc$yi==yyi & dt.acc$week==4 & dt.acc$region=="SACZ"]=mew4_sacz

  dt.acc$ME[dt.acc$yi==yyi & dt.acc$week==1 & dt.acc$region=="SESA"]=mew1_sesa
  dt.acc$ME[dt.acc$yi==yyi & dt.acc$week==2 & dt.acc$region=="SESA"]=mew2_sesa
  dt.acc$ME[dt.acc$yi==yyi & dt.acc$week==3 & dt.acc$region=="SESA"]=mew3_sesa
  dt.acc$ME[dt.acc$yi==yyi & dt.acc$week==4 & dt.acc$region=="SESA"]=mew4_sesa

  yye=yyi+1

  g5 <- ggplot(dt.verif2[(((year(startdate)==yyi & month(startdate)>=10) | (year(startdate)==yye & month(startdate)<=4)) & center=="sacz"),], aes(targetdate, L)) +
    geom_raster(aes(fill = cut(me,(40*c(-1,-0.75,-0.6,-0.45,-0.3,-0.15,0.15,0.3,0.45,0.6,0.75,1))))) +
    scale_fill_brewer("ME",palette = "RdBu",direction=-1) +
    scale_x_date(name="Target date",breaks = satsL1[((year(satsL1)==yyi & month(satsL1)>=10) | (year(satsL1)==yye & month(satsL1)<=4))],minor_breaks = waiver(), date_labels = "%d%b", expand = expand_scale(mult = 0,add = 0))+
    scale_y_continuous(name="Lead (days)",breaks=(seq(1,28,2)),limits=c(0.5,28.5),expand=c(0,0))+
    labs(title=paste0("SACZ region - Oct-Apr ",as.character(yyi),"-",as.character(yye)),tag="a)")+
    theme_bw()+  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
    geom_hline(yintercept = 7.5,size=0.25)+
    geom_hline(yintercept = 14.5,size=0.25)+
    geom_hline(yintercept = 21.5,size=0.25)+
    annotate("text",x=as.Date(paste0(as.character(yyi),"-10-13")), y=18, label=as.character(round(mew3_sacz,digits=2)),size=4)+
    annotate("text",x=as.Date(paste0(as.character(yyi),"-10-13")), y=25, label=as.character(round(mew4_sacz,digits=2)),size=4)+
    annotate("text",x=as.Date(paste0(as.character(yye),"-05-15")), y=4, label=as.character(round(mew1_sacz,digits=2)),size=4)+
    annotate("text",x=as.Date(paste0(as.character(yye),"-05-15")), y=11, label=as.character(round(mew2_sacz,digits=2)),size=4)


  g6 <-ggplot(dt.verif2[(((year(startdate)==yyi & month(startdate)>=10) | (year(startdate)==yye & month(startdate)<=4)) & center=="sesa"),], aes(targetdate, L)) +
    geom_raster(aes(fill = cut(me,(40*c(-1,-0.75,-0.6,-0.45,-0.3,-0.15,0.15,0.3,0.45,0.6,0.75,1))))) +
    scale_fill_brewer("ME",palette = "RdBu",direction=-1) +
    scale_x_date(name="Target date",breaks = satsL1[((year(satsL1)==yyi & month(satsL1)>=10) | (year(satsL1)==yye & month(satsL1)<=4))],minor_breaks = waiver(), date_labels = "%d%b", expand = expand_scale(mult = 0,add = 0))+
    scale_y_continuous(name="Lead (days)",breaks=(seq(1,28,2)),limits=c(0.5,28.5),expand=c(0,0))+
    labs(title=paste0("SESA region - Oct-Apr ",as.character(yyi),"-",as.character(yye)),tag="b)")+
    theme_bw()+  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
    geom_hline(yintercept = 7.5,size=0.25)+
    geom_hline(yintercept = 14.5,size=0.25)+
    geom_hline(yintercept = 21.5,size=0.25)+
    annotate("text",x=as.Date(paste0(as.character(yyi),"-10-13")), y=18, label=as.character(round(mew3_sesa,digits=2)),size=4)+
    annotate("text",x=as.Date(paste0(as.character(yyi),"-10-13")), y=25, label=as.character(round(mew4_sesa,digits=2)),size=4)+
    annotate("text",x=as.Date(paste0(as.character(yye),"-05-15")), y=4, label=as.character(round(mew1_sesa,digits=2)),size=4)+
    annotate("text",x=as.Date(paste0(as.character(yye),"-05-15")), y=11, label=as.character(round(mew2_sesa,digits=2)),size=4)


  fig <- grid.arrange(g5,g6, ncol = 1,top = textGrob("SubX MME rluta",gp=gpar(fontsize=13,font=3)))
  ggsave(filename=paste0("/home/alvarez/SubX_processed_Rdata/MeanError_regionxtarget_MME_",as.character(yyi),"-",as.character(yye),".png"),plot=fig,width = 14, height = 10)

}



# # Armo los box plot con los resultados de los scores

library("ggsci")


dt.acc$week=as.character(dt.acc$week)
dt.acc$combi=interaction(dt.acc$region,dt.acc$week)

# ACC

accsacz=dt.acc[region=="SACZ",c("ACC")]
accsesa=dt.acc[region=="SESA",c("ACC")]
accweek=dt.acc[region=="SESA",c("week")]

scatacc=data.frame(accsacz,accsesa,accweek)

(scat <-
ggplot(scatacc,aes(x=ACC,y=ACC.1,color=week))+
  geom_point()+
  scale_y_continuous(limits=c(-0.025,0.3))+
  scale_x_continuous(limits=c(-0.025,0.3))+
  theme(line = element_blank(),
        panel.background = element_rect(fill = "transparent"))+
    labs(x="ACC SACZ",y="ACC SESA")+
  geom_segment(x=-0.025,y=-0.025,xend=0.3,yend=0.3,color="gray60",alpha=0.1)

)

boxacc <- ggplot(dt.acc, aes(x = combi, y = ACC,color=region)) +
  coord_flip() +
  scale_y_continuous(limits = c(-0.05, 0.35))+
  labs(x = element_blank(), y = "ACC") +
  geom_boxplot(color = "gray60", outlier.alpha = 0) +
  geom_jitter(size = 2, alpha = 0.35, width = 0.2) +
  scale_color_manual(values=c('#8c510a','#01665e'))+
  annotation_custom(ggplotGrob(scat), xmin = 4, xmax = 8.5, ymin = 0.2, ymax = 0.365)+
  theme_bw()+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

ggsave(filename=paste0("/home/alvarez/SubX_processed_Rdata/ACC_centers_boxplot.png"),plot=boxacc,width = 14, height = 7)


# ME

mesacz=dt.acc[region=="SACZ",c("ME")]
mesesa=dt.acc[region=="SESA",c("ME")]
meweek=dt.acc[region=="SESA",c("week")]

scatme=data.frame(mesacz,mesesa,meweek)

(scat <-
    ggplot(scatme,aes(x=ME,y=ME.1,color=week))+
    geom_point()+
    scale_y_continuous(limits=c(-4.25,2.25))+
    scale_x_continuous(limits=c(-4.25,2.25))+
    theme(line = element_blank(),
          panel.background = element_rect(fill = "transparent"))+
    labs(x="ME SACZ",y="ME SESA")+
    geom_segment(x=0,y=-4.25,xend=0,yend=2.25,color="gray60",alpha=0.1)+
    geom_segment(x=-4.25,y=0,xend=2.25,yend=0,color="gray60",alpha=0.1)

)

boxme <- ggplot(dt.acc, aes(x = combi, y = ME,color=region)) +
  coord_flip() +
  scale_y_continuous(limits = c(-4.25, 7))+
  labs(x = element_blank(), y = "ME") +
  geom_boxplot(color = "gray60", outlier.alpha = 0) +
  geom_jitter(size = 2, alpha = 0.35, width = 0.2) +
  scale_color_manual(values=c('#8c510a','#01665e'))+
  annotation_custom(ggplotGrob(scat), xmin = 4, xmax = 8.5, ymin = 2.75, ymax = 7.5)+
  theme_bw()+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12))

ggsave(filename=paste0("/home/alvarez/SubX_processed_Rdata/ME_centers_boxplot.png"),plot=boxme,width = 14, height = 7)
