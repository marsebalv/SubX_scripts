# SubX data - olr anomaly verification for GMAO model (SubX database)
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
dt.verif = readRDS("/home/alvarez/SubX_processed_Rdata/toverif_GMAO_ONDEFMA.rds")
dt.verif=as.data.table(dt.verif)

# Comienzo eliminando las filas con week mayor a 4 (sólo había renombrado L en función de la week hasta 35)
dt.verif=dt.verif[week<5,]

test=dt.verif[,rlutaem.w:=mean(rlutaem,na.rm=TRUE),by=.(lat,lon,startdate,week)] # Funcionó!
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

fig <- grid.arrange(g1,g2,g4, ncol = 1,top = textGrob("SubX GMAO-GEOS_V2p1 rluta (99-15, Oct-Apr)",gp=gpar(fontsize=13,font=3)))
ggsave(filename="/home/alvarez/SubX_processed_Rdata/scores_map_GMAO.png",plot=fig,width = 10, height = 11)


