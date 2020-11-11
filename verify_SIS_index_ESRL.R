# SubX data - SIS index verification for ESRL model (SubX database)
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
dt.verif = as.data.table(readRDS("/home/alvarez/SubX_processed_Rdata/toverif_ESRL_ONDEFMA.rds"))

# Comienzo eliminando las filas con week mayor a 4 (sólo había renombrado L en función de la week hasta 35)
# Como ahora voy a tener promedios móviles de 5 días, voy a necesitar hasta 4 días más de lead (hasta 32)
#dt.verif=dt.verif[week<5,]

# Reduzco hasta el 2014 como startdates para que estén siempre los mismos 6 modelos
dt.verif=dt.verif[year(startdate)<=2014,]
dt.verif$week=NULL
dt.verif$anom=NULL

#***********************************************************
# Ahora lo que necesito es alinear las startdates con las del MME
# Create date vector of Saturdays in 1999-2015
fridays=seq(as.Date("1999-01-08"),as.Date("2015-12-25"),7)
prevsats=seq(as.Date("1999-01-02"),as.Date("2015-12-19"),7)

dt.verif$gatherweek=NA_real_
dt.verif$mmestartdate=as.Date("1800-01-01")
for(i in 1:886){
  dt.verif$gatherweek[(dt.verif$startdate>=prevsats[i] & dt.verif$startdate<=fridays[i])]=i
  dt.verif$mmestartdate[(dt.verif$startdate>=prevsats[i] & dt.verif$startdate<=fridays[i])]=fridays[i]
}

dt.verif$L.MME=NA_real_
dt.verif$L.MME=as.numeric(dt.verif$targetdate-dt.verif$mmestartdate)

# Ahora L startdate y gatherweek no me sivern más. Me debería quedar además sólo con las L.MME >= 1
dt.verif=dt.verif[(L.MME>0 & L.MME<50),!c("L","startdate","gatherweek")]

#***********************************************************

# Cargo el RTSIS para la verificación
load("~/SubX_scripts/RTSIS_OA9899_1415.rda")

# Cargo las anomalías de OLR suavizadas en RT-IS
load("~/SubX_scripts/lfremoved_RT_rluta.rda")
dt.olra$anom=NULL

# dt.verif: lon, lat, targetdate, startdate, L, MME
# dt.olra: lon, lat, targetdate, smt2.olra

# Ahora voy a a buscar agregar a la data table 4 valores por cada targetdate. Creo que me conviene sacar L y después volver a agregarlo
dt.verif$L.MME=NULL
# Renombro lf.olra a rlutaem aunque no sea modelo
setnames(dt.olra,old="lf.olra",new="rlutaem")
setnames(dt.verif,old="mmestartdate",new="startdate")

# Debería agregar una fila con lon, lat, targetdate, startdate y MME

#1. Necesito un vector con todas las startdates
starts=unique(dt.verif$startdate,incomparables = FALSE)

# Armo una data table con las observaciones para unir (Gracias Elio!)

obs_append <- list()
for (i in 1:length(starts)) {
  date = starts[i]
  dates <- seq(date - 3, date ,1) #necesito la del startdate y las 3 anteriores
  i <- length(obs_append) + 1
  obs <- dt.olra[targetdate %in% dates]
  obs <- obs[, startdate := date]
  # acá agregarle las columnas necesarias (targetdate, startdate..)
  obs_append[[i]] <- obs
} 
obs_append <- rbindlist(obs_append)
setcolorder(obs_append,neworder=c("lat","lon","targetdate","startdate","rlutaem"))

obs_append = data.table(obs_append)
# Uno a dt.verif

dt.verif = rbind(dt.verif,obs_append)

# Agrego la variable L
dt.verif[,L := targetdate - startdate]

#*********************************
# Ahora debería cargar el patrón SIS y combinar los datos
load("~/SubX_scripts/SISpattern_OA.rda")

# Esta función aplica un promedio móvil sin peso centrado de 5tp orden
smoothHF <- function(data){
  times = length(data)
  filterhf = rep(NA_real_,times)
  for (t in 3:(times-2)){
    aux = data[(t-2):(t+2)]
    filterhf[t] = mean(aux,na.rm=FALSE)
  }
  return(filterhf)
}

# Antes de aplicar la función me tengo que asegurar que la dt esté ordenada cronológicamente por lon y lat
dt.verif = dt.verif[order(targetdate),]

# Suavizo alta frecuencia
dt.verif = dt.verif[,smt1.olra := smoothHF(rlutaem),by=.(lon,lat,startdate)]
dt.verif = dt.verif[,smt2.olra := smoothHF(smt1.olra),by=.(lon,lat,startdate)]

# Bien, ahora descarto lo que no me sirve
dt.verif$rlutaem=NULL
dt.verif$smt1.olra=NULL
setnames(dt.verif,old="smt2.olra",new="EM")

rm(dt.olra,obs,obs_append)

# Elimino los L que no son pronóstico (pero retengo el 0 que me servirá para empalmar)
dt.verif = dt.verif[L>=0,]

dt.verif=merge(dt.verif,SIS.pattern.destandardized)

# Renombro columnas de la data table
setnames(dt.verif,old="value",new="stdSISpattern")

# Me quedo sólo con las startdates de ONDEFMA
dt.verif$startmonth=month(dt.verif$startdate)
OA = c(1,2,3,4,10,11,12);

dt.verif=dt.verif[dt.verif$startmonth %in% OA,]
dt.verif$startmonth=NULL

# Calculo el índice SIS para ONDEFMA con metodología RT
dt.verif = dt.verif[,SIS.index.fcst := (((sum(EM*stdSISpattern))*SISscale.OA-SISmean.OA)/SISstd.OA),by=.(targetdate,startdate)]

# Armo una dt para el SIS
dt.fcstSIS=dt.verif[lat==-40 & lon==290,.(targetdate,startdate,L,SIS.index.fcst)]

# Para el empalme en el gráfico
dt.union = dt.fcstSIS[L<2,]

# Listo, ya tengo calculado el pronóstico del SIS!!!
unionsttdates=unique(dt.union$startdate,incomparables = FALSE)
for(j in 1:length(unionsttdates)){
  dt.union[L==0 & startdate==unionsttdates[j]]$SIS.index.fcst=dt.RTSIS[targetdate == unionsttdates[j]]$SIS.index
}


# Voy a eliminar los leads mayores a 28
dt.fcstSIS=dt.fcstSIS[L<29,]

# Guardo para luego unir con los otros modelos y generar el spread

save(dt.fcstSIS, file = "~/SubX_scripts/SubX_SIS_index_forecast_ESRL.rda")
#--------------------------------------------------------------------
# Gráfico

# Selecciono la temporada a graficar
for(yi in 1999:2013){
#yi=2009

rtsis=dt.RTSIS[((month(targetdate) %in% c(10,11,12)) & year(targetdate)==yi) | ((month(targetdate) %in% c(1,2,3,4)) & year(targetdate)==(yi+1)) ,]
fcsis=dt.fcstSIS[((month(targetdate) %in% c(10,11,12)) & year(targetdate)==yi) | ((month(targetdate) %in% c(1,2,3,4)) & year(targetdate)==(yi+1)),]

# Determino las startdates
sttdates=unique(fcsis$startdate,incomparables = FALSE)
# Busco valores del RTSIS (obs) para las fechas de inicio (para marcarlos en el gráfico)
sttrtsis=rep(NA_real_,length(sttdates))
for(i in 1:length(sttdates)){
  sttrtsis[i]=rtsis[targetdate == sttdates[i],SIS.index]
}
sttvals=data.frame(sttdates,sttrtsis)

thecolors=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#cab2d6","#a65628","#f781bf")

#Arranco desde la primera startdate
s=1 
ss=seq(s,length(sttdates),4)

g1 <- ggplot()+
  theme_bw()+
  geom_line(data=rtsis,aes(targetdate,SIS.index),alpha=0.5)+
  scale_x_date(breaks=sttdates[seq(1,length(sttdates),2)],date_labels = "%d-%b",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1.),limits=c(-3.3,3.3),expand = c(0., 0.)) +
  geom_point(data=sttvals,aes(sttdates,sttrtsis),shape=18,color="#f781bf",size=3)

for(f in 1:length(ss)){
  g1 <- g1 + geom_line(data=fcsis[startdate == sttdates[ss[f]]],aes(x=targetdate,y=SIS.index.fcst),color=thecolors[f])
  g1 <- g1 + geom_line(data=dt.union[startdate == sttdates[ss[f]]],aes(x=targetdate,y=SIS.index.fcst),color=thecolors[f],linetype="dashed")
}

#Segunda
s=2
ss=seq(s,length(sttdates),4)

g2 <- ggplot()+
  theme_bw()+
  geom_line(data=rtsis,aes(targetdate,SIS.index),alpha=0.5)+
  scale_x_date(breaks=sttdates[seq(1,length(sttdates),2)],date_labels = "%d-%b",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1.),limits=c(-3.3,3.3),expand = c(0., 0.)) +
  geom_point(data=sttvals,aes(sttdates,sttrtsis),shape=18,color="#f781bf",size=3)

for(f in 1:length(ss)){
  g2 <- g2 + geom_line(data=fcsis[startdate == sttdates[ss[f]]],aes(x=targetdate,y=SIS.index.fcst),color=thecolors[f])
  g2 <- g2 + geom_line(data=dt.union[startdate == sttdates[ss[f]]],aes(x=targetdate,y=SIS.index.fcst),color=thecolors[f],linetype="dashed")
}


#Tercera
s=3
ss=seq(s,length(sttdates),4)

g3 <- ggplot()+
  theme_bw()+
  geom_line(data=rtsis,aes(targetdate,SIS.index),alpha=0.5)+
  scale_x_date(breaks=sttdates[seq(1,length(sttdates),2)],date_labels = "%d-%b",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1.),limits=c(-3.3,3.3),expand = c(0., 0.)) +
  geom_point(data=sttvals,aes(sttdates,sttrtsis),shape=18,color="#f781bf",size=3)

for(f in 1:length(ss)){
  g3 <- g3 + geom_line(data=fcsis[startdate == sttdates[ss[f]]],aes(x=targetdate,y=SIS.index.fcst),color=thecolors[f])
  g3 <- g3 + geom_line(data=dt.union[startdate == sttdates[ss[f]]],aes(x=targetdate,y=SIS.index.fcst),color=thecolors[f],linetype="dashed")
}

#Segunda
s=4
ss=seq(s,length(sttdates),4)

g4 <- ggplot()+
  theme_bw()+
  geom_line(data=rtsis,aes(targetdate,SIS.index),alpha=0.5)+
  scale_x_date(breaks=sttdates[seq(1,length(sttdates),2)],date_labels = "%d-%b",expand = c(0, 0))+
  theme(text = element_text(size=14))+
  scale_y_continuous(breaks = seq(-3,3,1.),limits=c(-3.3,3.3),expand = c(0., 0.)) +
  geom_point(data=sttvals,aes(sttdates,sttrtsis),shape=18,color="#f781bf",size=3)

for(f in 1:length(ss)){
  g4 <- g4 + geom_line(data=fcsis[startdate == sttdates[ss[f]]],aes(x=targetdate,y=SIS.index.fcst),color=thecolors[f])
  g4 <- g4 + geom_line(data=dt.union[startdate == sttdates[ss[f]]],aes(x=targetdate,y=SIS.index.fcst),color=thecolors[f],linetype="dashed")
}

fig <- grid.arrange(g1,g2,g3,g4, ncol = 1,top = textGrob(paste0("SubX ESRL EM SIS index forecast season ",as.character(yi),"/",as.character(yi+1)),gp=gpar(fontsize=13,font=3)))
ggsave(filename=paste0("/home/alvarez/SubX_processed_Rdata/SIS_index_forecast_ESRL_",as.character(yi),"_",as.character(yi+1),".png"),plot=fig,width = 11, height = 14)

}


