# SubX download data
#
# M. Alvarez - 2020
#-----------------------------------------------------------------------------------------------------------------------
rm(list=ls())

setwd("/home/alvarez/")

# Call libraries to be used
library("ncdf4")
library("metR")
library('pracma')

#--------------------------------------------------------------------------------------
#  Variables to be modified by user
#--------------------------------------------------------------------------------------
  
#outPath='/pikachu/datos4/SubX/'
outPath='/datos/SubX/'
#varnames=c('ua','ua','rlut','tas','ts','zg','va','va','pr','zg')          # Variable names
#plevstrs=c('850','200','toa','2m','sfc','500','200','850','sfc','200')    # Must be same size as varname, for variables with no plevels, use sfc or 10m
groups=c('GMAO','RSMAS','ESRL','ECCC','NRL','EMC')                        # Modeling Groups (must be same # elements as models below)
models=c('GEOS_V2p1','CCSM4','FIMr1p1','GEM','NESM','GEFS')               # Model Name (must be same # of elements as groups above)
# groups=c('RSMAS','ESRL','ECCC','NRL','EMC')                        # Modeling Groups (must be same # elements as models below)
# models=c('CCSM4','FIMr1p1','GEM','NESM','GEFS')               # Model Name (must be same # of elements as groups above)

# groups=c('GMAO')
# models=c('GEOS_V2p1')
varnames=c('rlut')
plevstrs=c('toa')
                                                                                                                                                                
dfv=-9.99e-8;                                                             # Default missing_value or FillValue if not specified in input data file
type='hindcast'                                                           # hindcast or forecast
                
#--------------------------------------------------------------------------------------
# DO NOT MODIFY
#--------------------------------------------------------------------------------------
url='http://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/' # IRI URL
nvars=numel(varnames)
nmodels=numel(models)      

#---------------------------------------------------------------------------------------
#  Main Program  
#---------------------------------------------------------------------------------------
 
for(ivar in 1:nvars){
  # Define variable name and level
  varname=varnames[ivar]
  plevstr=plevstrs[ivar]
  
  for(imodel in 1:nmodels){
    # Define the model and group
    model=models[imodel]
    group=groups[imodel]
    
    fprintf('%s%s%s%s%s%s\n','Getting: ',group,'-',model,' ',varname)
    
    # Open File
    inFname=paste0(url,'.',group,'/.',model,'/.',type,'/.',varname,'/dods') # input filename
    ncid=nc_open(inFname, write=FALSE, readunlim=TRUE, verbose=FALSE,auto_GMT=TRUE, suppress_dimvals=FALSE )
    
    # Determine number of variables, dimensions, etc. in file
    ndims = ncid$ndims
    nvars = ncid$nvars
    ngatts = ncid$natts
    unlimdimid = ncid$unlimdimid
    
    # Determine dimensions and size of dimensions
    # (EVITO EL FOR) PUEDE SER QUE ME FALTE ALGUNA VARIABLE ACÁ
    lon=ncid$dim$X$vals
    lat=ncid$dim$Y$vals
    leads=ncid$dim$L$vals
    ens=ncid$dim$M$vals
    nens=ncid$dim$M$len
    ics=ncid$dim$S$vals
    unitsic=ncid$dim$S$units 
    nics=ncid$dim$S$len
    if(ndims==6){
    # There are pressure levels
      levs=ncid$dim$P$vals
      plev=which(levs==as.numeric(plevstr))
      }
    
    # Convert netCDF days since values to calendar dates
    origin=substr(ncid$dim$S$units,12,21)
    startdates=as.character(as.Date(ncid$dim$S$vals,origin=as.Date(origin)))

    # Get information about specified variable
    # Set the dimensions to read - [nlons,nlats,nleads,nens,nlevs,nics]
    vardimlen=ncid$var$rlut$size   # EL TEMA ACÁ ES QUE ME QUEDA EL NOMBRE DE LA VARIABLE EN EL MEDIO DE VARDIMLEN
    
    #if(ndims == 6){  # multilevel data
    #count=c(vardimlen[1],vardimlen[2],vardimlen[3],1,1,1)
    #}else{ # single level data
    #count=c(vardimlen[1],vardimlen[2],vardimlen[3],1,1)
    #}
    
    # Create Output directory if needed
    outDir=paste0(outPath,type,'/',varname,plevstr,'/daily/full/',group,'-',model,'/');
    dir.create(outDir,recursive=TRUE)
    
    # COrregir pues rlut en el medio
    # Set Global Attribute Values
    units=ncid$var$rlut$units       
    longtitle=ncid$var$rlut$longname
    title=longtitle
    fillValue=ncid$var$rlut$missval #No lo encontré así que pongo este por ahora. Habrá que buscar
    comments='SubX project http://cola.gmu.edu/~kpegion/subx/'
    source='SubX IRI'
    institution='IRI'
    
    # Ahora empezaría el loop para guardar. Todavia me 
    # falta corregir lo de los tiempos y todos los lados donde se intercala rlut como nombre
    
    # Loop over all ensemble members
    for(iens in 0:(nens-1)){
      
      # Set the ensemble member string for output file
      ee=as.character(as.integer(iens+1))
      
      for(i in 0:(nics-1)){
        
        # Loop over and read all start dates  REVISAR
        imn=as.numeric(substr(startdates[i+1],6,7))
        iyr=as.numeric(substr(startdates[i+1],1,4))
        idy=as.numeric(substr(startdates[i+1],9,10))
        
        # Construct output filename
        ofname=paste0(outDir,varname,'_',as.character(plevstr),'_',group,'-',model,'_',substr(startdates[i+1],1,4),substr(startdates[i+1],6,7),substr(startdates[i+1],9,10),'.e',ee,'.SISdom.daily.nc')

        # Read all leads, lats, lons for a single start date, ensemble member, and pressure level
        # Resulting array is of size data[nx,ny,nlead] 
        # Define subdomain
        subdomX=seq(287,330,1)
        subdomY=seq(-43,8,1)
        
        if(ndims == 6){ # multilevel data
          data = metR::ReadNetCDF(inFname, subset = list(M = (iens+1), P = plev, S = as.Date(startdates[i+1]), X = subdomX, Y=subdomY),out='array')
        } else{ # single level data
          data = metR::ReadNetCDF(inFname, subset = list(M = (iens+1), S = as.Date(startdates[i+1]),X = subdomX, Y=subdomY),out='array')
        }
        
        data=data[[1]]
        data=drop(data)
        
        fprintf('%s%s\n','Writing File: ',ofname)
        
        # Write dat
        #1. Define dimensions
        #londim <- ncdim_def("X",ncid$dim$X$units,as.double(ncid$dim$X$vals)) 
        #latdim <- ncdim_def("Y",ncid$dim$Y$units,as.double(ncid$dim$Y$vals)) 
        londim <- ncdim_def("X",ncid$dim$X$units,as.double(subdomX)) 
        latdim <- ncdim_def("Y",ncid$dim$Y$units,as.double(subdomY))       
        timedim <- ncdim_def("L",ncid$dim$L$units,as.double(ncid$dim$L$vals))
        # 2. Define variables
        var_def <- ncvar_def(varname,ncid$var$rlut$units,list(londim,latdim,timedim),fillValue,varname,prec="float")
        # Create netCDF file and put arrays
        ncout <- nc_create(ofname,list(var_def),force_v4=TRUE)
        # Put variables
        ncvar_put(ncout,var_def,data)
        # Add global attributess
        ncatt_put(ncout,0,"title",title)
        ncatt_put(ncout,0,"institution",institution)
        ncatt_put(ncout,0,"source",source)
        ncatt_put(ncout,0,"CreationDate",as.character(Sys.Date()))
        ncatt_put(ncout,0,"CreatedBy","Mariano S. Alvarez")

        # close the file, writing data to disk
        nc_close(ncout)

      } # end i (number of start dates)
      
      
    } # end iens (ensemble members)
    
    # Close input netcdf File
    nc_close(ncid);
    
  } # end imodel
  
} # end ivar