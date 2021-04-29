#SDM for Dracocephalum ruyshiana

library(raster)
library(rasterVis)
library(rgdal)
library(sp)
library(sdm)
library(gridExtra)
library(readxl)#Read in xls
library(RColorBrewer)

#Norway outline
norway<-getData('GADM',level=0,country='NOR')

#Environmental variables
list.files('Data/',pattern='tif')

#Stack at 1km resolution
predvars<-stack(list.files('Data/',pattern='tif',full.names=T))
predvars
#Tidy layer names
names(predvars)<-sub("X_","",names(predvars))
plot(predvars)

#Convert temp to degrees
predvars$bio10_16<-predvars$bio10_16/10
#Convert soil pH to decimal ph
predvars$SoilpH<-predvars$SoilpH/10


norwayP<-spTransform(norway,CRS=crs(predvars))


#Importing Species data
darcep<- read_excel("Data/Dracocephalum ruyschiana_GBIF_data_edited.xlsx",skip=4)
View(darcep)

#Removing points without coordinates
darcep<-darcep[!is.na(darcep$decimalLatitude) & !is.na(darcep$decimalLongitude),]

#Removing points marked as outliers
darcep<-darcep[is.na(darcep$Outliers),]

#Removing fossil records
summary(as.factor(darcep$basisOfRecord))
darcep<-darcep[darcep$basisOfRecord!='FOSSIL_SPECIMEN',]

dim(darcep)

#Convert to spatial points dataframe
darcep_sp<-SpatialPointsDataFrame(cbind(as.numeric(darcep$decimalLongitude),as.numeric(darcep$decimalLatitude)),darcep,proj4string=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
#Project to utm grid
darcep_utm<-spTransform(darcep_sp,CRS=crs(predvars))

plot(norway)
points(darcep_sp,cex=0.1,pch=16,col=2)

#Simple df
#darcep_simple<-data.frame(darcep_utm[names(darcep_utm)%in%c('decimalLongitude','decimalLatitude','species')])
darcep_simple<-data.frame(cbind(darcep_utm@coords,species=darcep$species))

#Plots of env vars with occurrence records
pt<-levelplot(predvars$bio10_16,scales=list(draw=F),margin=F,par.settings='YlOrRdTheme',main=expression('Mean summer temperature'~(degree~C)),
              colorkey=list(space='bottom',title=''))+
  latticeExtra::layer(sp.polygons(norwayP,col=grey(0.5)))+ 
  latticeExtra::layer(sp.points(darcep_utm,col=1,pch=16,cex=0.5))

pps<-levelplot(predvars$bio15_16,scales=list(draw=F),margin=F,par.settings='RdBuTheme',main='Precipitation seasonality',
               colorkey=list(space='bottom',title=expression('')))+
  latticeExtra::layer(sp.polygons(norwayP,col=grey(0.5)))+ 
  latticeExtra::layer(sp.points(darcep_utm,col=1,pch=16,cex=0.5))

#pph<-levelplot(mask(predvars$SoilpH,predvars$bio10_16),scales=list(draw=F),margin=F,par.settings='RdBuTheme',main='Soil pH',
#               colorkey=list(space='bottom',title=expression('')))+
#  layer(sp.polygons(norwayP,col=grey(0.5)))+ 
#  layer(sp.points(darcep_utm,col=1,pch=16,cex=0.5))

pp<-levelplot(predvars$bio12_16,scales=list(draw=F),margin=F,par.settings='RdBuTheme',main=expression('Mean annual precipitation'~(mm)),
              colorkey=list(space='bottom'))+
  latticeExtra::layer(sp.polygons(norwayP,col=grey(0.5)))+ 
  latticeExtra::layer(sp.points(darcep_utm,col=1,pch=16,cex=0.5))


#Background data
background_dat<-read.table('Data/BackgroundBiasCorrected.txt',header=T)

background_utm<-SpatialPoints(background_dat,proj4string=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

pB<-levelplot(predvars[[1]]*0,col.regions=grey(1),colorkey=F,scales=list(draw=F),margin=F,par.settings='YlOrRdTheme',main='Background data')+
  latticeExtra::layer(sp.polygons(norwayP,col=grey(0.5)))+
  latticeExtra::layer(sp.points(background_utm,pch=1,cex=0.1,col=1))
pB

tiff('Figures/AllData.tif',width = 1600,height = 1000,pointsize = 20,units = 'px')
#pdf('Figures/AllData.pdf',width = 16,height = 10,pointsize = 20)
grid.arrange(pt,pp,pps,ncol=3)
#dev.off()
dev.off()

#SDM


#Make the sdm dataset

#Single object
darcep_simp <- darcep_utm[,names(darcep_utm)%in%c("decimalLatitude","decimalLongitude")]
darcep_simp$species<-rep('Dracocephalum_ruyschiana',times=nrow(darcep_utm))

sdmdataset<-sdmData(species~
                      +bio10_16+bio12_16+bio15_16,
                    train=darcep_simp,
                    predictors=predvars,
                    bg=list(sample(background_utm,1000),remove=T)
                    #bg=sample(background_utm,1000)
                                        )
sdmdataset

#SDM
sdm1<-sdm(Dracocephalum_ruyschiana~
                 bio10_16+bio12_16+bio15_16,
               data=sdmdataset,
               methods=c('glm','gam','rf','gbm','mda','fda','brt'),
               replication=c('cv'),cv.folds=5)
sdm1

modeval1<-cbind(sdm1@run.info,getEvaluation(sdm1))
with(modeval1,tapply(AUC,species,mean))
with(modeval1,tapply(AUC,species,sd))

#Variables importances
varimplist<-list()
for (i in 1:max(sdm1@run.info$modelID)){
  ifelse(sdm1@run.info$success[i]==TRUE,
         {
           ifelse(!is.null(getVarImp(sdm1,id=i)),
                  varimplist[[i]]<-getVarImp(sdm1,id=i)@varImportance,
                  varimplist[[i]]<-df1)
           varimplist[[i]]$species<-sdm1@run.info$species[i]
           varimplist[[i]]$method<-sdm1@run.info$method[i]
           varimplist[[i]]$repid<-sdm1@run.info$replicationID[i]}
         ,print(paste('Model failiure run ',i)))
}
VarImpdf<-do.call('rbind',varimplist)

#Summarise by species
sem<-function(x)sd(x,na.rm=T)/sqrt(length(!is.na(x)))
varimp_mean<-with(VarImpdf,tapply(corTest,list(species,variables),mean))  
varimp_sem<-with(VarImpdf,tapply(corTest,list(species,variables),sem))  

varimp_mean
varimp_sem


tiff('Figures/VarImp.tif',width=800,height=600,pointsize = 20)
#pdf('Figures/VarImp.pdf',width=9,height=7,pointsize = 20)
par(mar=c(5,5,1,1))
b1<-barplot(varimp_mean,beside=T,ylab='Variable importance',names.arg=c('MST','MAP','Precipitation \n seasonality'),las=1,ylim=c(0,0.80),
            legend.text=sub("_"," ",rownames(varimp_mean)))     
arrows(b1,varimp_mean+varimp_sem,b1,varimp_mean-varimp_sem,length=0.05,code=3,angle=90)
dev.off()
dev.off()

#Response curves
responsecurvelist<-list()
for (i in 1:length(levels(as.factor(sdm1@run.info$species)))){
  responsecurvelist[[i]]<-getResponseCurve(sdm1,id=sdm1@run.info$modelID[sdm1@run.info$species==levels(as.factor(sdm1@run.info$species))[i]
                                                                                   &sdm1@run.info$method%in% c('glm','gam','brt')]
                                           ,mean=T,main=levels(as.factor(sdm1@run.info$species))[i])
}

tiff('Figures/ResponseCurves.tif',width = 1200,height=800,units='px',pointsize = 20)
#pdf('Figures/ResponseCurves.pdf',width = 13,height=13,pointsize = 20)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
plot(responsecurvelist[[1]]@response$bio10_16[,1],apply(responsecurvelist[[1]]@response$bio10_16[,2:ncol(responsecurvelist[[1]]@response$bio10_16)],1,mean)
     ,type='l',xlab=expression('MST'~(degree~C)),ylab='Response',las=1,ylim=c(0,1))
lines(responsecurvelist[[1]]@response$bio10_16[,1],apply(responsecurvelist[[1]]@response$bio10_16[,2:ncol(responsecurvelist[[1]]@response$bio10_16)],1,mean)
      +apply(responsecurvelist[[1]]@response$bio10_16[,2:ncol(responsecurvelist[[1]]@response$bio10_16)],1,sd),lty=2)
lines(responsecurvelist[[1]]@response$bio10_16[,1],apply(responsecurvelist[[1]]@response$bio10_16[,2:ncol(responsecurvelist[[1]]@response$bio10_16)],1,mean)
      -apply(responsecurvelist[[1]]@response$bio10_16[,2:ncol(responsecurvelist[[1]]@response$bio10_16)],1,sd),lty=2)
#legend('topl',lty=1,col=c(1,2),c('Carex lepidocarpa','Carex jemtlandica'))

plot(responsecurvelist[[1]]@response$bio12_16[,1],apply(responsecurvelist[[1]]@response$bio12_16[,2:ncol(responsecurvelist[[1]]@response$bio12_16)],1,mean)
     ,type='l',xlab='MAP (mm)',ylab='Response',las=1,ylim=c(0,1))
lines(responsecurvelist[[1]]@response$bio12_16[,1],apply(responsecurvelist[[1]]@response$bio12_16[,2:ncol(responsecurvelist[[1]]@response$bio12_16)],1,mean)
      +apply(responsecurvelist[[1]]@response$bio12_16[,2:ncol(responsecurvelist[[1]]@response$bio12_16)],1,sd),lty=2)
lines(responsecurvelist[[1]]@response$bio12_16[,1],apply(responsecurvelist[[1]]@response$bio12_16[,2:ncol(responsecurvelist[[1]]@response$bio12_16)],1,mean)
      -apply(responsecurvelist[[1]]@response$bio12_16[,2:ncol(responsecurvelist[[1]]@response$bio12_16)],1,sd),lty=2)


plot(responsecurvelist[[1]]@response$bio15_16[,1],apply(responsecurvelist[[1]]@response$bio15_16[,2:ncol(responsecurvelist[[1]]@response$bio15_16)],1,mean)
     ,type='l',xlab='Precipitation seasonality',ylab='Response',las=1,ylim=c(0,1))
lines(responsecurvelist[[1]]@response$bio15_16[,1],apply(responsecurvelist[[1]]@response$bio15_16[,2:ncol(responsecurvelist[[1]]@response$bio15_16)],1,mean)
      +apply(responsecurvelist[[1]]@response$bio15_16[,2:ncol(responsecurvelist[[1]]@response$bio15_16)],1,sd),lty=2)
lines(responsecurvelist[[1]]@response$bio15_16[,1],apply(responsecurvelist[[1]]@response$bio15_16[,2:ncol(responsecurvelist[[1]]@response$bio15_16)],1,mean)
      -apply(responsecurvelist[[1]]@response$bio15_16[,2:ncol(responsecurvelist[[1]]@response$bio15_16)],1,sd),lty=2)

dev.off()
dev.off()

#Predictions
drapred<-predict(sdm1,predvars,filename='Outputs/Model Predictions/dracep',
                 species='Dracocephalum_ruyshiana',mean=T,overwrite=T)
drapred  
drapred_mean<-calc(drapred,mean)
levelplot(drapred_mean,scales=list(draw=F),margin=F)

#Ensembling
ens_dracep<-ensemble(sdm1,newdata = predvars,filename = 'Outputs/Model Predictions/dracep_ensemble',
                  setting=list(method='weighted',stat='AUC'),overwrite=T)
                           
ens_dracep<-raster('Outputs/Model Predictions/dracep_ensemble')

tiff('Figures/HSMMap.tif',width=1200,height=800,units='px',pointsize = 20,res=300)
#pdf('Figures/HSMMap.pdf',width=1200,height=800,pointsize = 20)
levelplot(ens_dracep,scales=list(draw=F),margin=F,par.settings='YlOrRdTheme')+
  latticeExtra::layer(sp.polygons(norwayP,col=grey(0.5),lwd=0.02))+
  latticeExtra::layer(sp.points(darcep_utm,pch=1,cex=0.1,col=1))
dev.off()
dev.off()

kmlcols<- brewer.pal(9, 'YlOrRd')
KML(projectRaster(ens_dracep,res=0.008333333,crs=norway@proj4string),
    'Outputs/Dracocephalum_ruyschiana_climatesuitabilitymodel',overwrite=T,
    col=kmlcols)

#Niche (MST and precip season)
niche_temp_pc_dracep<-sdm::niche(predvars,ens_dracep,c('bio10_16','bio15_16'),plot=F)

#Stack up niches

levelplot(niche_temp_pc_dracep,par.settings='YlOrRdTheme')
#Set extent as actual climate variables
extent(niche_temp_pc_dracep)<-stack(niche_temp_pc_dracep@scaleParams)[,1]
tiff('Figures/Niches.tif',width=800,height = 800,res=150)
pdf('Figures/Niches.pdf',width=8,height = 8)
levelplot(niche_temp_pc_dracep,par.settings='YlOrRdTheme',
          xlab=expression('MST'~(degree~C)),ylab='Precipitation seasonality',
          cex=0.8)
dev.off()
dev.off()

niche_temp_map_dracep<-sdm::niche(predvars,ens_dracep,c('bio10_16','bio12_16'),plot=F)
plot(niche_temp_map_dracep)
levelplot(raster(niche_temp_map_dracep),par.settings='YlOrRdTheme')
#Set extent as actual climate variables
extent(niche_temp_map_dracep)<-stack(niche_temp_map_dracep@scaleParams)[,1]
levelplot(niche_temp_map_dracep,par.settings='YlOrRdTheme',
          xlab=expression('MST'~(degree~C)),ylab='Precipitation (mm)',
          cex=0.8)

lp1<-levelplot(niche_temp_map_dracep,par.settings='YlOrRdTheme',
               xlab=expression('MST'~(degree~C)),ylab='Precipitation (mm)',
               names.attr=c('C. lepidocarpa','C. jemtlandica'),cex=0.8)
lp2<-levelplot(niche_temp_ps_dracep,par.settings='YlOrRdTheme',
               xlab=expression('MST'~(degree~C)),ylab='Precipitation seasonality',
               names.attr=c('C. lepidocarpa','C. jemtlandica'),cex=0.8)
#tiff('Figures/Niches2.tif',height = 1200,width = 1000,units='px',res=100)
grid.arrange(lp1,lp2)
#dev.off()