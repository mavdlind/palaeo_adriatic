####Adriatic Pleistocene####
####load libraries####
library("tidyverse") # data manipulation and plotting
library("rcarbon") # 14C calibration and modelling

####load data####
adriatic<-read.csv("adriatic.csv", header=TRUE) # read data csv file 
adriatic$LON<-as.numeric(as.character(adriatic$LON)) # turning LON values into numeric
adriatic$LAT<-as.numeric(as.character(adriatic$LAT)) # turning LAT values into numeric
adriatic$C14AGE<-as.numeric(as.character(adriatic$C14AGE)) # turning C14AGE values into numeric
adriatic$C14STD<-as.numeric(as.character(adriatic$C14STD)) # turning C14AGE values into numeric
adriatic<-drop_na(adriatic,C14AGE,LAT,LON) # drop rows with no C14AGE, or lat/lon values
adriatic.good<-filter(adriatic,quality=="1"|quality=="2"|quality=="3") # only keep data with appropriate quality
adriatic.good<-filter(adriatic.good,COUNTRY!="Serbia")
adriatic.good<-mutate(adriatic.good,DateID=row_number()) # add id number to each entry
adriatic.tc<-filter(adriatic.good,TC_ID=="OK")
# setup
ncores <- 1
nsim<-250
realstartBCAD <- -52000
realendBCAD <- -9000
bracket <- 100
workingstartBP <- abs(realstartBCAD-1950)+bracket
workingendBP <- abs(realendBCAD-1950)-bracket
if (workingendBP<0){ workingendBP <- 0 }
workingstartCRA <- uncalibrate(workingstartBP)$ccCRA
workingendCRA <- uncalibrate(workingendBP)$ccCRA
if (workingendCRA<0){ workingendCRA <- 0 }
xlim<-c(40000,10000)
#calibrate
adriatic.cal<-calibrate(x=adriatic.good$C14AGE,errors=adriatic.good$C14STD,calCurves = 'intcal20',ids=adriatic.good$DateID)
adriatic.tc.cal<-calibrate(x=adriatic.tc$C14AGE,errors=adriatic.tc$C14STD,calCurves = 'intcal20',ids=adriatic.tc$DateID)
adriatic.bins <- binPrep(sites=adriatic.good$SITENAME, ages=adriatic.good$C14AGE, h=500)
adriatic.spd <- spd(x=adriatic.cal, bins=adriatic.bins, timeRange=c(workingstartBP,workingendBP), datenormalised=FALSE)
plot(adriatic.spd,xlim=xlim)
plot(adriatic.spd,runm=200,add=TRUE,type="simple",col="red",lwd=2,lty=2)
#permutation test based on E-W Adriatic
set.seed(123)
adriatic.permregion <- permTest(x=adriatic.cal, bins=adriatic.bins, marks=adriatic.good$REGION,
                             timeRange=c(workingstartBP,workingendBP), runm=100, nsim=nsim)
par(mfrow=c(2,1))
plot(adriatic.permregion,focalm=1,main="W Adriatic",xlim=xlim)
plot(adriatic.permregion,focalm=2,main="E Adriatic",xlim=xlim)
sink(file="perm_region.txt")
summary(adriatic.permregion)
sink()
#permutation test based on date quality
set.seed(123)
adriatic.permquality <- permTest(x=adriatic.cal, bins=adriatic.bins, marks=adriatic.good$quality,
                                 timeRange=c(workingstartBP,workingendBP), runm=100, nsim=nsim)
par(mfrow=c(3,1))
plot(adriatic.permquality,focalm=1,main="Highly reliable",xlim=xlim)
plot(adriatic.permquality,focalm=2,main="Reliable",xlim=xlim)
plot(adriatic.permquality,focalm=3,main="Slightly reliable",xlim=xlim)
sink(file="perm_quality.txt")
summary(adriatic.permquality)
sink()

####plot final version permutation####
## Plot
xlim <- c(-36000,-12000)
xts <- seq(-36000,2000,1000)
xtl <- seq(-36000,2000,1000)
#dev.new(device=pdf, width=4, height=6.5)
tiff("/home/marc/Dropbox/m_pers/aitor/v2/20210915/fig02.tif",compression="lzw")
layout(matrix(c(1,2,3,4), 4, 1, byrow=TRUE), widths=4, heights=c(1.5,1.5,1.5,1.9))
par(mar=c(0, 2, 0.2, 1)) #c(bottom, left, top, right)
par(yaxs="i")
par(xaxs="i")
plotymax <- max(adriatic.permquality$envelope[["1"]][,2], adriatic.permquality$observed[["1"]]$PrDens)*1.1
tmp <- plot(adriatic.permquality, focalm="1", calendar="BP", col.obs="red", lwd.obs=1, xlim=c(36000,12000), ylim=c(0,plotymax), xaxt='n',yaxt='n',drawaxes=FALSE)
segments(x0=xts, y0=rep(0,length(xts)), x1=xts, y1=rep(0.01*plotymax,length(xts)), lwd=0.5)
segments(x0=xtl, y0=rep(0,length(xtl)), x1=xtl, y1=rep(0.03*plotymax,length(xtl)), lwd=1)
text(x=35900, y=plotymax*0.9, labels="A. Highly reliable", font=2, cex=0.9, adj=c(0,0.7))
plotymax <- par("usr")[4]
legend(x=35900, y=plotymax*0.7,legend=c("SPD","95% MC envelope","positive deviation","negative deviation"),col=c("red","lightgrey",rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)),lty=c(1,1,1,1),lwd=c(1,5,5,5),cex=0.8, bg="white")
axis(side=2, cex.axis=0.8, at=c(0,1))
par(mar=c(0, 2, 0, 1))
plotymax <- max(adriatic.permquality$envelope[["2"]][,2], adriatic.permquality$observed[["2"]]$PrDens)*1.1
plot(adriatic.permquality, focalm="2", calendar="BP", col.obs="orange", lwd.obs=1, xlim=c(36000,12000), ylim=c(0,plotymax*1.1), xaxt='n',yaxt='n',drawaxes=FALSE)
segments(x0=xts, y0=rep(0,length(xts)), x1=xts, y1=rep(0.01*plotymax,length(xts)), lwd=0.5)
segments(x0=xtl, y0=rep(0,length(xtl)), x1=xtl, y1=rep(0.03*plotymax,length(xtl)), lwd=1)
text(x=35900, y=plotymax*0.9, labels="B. Reliable", font=2, cex=0.9, adj=c(0,0.7))
axis(side=2, cex.axis=0.8, at=c(0,1))
par(mar=c(3, 2, 0, 1))
plotymax <- max(adriatic.permquality$envelope[["3"]][,2], adriatic.permquality$observed[["3"]]$PrDens)*1.2
plot(adriatic.permquality, focalm="3", calendar="BP", col.obs="blue", lwd.obs=1, xlim=c(36000,12000), ylim=c(0,plotymax), xaxt='n',yaxt='n',drawaxes=FALSE)
#axis(side=1, at=seq(-36000,12000,1000), labels=seq(-36000,10000,1000), las=2, cex.axis=0.8)
segments(x0=xts, y0=rep(0,length(xts)), x1=xts, y1=rep(0.01*plotymax,length(xts)), lwd=0.5)
segments(x0=xtl, y0=rep(0,length(xtl)), x1=xtl, y1=rep(0.03*plotymax,length(xtl)), lwd=1)
text(x=35900, y=plotymax*0.9, labels="C. Slightly reliable", font=2, cex=0.9, adj=c(0,0.7))
axis(side=2, cex.axis=0.8, at=c(0,1))
xticks <- seq(36000,12000,-2000)
axis(side=1, at=xticks, labels=xticks, las=2, cex.axis=0.8)
par(yaxs="r")
dev.off()

#permutation test based on date type
set.seed(123)
adriatic.permtype <- permTest(x=adriatic.cal, bins=adriatic.bins, marks=adriatic.good$DATETYPE,
                                 timeRange=c(workingstartBP,workingendBP), runm=100, nsim=nsim)
par(mfrow=c(2,1))
plot(adriatic.permtype,focalm=1,main="AMS",xlim=xlim)
plot(adriatic.permtype,focalm=2,main="Non-AMS",xlim=xlim)
sink(file="perm_type.txt")
summary(adriatic.permtype)
sink()

#create exponential model
set.seed(123)
adriatic.expnull <- modelTest(adriatic.cal, errors=adriatic.good$C14STD, bins=adriatic.bins, nsim=1000, runm=50,
                     timeRange=c(workingstartBP,workingendBP), model="exponential", ncores=ncores, datenormalised=FALSE)
sink(file="exponential_null.txt")
summary(adriatic.expnull)
sink()
adriatic.uni<-modelTest(adriatic.cal, errors=adriatic.good$C14STD, bins=adriatic.bins, nsim=1000, runm=50,
                        timeRange=c(workingstartBP,workingendBP), model="uniform", ncores=ncores, datenormalised=FALSE)
sink(file="uniform_null.txt")
summary(adriatic.uni)
sink()
#plot spd vs uniform null and exponential null
layout(matrix(c(1,2,3,4), 4, 1, byrow=TRUE), widths=4.5, heights=c(1.65,1.65,1.65,2))
par(mar=c(0, 1, 1, 1)) #c(bottom, left, top, right)
yMax <- max(adriatic.spd$grid$PrDens)
plot(adriatic.spd, xlim=c(36000,12000), ylim=c(0,yMax),xaxt='n',yaxt='n')
text(x=35900, y=yMax*0.9, labels="A. SPDs", font=2, cex=0.9, adj=c(0,0.7))
legend(x=35900, y=yMax*0.7,legend=c("SPD","95% MC envelope","positive deviation","negative deviation"),col=c("red","lightgrey",rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)),lty=c(1,1,1,1),lwd=c(1,5,5,5),cex=0.8, bg="white")
box()
par(mar=c(0, 1, 0, 1)) #c(bottom, left, top, right)
yMax <- max(adriatic.uni$result$PrDens)*1.1
plot(adriatic.uni, col.obs="darkred", lwd.obs=1, ylim=c(0,yMax), xlim=c(36000,12000), drawaxes=FALSE)
lines(adriatic.uni$fit$calBP,adriatic.uni$fit$PrDens, col="black", lty="dashed", lwd=0.5)
text(x=35900, y=yMax*0.9, labels="B. Uniform model", font=2, cex=0.9, adj=c(0,0.7))
box()
par(mar=c(0, 1, 0, 1)) #c(bottom, left, top, right)
yMax <- max(adriatic.expnull$result$PrDens)*1.1
plot(adriatic.expnull, col.obs="darkred", lwd.obs=1, ylim=c(0,yMax), xlim=c(36000,12000), drawaxes=FALSE)
lines(adriatic.expnull$fit$calBP,adriatic.expnull$fit$PrDens, col="black", lty="dashed", lwd=0.5)
text(x=35900, y=yMax*0.9, labels="C. Exponential growth model", font=2, cex=0.9, adj=c(0,0.7))
xticks <- seq(36000,12000,-2000)
axis(side=1, at=xticks, labels=xticks, las=2, cex.axis=0.8)
box()

####for mapping####
wgs84<-"+proj=longlat +datum=WGS84 +no_defs"


####mapping based on Bayesian modelling (using dates with adequate TC)####
#start gravettian
st.gravettian<-subset(adriatic.cal,BP<=35342 & BP>=33594,p=0.682)
st.gravettian.id<-tibble(DateID=as.integer(st.gravettian$metadata$DateID)) # tibble with id numbers of Caldates
st.gravettian.georef<-inner_join(adriatic.good,st.gravettian.id) # inner join to attach all necessary info to CalDates
points.st.gravettian<-st_as_sf(st.gravettian.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#overlap gravettian - epigravettian
overlap.grav.epig<-subset(adriatic.cal,BP<=26474 & BP>=24232,p=0.682)
overlap.grav.epig.id<-tibble(DateID=as.integer(overlap.grav.epig$metadata$DateID)) # tibble with id numbers of Caldates
overlap.grav.epig.georef<-inner_join(adriatic.good,overlap.grav.epig.id) # inner join to attach all necessary info to CalDates
overlap.grav.epig.georef<-filter(overlap.grav.epig.georef,technocomplex=="Gravettian"|technocomplex=="Early Epigravettian")
points.overlap.grav.epi<-st_as_sf(overlap.grav.epig.georef,coords=c("LON","LAT"),crs=wgs84)
#overlap end early epigravettian - start late epigravettian
overlap.eepig.lepig<-subset(adriatic.cal,BP<=18165 & BP>=17106,p=0.682)
overlap.eepig.lepig.id<-tibble(DateID=as.integer(overlap.eepig.lepig$metadata$DateID)) # tibble with id numbers of Caldates
overlap.eepig.lepig.georef<-inner_join(adriatic.good,overlap.eepig.lepig.id) # inner join to attach all necessary info to CalDates
overlap.eepig.lepig.georef<-filter(overlap.eepig.lepig.georef,technocomplex=="Early Epigravettian"| technocomplex=="Late Epigravettian")
points.overlap.eepig.lepig<-st_as_sf(overlap.eepig.lepig.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
####mapping based on climate events (Rasmussen et al. 2014) (using all 1-3 quality dates regardless of TC)####
#Long = using longest possible error countings
#mean = using mean dates only
#GS3 (long)
gs3.long<-subset(adriatic.cal,BP<=28312 & BP>=22694,p=0.682)
gs3.long.id<-tibble(DateID=as.integer(gs3.long$metadata$DateID)) # tibble with id numbers of Caldates
gs3.long.georef<-inner_join(adriatic.good,gs3.long.id) # inner join to attach all necessary info to CalDates
points.gs3.long<-st_as_sf(gs3.long.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GS3 (mean)
gs3.mean<-subset(adriatic.cal,BP<=27490 & BP>=23290,p=0.682)
gs3.mean.id<-tibble(DateID=as.integer(gs3.mean$metadata$DateID)) # tibble with id numbers of Caldates
gs3.mean.georef<-inner_join(adriatic.good,gs3.mean.id) # inner join to attach all necessary info to CalDates
points.gs3.mean<-st_as_sf(gs3.mean.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GS2.1c (long)
gs2.1c.long<-subset(adriatic.cal,BP<=23483 & BP>=20368,p=0.682)
gs2.1c.long.id<-tibble(DateID=as.integer(gs2.1c.long$metadata$DateID)) # tibble with id numbers of Caldates
gs2.1c.long.georef<-inner_join(adriatic.good,gs2.1c.long.id) # inner join to attach all necessary info to CalDates
points.gs2.1c.long<-st_as_sf(gs2.1c.long.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GS2.1c (mean)
gs2.1c.mean<-subset(adriatic.cal,BP<=23290 & BP>=20850,p=0.682)
gs2.1c.mean.id<-tibble(DateID=as.integer(gs2.1c.mean$metadata$DateID)) # tibble with id numbers of Caldates
gs2.1c.mean.georef<-inner_join(adriatic.good,gs2.1c.mean.id) # inner join to attach all necessary info to CalDates
points.gs2.1c.mean<-st_as_sf(gs2.1c.mean.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GS2.1b (long)
gs2.1b.long<-subset(adriatic.cal,BP<=21332 & BP>=17100,p=0.682) # tibble with id numbers of Caldates
gs2.1b.long.id<-tibble(DateID=as.integer(gs2.1b.long$metadata$DateID)) # tibble with id numbers of Caldates
gs2.1b.long.georef<-inner_join(adriatic.good,gs2.1b.long.id) # inner join to attach all necessary info to CalDates
points.gs2.1b.long<-st_as_sf(gs2.1b.long.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GS2.1b (mean)
gs2.1b.mean<-subset(adriatic.cal,BP<=20850 & BP>=17430,p=0.682)
gs2.1b.mean.id<-tibble(DateID=as.integer(gs2.1b.mean$metadata$DateID)) # tibble with id numbers of Caldates
gs2.1b.mean.georef<-inner_join(adriatic.good,gs2.1b.mean.id) # inner join to attach all necessary info to CalDates
points.gs2.1b.mean<-st_as_sf(gs2.1b.mean.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GS2.1a (long)
gs2.1a.long<-subset(adriatic.cal,BP<=17760 & BP>=14456,p=0.682)
gs2.1a.long.id<-tibble(DateID=as.integer(gs2.1a.long$metadata$DateID)) # tibble with id numbers of Caldates
gs2.1a.long.georef<-inner_join(adriatic.good,gs2.1a.long.id) # inner join to attach all necessary info to CalDates
points.gs2.1a.long<-st_as_sf(gs2.1a.long.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GS2.1a (mean)
gs2.1a.mean<-subset(adriatic.cal,BP<=17430 & BP>=14642,p=0.682)
gs2.1a.mean.id<-tibble(DateID=as.integer(gs2.1a.mean$metadata$DateID)) # tibble with id numbers of Caldates
gs2.1a.mean.georef<-inner_join(adriatic.good,gs2.1a.mean.id) # inner join to attach all necessary info to CalDates
points.gs2.1a.mean<-st_as_sf(gs2.1a.mean.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GI1 (long)
gi1.long<-subset(adriatic.cal,BP<=14828 & BP>=12708,p=0.682) # tibble with id numbers of Caldates
gi1.long.id<-tibble(DateID=as.integer(gi1.long$metadata$DateID)) # tibble with id numbers of Caldates
gi1.long.georef<-inner_join(adriatic.good,gi1.long.id) # inner join to attach all necessary info to CalDates
points.gi1.long<-st_as_sf(gi1.long.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GI1 (mean)
gi1.mean<-subset(adriatic.cal,BP<=14642 & BP>=12846,p=0.682)
gi1.mean.id<-tibble(DateID=as.integer(gi1.mean$metadata$DateID)) # tibble with id numbers of Caldates
gi1.mean.georef<-inner_join(adriatic.good,gi1.mean.id) # inner join to attach all necessary info to CalDates
points.gi1.mean<-st_as_sf(gi1.mean.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GS1 (long)
gs1.long<-subset(adriatic.cal,BP<=12708,BP>=11554,p=0.682)
gs1.long.id<-tibble(DateID=as.integer(gs1.long$metadata$DateID)) # tibble with id numbers of Caldates
gs1.long.georef<-inner_join(adriatic.good,gs1.long.id) # inner join to attach all necessary info to CalDates
points.gs1.long<-st_as_sf(gs1.long.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping
#GS1 (mean)
gs1.mean<-subset(adriatic.cal,BP<=12846,BP>=11703,p=0.682)
gs1.mean.id<-tibble(DateID=as.integer(gs1.mean$metadata$DateID)) # tibble with id numbers of Caldates
gs1.mean.georef<-inner_join(adriatic.good,gs1.mean.id) # inner join to attach all necessary info to CalDates
points.gs1.mean<-st_as_sf(gs1.mean.georef,coords=c("LON","LAT"),crs=wgs84) # sf object for mapping

