########################################################################################################
#CLOSED-POPULATION MARK-RECAPTURE MODEL: ESTIMATING NUMBER OF CARCASSES FROM VULTURE SATELLITE DATA    #
########################################################################################################

library(rgdal)
library(raster)
library(R2WinBUGS)
library(logger)

#################################################################################

#### No data provided for first part of script as it contains sensitive information. Code is still shown for understanding, and sample data is provided later in the script for use


#1. Load data and prep for creating dataframes:
## Load Carcass data from LDA (with XY coordinates, spatial variables and region/protected area) now ##

# # Which vultures were in which regions on a given day? 
GPSdata$Date<-as.Date(GPSdata$timestamp)
GPSdata_spdf<-GPSdata
coordinates(GPSdata_spdf)<-~X+Y
proj4string(GPSdata_spdf)<-CRS("+init=epsg:4326")
GPSdata_spdf_t<-spTransform(GPSdata_spdf,crs(regions))
regions$combo_name<-regions$Region_Nam
regions$combo_name[is.na(regions$ORIG_NAME)==F]<-regions$ORIG_NAME[is.na(regions$ORIG_NAME)==F]

#Use national park name (rather than the region name) where it exists
GPSdata_spdf_t$Region_2<-GPSdata$Region_2<-over(GPSdata_spdf_t,regions)$combo_name #Much quicker than looping 'over'!

#Omit regions where too few data and the convergence is poor
cc_regions_to_omit<-c("Kigoma","Kizigo G.R. (C)","Mikumi National Park","Rukwa","Tabora")

Carcass.cleaned<-Carcass.cleaned[(Carcass.cleaned$Region %in% cc_regions_to_omit)==F,]
GPSdata_cc<-GPSdata[(GPSdata$Region %in% cc_regions_to_omit)==F,]

#Remove birds from GPSdata which are not in Carcass.cleaned (did not find identified carcasses while tagged)
unccb<-unique(Carcass.cleaned$MovebankID)
ungpsb_cc<-unique(GPSdata_cc$MovebankID)
birds_missing_from_GPSdata_cc<-ungpsb_cc[which((ungpsb_cc %in% unccb)==F)]
GPSdata_cc<-GPSdata_cc[which((GPSdata_cc$MovebankID %in% birds_missing_from_GPSdata_cc)==F),]
unccb[which((unccb %in% ungpsb_cc)==F)] #Check no birds in Carcass.cleaned that are not in GPSdata

#Make bird and carcass names unambiguously *character* strings so as to make querying capture history easier
GPSdata_cc$bird_name<-paste0("b_",GPSdata_cc$MovebankID)

Carcass.cleaned$bird_name<-paste0("b_",Carcass.cleaned$MovebankID)
Carcass.cleaned$carcass_name<-paste0("c_",Carcass.cleaned$clusterId)
un_birds_cc<-unique(as.character(Carcass.cleaned$bird_name))
un_carcs_cc<-unique(as.character(Carcass.cleaned$carcass_name))

################################################
##2. Capture history and availability history.
################################################

#Full carcass dataset first
#First, get the region of each carcass
Carcass.cleaned_spdf<-Carcass.cleaned
coordinates(Carcass.cleaned_spdf)<-~X+Y
proj4string(Carcass.cleaned_spdf)<-CRS("+init=epsg:4326")
Carcass.cleaned_spdf_t<-spTransform(Carcass.cleaned_spdf,crs(regions))
Carcass.cleaned$Region_2<-character(nrow(Carcass.cleaned))

#Fill in regions
Carcass.cleaned$Region_2<-over(Carcass.cleaned_spdf_t,regions)$combo_name

#Make full capture history. Which carcasses were detected by which birds?
#Rows are carcasses, columns are birds
cap_hist<-data.frame(matrix(0,ncol=length(un_birds_cc),nrow=length(un_carcs_cc)))
names(cap_hist)<-un_birds_cc
row.names(cap_hist)<-un_carcs_cc
for(i in 1:nrow(Carcass.cleaned)){
  cap_hist[Carcass.cleaned[i,]$carcass_name,Carcass.cleaned[i,]$bird_name]<-1
}

#Availability history: which carcasses could have been detected by which birds?
#If a vulture was present in the same region as the carcass during the period for which the carcass was known, the carcass could have been detected by that vulture.

#Create table of dates when each carcass was known.
carc_date_tab<-data.frame(un_carcs_cc)
carc_date_tab[,2:4]<-character(nrow(carc_date_tab))
names(carc_date_tab)<-c("carcass","min_date","max_date","region")
for(i in 1:nrow(carc_date_tab)){
  temp_first_dates<-Carcass.cleaned[Carcass.cleaned$carcass_name==carc_date_tab[i,]$carcass,]$Date
  carc_date_tab$min_date[i]<-as.character(min(temp_first_dates))
  temp_last_dates<-Carcass.cleaned[Carcass.cleaned$carcass_name==carc_date_tab[i,]$carcass,]$LAST_date
  carc_date_tab$max_date[i]<-as.character(max(temp_last_dates))
  carc_date_tab$region[i]<-unique(Carcass.cleaned[Carcass.cleaned$carcass_name==carc_date_tab[i,]$carcass,]$Region_2)
}
carc_date_tab$min_poss_date<-as.Date(carc_date_tab$min_date) 
carc_date_tab$max_poss_date<-as.Date(carc_date_tab$max_date)

#Populate availability history
avail_hist<-data.frame(matrix(0,dim(cap_hist)[1],dim(cap_hist)[2]))
names(avail_hist)<-names(cap_hist)
row.names(avail_hist)<-row.names(cap_hist)

#For each carcass, which birds were in that region within the time it was detectable?
GPSdata_cc$region_date<-paste0(GPSdata_cc$Region_2,GPSdata_cc$Date)

for(i in 1:nrow(carc_date_tab)){
  temp_region_dates<-paste0(carc_date_tab[i,]$region,as.Date(carc_date_tab[i,]$min_poss_date:carc_date_tab[i,]$max_poss_date,origin="1970-01-01"))
  temp_avail_birds<-unique(GPSdata_cc[which(GPSdata_cc$region_date %in% temp_region_dates),]$bird_name)
  avail_hist[carc_date_tab$carcass[i],temp_avail_birds]<-1
}

table(apply(avail_hist,1,sum)) #Check no zeros
table(cap_hist[which(avail_hist==0,arr.ind=T)]) #Check no 1s

#Full carcass dataset regions
un_regions_cc<-unique(Carcass.cleaned$Region_2)

###############################################
####DATA PROVIDED FROM THIS SECTION ONWARDS
###############################################
load("CRdata.Rdata")

###############################################
## 3.  Analyse all regions at the same time. 
###############################################

#Need the number of rows to be the same for all regions' capture histories. So set this as nz + the maximum number of known carcasses in any region
nz<- 1000 
max_carcs<-max(table(carc_date_tab$region))
cap_hist_allregions<-avail_hist_allregions<-array(dim=c(length(un_regions_cc),nz+max_carcs,dim(cap_hist)[2])) #x regions, y carcasses, z vultures
n_all_carcasses<-l_ever_in_region<-numeric(length(un_regions_cc)) #M and T in the model. n_all_carcasses = nz + the maximum number of known carcasses in a region.
#We will only fill in the first n_all_carcasses row for each capture/availability history

n_carcasses<-numeric(length(un_regions_cc)) #For results table
vult_ind<-matrix(nrow=length(un_regions_cc),ncol=length(un_birds_cc))

for(i in 1:length(un_regions_cc)){
  temp_carcs<-unique(Carcass.cleaned[Carcass.cleaned$Region_2==un_regions_cc[i],]$carcass_name)
  temp_cap_hist<-cap_hist[temp_carcs,]
  temp_avail_hist<-avail_hist[temp_carcs,]
  
  ever_in_region<-as.numeric(which(apply(temp_avail_hist,2,sum)>0))
  never_in_region<-as.numeric(which(apply(temp_avail_hist,2,sum)==0))
  l_ever_in_region[i]<-length(ever_in_region)
  
  temp_cap_hist_resorted<-cbind(temp_cap_hist[,ever_in_region],temp_cap_hist[,never_in_region])
  temp_avail_hist_resorted<-cbind(temp_avail_hist[,ever_in_region],temp_avail_hist[,never_in_region])
  
  #Prep for MR analysis: augment dataset for (known and potential) carcasses. 
  #WinBUGS needs the same number of rows for capture histories across regions
  temp_y<-as.matrix(temp_cap_hist_resorted)
  n_carcasses[i]<-nrow(temp_y)
  temp_ah<-as.matrix(temp_avail_hist_resorted)
  names(temp_y)<-names(temp_ah)<-NULL
  row.names(temp_y)<-row.names(temp_ah)<-NULL
  temp_yaug<-rbind(temp_y,array(0,dim=c(nz,ncol(temp_y))))
  temp_ahaug<-rbind(temp_ah,array(NA,dim=c(nz,ncol(temp_ah)))) #Has to be NAs, because we don't *know* whether the other birds were in the area at the time that the imaginary carcasses were there
  n_all_carcasses[i]<-nrow(temp_yaug)
  
  cap_hist_allregions[i,1:n_all_carcasses[i],]<-temp_yaug 
  avail_hist_allregions[i,1:n_all_carcasses[i],]<-temp_ahaug 
  
  #Work out the index of each vulture in the reshuffled capture history
  vult_ind[i,]<-c(ever_in_region,never_in_region)
}

#Calculate with Separate Capture Probabilities by Region
sink("model.bug")
cat("
    model{
    #Priors
    # mean.p~dnorm(0,0.0001) #For global intercept
    for(h in 1:nr){
    omega[h]~dunif(0,1)
    # p.off[h]~dnorm(0,0.0001) #For global intercept
    p[h]~dunif(0,1)
    for(i in 1:(M[h])){
    for(j in 1:(T[h])){
    ahaug[h,i,j]~dunif(0,1)
    }
    }
    }
    
    #Likelihood
    for(h in 1:nr){
    # logit(p[h])<-mean.p+p.off[h]
    for(i in 1:(M[h])){ 
    z[h,i]~dbern(omega[h]) #Inclusion indicators
    for(j in 1:(T[h])){ 
    yaug[h,i,j]~dbern(p.eff[h,i,j])
    p.eff[h,i,j]<-z[h,i]*p[h]*ahaug[h,i,j]
    } #j
    } #i
    } #h
    
    #Derive population size estimates for each region
    for(h in 1:nr){
    N[h]<-sum(z[h,1:(M[h])])
    }
    }
    ",fill=T)
sink()

#Bundle data
win.data<-list(yaug=cap_hist_allregions,ahaug=avail_hist_allregions,M=n_all_carcasses,T=l_ever_in_region,nr=length(un_regions_cc))

#Initial values
inits<-function() list(p=runif(length(un_regions_cc),0,1))

#Parameters monitored
params<-c("p","omega","N")

#MCMC settings
ni<-100000 #Runs at about 25k iterations/hr
nt<-20
nb<-50000
nc<-3


### INSERT YOUR DIRECTORY HERE
bugs.dir<- ##

# This takes considerable time to run: Output data provided below
out<-bugs(win.data,inits,params,"model.bug",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb,debug=F,bugs.directory=bugs.dir,working.directory=getwd())

###################################
## Examine Output
###################################

load("Sample_CR_Output.RData")
load("CRdata2.Rdata")

#View output - check Rhat is less than 1.1
out$summary

#Examine chains for all parameters
params<-row.names(data.frame(out$summary))
N_params_ind<-which(substr(params,1,1)=="N")

# Add file directory to output figures into
for(i in 1:length(un_regions_cc)){
  png(paste0("YOUR_FILE_DIRECTORY/Chains_N_",un_regions_cc[i],".png"),8,6,res=500,units="in")
  plot(out$sims.array[,1,N_params_ind[i]],type="l",ylim=c(min(out$sims.array[,,N_params_ind[i]])-5,max(out$sims.array[,,N_params_ind[i]])+5),col="red",xlab="Iterations (thinned) after burn-in",ylab="True number of carcasses",main=paste0("Convergence: ",un_regions_cc[i]))
  lines(out$sims.array[,2,N_params_ind[i]],col="blue")
  lines(out$sims.array[,3,N_params_ind[i]],col="green")
  dev.off()
}

# Create summary table for results 
vulture_integ_MR_results<-data.frame(matrix(nrow=length(un_regions_cc),ncol=8))
names(vulture_integ_MR_results)<-c("region","n_vultures","obsd_carcs","estd_carcs_q.025","estd_carcs_mean","estd_carcs_q.975","p","omega")
vulture_integ_MR_results$region<-un_regions_cc
vulture_integ_MR_results$n_vultures<-l_ever_in_region
vulture_integ_MR_results$obsd_carcs<-n_carcasses
vulture_integ_MR_results$estd_carcs_mean<-out$mean$N
vulture_integ_MR_results$p<-out$mean$p
vulture_integ_MR_results$omega<-out$mean$omega
for(i in 1:nrow(vulture_integ_MR_results)){
  vulture_integ_MR_results$estd_carcs_q.025[i]<-quantile(out$sims.list$N[,i],0.025)
  vulture_integ_MR_results$estd_carcs_q.975[i]<-quantile(out$sims.list$N[,i],0.975)
}


#Add file directory to output posterior plots for all parameters
for(i in 1:length(un_regions_cc)){
  png(paste0("YOUR_FILE_DIRECTORY/N_",un_regions_cc[i],".png"),8,6,res=500,units="in")
  hist(out$sims.list$N[,i],xlab="Posterior true number of carcasses",xlim=c(n_carcasses[i]-2,max(out$sims.list$N[,i])),main=paste0(un_regions_cc[i]," (dashed = observed)"))
  abline(v=n_carcasses[i],lwd=2,lty=2)
  dev.off()
}


