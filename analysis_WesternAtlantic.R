######################################################################
##’ @title Area-normalized scaling of ReefBudget calcification and bioerosion rates for use with CoralNet
##’
##’ @author Travis A Courtney, PhD
##’ @contact travis.courtney@upr.edu
##’ @date 2024-07-30
##’ @log Version 2.0 - Please see README for a description of updates from Version 1.0

###########include or exclude bioerosion##########
#bioerosion rates are included by default, set include_bioerosion to FALSE to calculate gross carbonate production only
include_bioerosion=TRUE

###########require packages and scripts##########
library(fs)
library(tidyverse)
library(readxl)
library(triangle)

#turn off scientific notation
options(scipen=999)

#set seed so the code generates consistent rates between runs
set.seed(999)

#########CoralNet Western Atlantic Calcification Rates##########
#Import Required Data
setwd("Data_WesternAtlantic")
ReefBudget_Atlantic <- read_excel("Caribbean carbonate production template V2.3-June24.xlsx", sheet = "Calcification Rates", skip = 9) #import Indo-Atlantic ReefBudget rates
Gonzalez_Barrios <- read_csv("Gonzalez-Barios_Alvarez-Filip_2019_Rugosity.csv") #Import Gonzalez-Barrios taxa-level estimates of colony structural complexity
CoralNet_Atlantic <- read_csv("CoralNet_labelset.csv") #Import CoralNet labels
microbioerosion=read_excel("Caribbean carbonate production template V2.3-June24.xlsx", sheet = "Microbioerosion") #import Caribbean ReefBudget microbioerosion data
macrobioerosion=read_excel("Caribbean carbonate production template V2.3-June24.xlsx", sheet = "Macrobioerosion", skip = 3) #import Caribbean 
NCRMP_Florida <- read_csv("CRCP_Coral_Demographic_Survey_Florida_4525_37fd_301a.csv") %>% 
  select(SPECIES_NAME,MAX_DIAMETER)
NCRMP_PuertoRico <- read_csv("CRCP_Coral_Demographic_Survey_Puerto_Rico_f50a_5235_ae5f.csv") %>% 
  select(SPECIES_NAME,MAX_DIAMETER)
NCRMP_USVI <- read_csv("CRCP_Coral_Demographic_Survey_USVI_c475_5023_3b2c.csv") %>% 
  select(SPECIES_NAME,MAX_DIAMETER)
setwd("..")

#Reorganize ReefBudget dataframe
colnames(ReefBudget_Atlantic)=c("code","name","morphology","ext_mean","ext_sd","dens_mean","dens_sd","cf","coeff_mean","coeff_lwr","coeff_upr","int_mean","int_lwr","int_upr","notes")
ReefBudget_Atlantic <- ReefBudget_Atlantic[!is.na(ReefBudget_Atlantic$coeff_mean),]

#Merge rates with taxa-specific rugosity index from Gonzalez-Barrios data
rugosity_rates= merge(ReefBudget_Atlantic,Gonzalez_Barrios,by="name",all.x=TRUE)

#use Madracis formosa rugosity for Madracis carmabi
rugosity_rates[!is.na(str_match(rugosity_rates$name,"Madracis carmabi")),]$rugosity_mean=rugosity_rates[!is.na(str_match(rugosity_rates$name,"Madracis formosa")),]$rugosity_mean
rugosity_rates[!is.na(str_match(rugosity_rates$name,"Madracis carmabi")),]$rugosity_sd=rugosity_rates[!is.na(str_match(rugosity_rates$name,"Madracis formosa")),]$rugosity_sd

#find taxa with complete and missing rugosity data
complete_rugosity=rugosity_rates[!is.na(rugosity_rates$rugosity_mean),]
missing_rugosity=rugosity_rates[is.na(rugosity_rates$rugosity_mean),]

#replace missing SD data with zeros
complete_rugosity[!is.na(str_match(complete_rugosity$rugosity_sd,"nd")),]$rugosity_sd=0

#replace missing branching morphologies with mean±SD from Gonzalez-Barrios branching rugosities
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"bra")),]$rugosity_mean=mean(subset(Gonzalez_Barrios,morphology=="Branching")$rugosity_mean)
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"bra")),]$rugosity_sd=sd(subset(Gonzalez_Barrios,morphology=="Branching")$rugosity_mean)
#replace missing plating morphologies with mean±SD from Gonzalez-Barrios foliose rugosities
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"plat")),]$rugosity_mean=mean(subset(Gonzalez_Barrios,morphology=="Foliose")$rugosity_mean)
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"plat")),]$rugosity_sd=sd(subset(Gonzalez_Barrios,morphology=="Foliose")$rugosity_mean)
#replace missing massive and submassive morphologies with mean±SD from Gonzalez-Barrios massive rugosities
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"mass")),]$rugosity_mean=mean(subset(Gonzalez_Barrios,morphology=="Massive")$rugosity_mean)
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"mass")),]$rugosity_sd=sd(subset(Gonzalez_Barrios,morphology=="Massive")$rugosity_mean)
#replace missing free-living morphologies with mean±SD free-living morphology from Husband et al (2021)
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"free")),]$rugosity_mean=1.38
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"free")),]$rugosity_sd=0.29
#replace everything else with mean±SD encrusting morphology from Husband et al (2021)
missing_rugosity[is.na(missing_rugosity$rugosity_mean),]$rugosity_mean=1.64
missing_rugosity[is.na(missing_rugosity$rugosity_sd),]$rugosity_sd=0.52

#merge complete and missing rugosities with filled in rugosity data
rug_rates=rbind(complete_rugosity,missing_rugosity)
rug_rates$rugosity_sd=as.numeric(rug_rates$rugosity_sd)
rug_rates$genera=word(rug_rates$name, 1)

#load in mean colony sizes from the NCRMP data set
NCRMP_genus = rbind(NCRMP_Florida,NCRMP_PuertoRico,NCRMP_USVI) %>% 
  rename("name"=SPECIES_NAME) %>% 
  drop_na(name) %>% 
  mutate(genera=word(name, 1)) %>% 
  group_by(genera) %>% 
  summarize(length_mean=quantile(as.numeric(MAX_DIAMETER), probs = 0.5,na.rm=TRUE),length_lwr=quantile(as.numeric(MAX_DIAMETER), probs = 0.025,na.rm=TRUE),length_upr=quantile(as.numeric(MAX_DIAMETER), probs = 0.975,na.rm=TRUE)) %>% 
  ungroup()

NCRMP_species = rbind(NCRMP_Florida,NCRMP_PuertoRico,NCRMP_USVI) %>% 
  rename("name"=SPECIES_NAME) %>% 
  drop_na(name) %>% 
  group_by(name)  %>%
  summarize(length_mean=quantile(as.numeric(MAX_DIAMETER), probs = 0.5,na.rm=TRUE),length_lwr=quantile(as.numeric(MAX_DIAMETER), probs = 0.025,na.rm=TRUE),length_upr=quantile(as.numeric(MAX_DIAMETER), probs = 0.975,na.rm=TRUE)) %>% 
  ungroup()

#fill in missing species data with genera means
NCRMP_species_rates=merge(rug_rates,NCRMP_species,by="name",all.x=TRUE)

NCRMP_species_rates_na=NCRMP_species_rates[is.na(NCRMP_species_rates$length_mean),] %>% 
  select(-c(length_mean,length_lwr,length_upr))

NCRMP_genera_rates=merge(NCRMP_species_rates_na,NCRMP_genus,by="genera",all.x=TRUE)

NCRMP_sp_gen_rates=rbind(NCRMP_species_rates[!is.na(NCRMP_species_rates$length_mean),],NCRMP_genera_rates)

NCRMP_morphology=NCRMP_sp_gen_rates %>% 
  group_by(morphology.x) %>% 
  summarize(length_mean=mean(length_mean,na.rm=TRUE),length_lwr=mean(length_lwr,na.rm=TRUE),length_upr=mean(length_upr,na.rm=TRUE))

NCRMP_morphology[!is.na(str_match(NCRMP_morphology$morphology.x,"encrusting")),]$length_mean=
  NCRMP_morphology[!is.na(str_match(NCRMP_morphology$morphology.x,"encrusting/sub-massive")),]$length_mean

NCRMP_morphology[!is.na(str_match(NCRMP_morphology$morphology.x,"encrusting")),]$length_lwr=
  NCRMP_morphology[!is.na(str_match(NCRMP_morphology$morphology.x,"encrusting/sub-massive")),]$length_lwr

NCRMP_morphology[!is.na(str_match(NCRMP_morphology$morphology.x,"encrusting")),]$length_upr=
  NCRMP_morphology[!is.na(str_match(NCRMP_morphology$morphology.x,"encrusting/sub-massive")),]$length_upr

NCRMP_morphology_subs=merge(NCRMP_sp_gen_rates[is.na(NCRMP_sp_gen_rates$length_mean),],NCRMP_morphology,by="morphology.x",all.x=TRUE) %>% 
  select(!c(length_mean.x,length_lwr.x,length_upr.x)) %>% 
  rename(length_mean=length_mean.y) %>% 
  rename(length_lwr=length_lwr.y) %>% 
  rename(length_upr=length_upr.y)

NCRMP_lengths=rbind(NCRMP_sp_gen_rates[!is.na(NCRMP_sp_gen_rates$length_mean),],NCRMP_morphology_subs)

#create genera means
genera_rates = NCRMP_lengths %>%
  group_by(genera) %>%
  summarize(coeff_mean=mean(coeff_mean),coeff_lwr=mean(coeff_lwr),coeff_upr=mean(coeff_upr),int_mean=mean(int_mean),int_lwr=mean(int_lwr),int_upr=mean(int_upr),rugosity_mean=mean(rugosity_mean),rugosity_sd=mean(rugosity_sd),length_mean=mean(length_mean),length_lwr=mean(length_lwr),length_upr=mean(length_upr))

#remove macroalgae genera because this only applies to Macroalgae with CCA
genera_rates=genera_rates[is.na(str_match(genera_rates$genera,"Macroalgae")),]
#remove Peysonellid genera to avoid duplicates
genera_rates=genera_rates[is.na(str_match(genera_rates$genera,"Peysonellid")),]

#merge genera means and taxa-specific rates
rug_means=as_tibble(with(NCRMP_lengths,cbind(name,coeff_mean,coeff_lwr,coeff_upr,int_mean,int_lwr,int_upr,rugosity_mean,rugosity_sd,length_mean,length_upr,length_lwr)))
genera_means=as_tibble(with(genera_rates,cbind(name=genera,coeff_mean,coeff_lwr,coeff_upr,int_mean,int_lwr,int_upr,rugosity_mean,rugosity_sd,length_mean,length_upr,length_lwr)))
rug_all=as_tibble(rbind(genera_means,rug_means))
rug_all=transform(rug_all, coeff_mean = as.numeric(coeff_mean),coeff_lwr = as.numeric(coeff_lwr),coeff_upr = as.numeric(coeff_upr),int_mean = as.numeric(int_mean),int_lwr = as.numeric(int_lwr),int_upr = as.numeric(int_upr),rugosity_mean=as.numeric(rugosity_mean),rugosity_sd=as.numeric(rugosity_sd),length_mean=as.numeric(length_mean),length_lwr=as.numeric(length_lwr),length_upr=as.numeric(length_upr))

#correct intercept rounding from ReefBudget import
rug_all$int_mean=round(rug_all$int_mean,digits=4)
rug_all$int_lwr=round(rug_all$int_lwr,digits=4)
rug_all$int_upr=round(rug_all$int_upr,digits=4)

#include microbioerosion rate based on initial logical
microbioerosion_atlantic=ifelse(include_bioerosion==TRUE,-1,0)*as.numeric(microbioerosion[2,5])

#determines mean, lower, and upper rate in kg m^-2 yr^-1 for 100% cover of the respective coral using NCRMP colony sizes and species-specific rugosity from Gonzalez-Barrios
scaled_calcification_rates=NULL
for (k in 1:nrow(rug_all))
{
  nsim=100000
  name=rug_all$name[k]
  s=rtriangle(n=nsim,c=rug_all$length_mean[k],a=rug_all$length_lwr[k],b=rug_all$length_upr[k])
  c=rtriangle(n=nsim,c=rug_all$coeff_mean[k],a=rug_all$coeff_lwr[k],b=rug_all$coeff_upr[k])
  i=rtriangle(n=nsim,c=rug_all$int_mean[k],a=rug_all$int_lwr[k],b=rug_all$int_upr[k])
  r=rnorm(n=nsim,mean=rug_all$rugosity_mean[k],sd=rug_all$rugosity_sd[k])
  n=100/s
  g=(n*(c*s*r+i))/10
  calc=round(quantile(g,probs=0.5,na.rm=TRUE),digits=2)
  calc_lower=round(quantile(g,probs=0.25,na.rm=TRUE),digits=2)
  calc_upper=round(quantile(g,probs=0.75,na.rm=TRUE),digits=2)
  scaled_calcification_rates=rbind(scaled_calcification_rates,data.frame(name,calc,calc_lower,calc_upper))
}

#Modify CoralNet label columns
colnames(CoralNet_Atlantic)=c("name","functional_group","label")

#merge calcification rates with CoralNet Labels
CoralNet_Atlantic_Rates = merge(scaled_calcification_rates,CoralNet_Atlantic,by="name",all=TRUE)

#Correct mismatched names (including typos) with CoralNet to ensure correct mapping of calcification rate data
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Crustose coralline algae")),]$name="CCA (crustose coralline algae)"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Dichocoenia stokesii")),]$name="Dichocoenia stokesi"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"^Hard$")),]$name="Hard coral"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Hard coral \\(branched\\)")),]$name="Hard Coral (branching)"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Hard coral \\(encrusting\\)")),]$name="Hard Coral (encrusting)"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Hard coral \\(massive\\)")),]$name="Hard Coral (massive)"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Hard coral \\(plate/foliose\\)")),]$name="Hard Coral (foliose)"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Mycetophyllia danae")),]$name="Mycetophyllia danaana"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Mycetophyllia lamarckiana")),]$name="Mycetophyllia lamarckana"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Peysonellid")),]$name="Macroalgae: Laminate: red: peysonnelia"

#Import microbioerosion data and apply to "rock" substrate labels
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock")),]$calc=microbioerosion_atlantic
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock")),]$calc_lower=microbioerosion_atlantic
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock")),]$calc_upper=microbioerosion_atlantic

#Import microbioerosion data and apply to parrotfish bite scar substrate labels
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$label,"BiteScar")),]$calc=microbioerosion_atlantic
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$label,"BiteScar")),]$calc_lower=microbioerosion_atlantic
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$label,"BiteScar")),]$calc_upper=microbioerosion_atlantic

#Replace "Rock Crustose Coralline Algae" label with "Crustose Coralline Algae" rate
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock Crustose Coralline Algae")),]$calc=scaled_calcification_rates[!is.na(str_match(scaled_calcification_rates$name,"Crustose coralline algae")),]$calc
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock Crustose Coralline Algae")),]$calc_lower=scaled_calcification_rates[!is.na(str_match(scaled_calcification_rates$name,"Crustose coralline algae")),]$calc_lower
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock Crustose Coralline Algae")),]$calc_upper=scaled_calcification_rates[!is.na(str_match(scaled_calcification_rates$name,"Crustose coralline algae")),]$calc_upper

#Replace "CCA" labels with "Crustose Coralline Algae" rate
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"CCA")),]$calc=scaled_calcification_rates[!is.na(str_match(scaled_calcification_rates$name,"Crustose coralline algae")),]$calc
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"CCA")),]$calc_lower=scaled_calcification_rates[!is.na(str_match(scaled_calcification_rates$name,"Crustose coralline algae")),]$calc_lower
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"CCA")),]$calc_upper=scaled_calcification_rates[!is.na(str_match(scaled_calcification_rates$name,"Crustose coralline algae")),]$calc_upper

#Import macrobioerosion data and apply to boring sponge substrate labels
macrobiorates=macrobioerosion[1:8,c(3,7)]
colnames(macrobiorates)=c("name","rate")
macrobiorates$rate=round(ifelse(include_bioerosion==TRUE,-1,0)*as.numeric(macrobiorates$rate),digits=2)
cliona_delitrix_rate=round(macrobiorates[!is.na(str_match(macrobiorates$name,"delitrix")),]$rate,digits=2)
cliona_mean=round(mean(macrobiorates[!is.na(str_match(macrobiorates$name,"C.")),]$rate),digits=2)
cliona_sd=round(sd(macrobiorates[!is.na(str_match(macrobiorates$name,"C.")),]$rate),digits=2)
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Clio")),]$calc=cliona_mean
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Clio")),]$calc_lower=(cliona_mean-cliona_sd)
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Clio")),]$calc_upper=(cliona_mean+cliona_sd)
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Cliona delitrix")),]$calc=cliona_delitrix_rate
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Cliona delitrix")),]$calc_lower=cliona_delitrix_rate
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Cliona delitrix")),]$calc_upper=cliona_delitrix_rate

#add atlantic region label
CoralNet_Atlantic_Rates$region="Atlantic"
#select only names that match CoralNet, remove extra columns, and drop missing calcification rates
CoralNet_Atlantic_Rates2=merge(CoralNet_Atlantic_Rates,CoralNet_Atlantic,by="name",all.y=TRUE) %>% 
  select('Region'=region,'Name'=name,'Mean'=calc,'Lower bound'=calc_lower,'Upper bound'=calc_upper) %>% 
  drop_na()
#order calcification rates by name
calc_rates_atlantic=CoralNet_Atlantic_Rates2[order(CoralNet_Atlantic_Rates2$Name),] %>% distinct()

file_name=ifelse(include_bioerosion==TRUE,"CoralNet_WesternAtlantic_Calcification_With_Bioerosion_Rates_v2.csv","CoralNet_WesternAtlantic_Calcification_Without_Bioerosion_Rates_v2.csv")
write.csv(calc_rates_atlantic,file_name,row.names=FALSE)
