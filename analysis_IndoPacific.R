######################################################################
##’ @title Area-normalized scaling of Indo-Pacific ReefBudget calcification and bioerosion rates for use with CoralNet
##’
##’ @author Travis A Courtney, PhD
##’ @contact travis.courtney@upr.edu
##’ @date 2024-06-25
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

#########CoralNet Indo-Pacific Calcification Rates#########
#Import Required Data
setwd("Data_IndoPacific")
setwd("NOAA_Pacific_Coral_Demography") #Import NOAA Pacific Coral Demography Data from multiple files
NOAA_demography_raw=
  dir_ls(regexp = "\\.csv$") %>% 
  map_dfr(read_csv, .id = "source")
setwd("..")
NOAA_CRED_CoralNet_labelset <- read_csv("NOAA_CRED_CoralNet_labelset2.csv") #Import NOAA CRED CoralNet labelset
Husband_et_al_rugosity <- read_excel("338_2022_2247_MOESM2_ESM.xlsx") #Husband et al estimates of colony structural complexity
ReefBudget_Pacific <- read_excel("Indo-Pacific_Carbonate_Production_v1.3-Jan23.xlsx", sheet = "Calcification Rates", skip = 9) #import Indo-Pacific ReefBudget rates
macromicro_bioerosion=read_excel("Indo-Pacific_Carbonate_Production_v1.3-Jan23.xlsx", sheet = "Macro- & Microbioerosion") #import Indo_Pacific ReefBudget bioerosion rates
duplicates=read_csv("CoralNetDuplicates.csv") #import coralnet duplicate files
setwd("..")

#remove missing data points
NOAA_demography = subset(NOAA_demography_raw,COLONYLENGTH!=-9)

#re-assign labels to match CoralNet ID's
NOAA_demography$MORPHOLOGY=str_replace(NOAA_demography$MORPHOLOGY,"Mounding","Massive")
NOAA_demography$MORPHOLOGY=str_replace(NOAA_demography$MORPHOLOGY,"Encrusting \\(columnar\\)","Encrusting")
NOAA_demography$MORPHOLOGY=str_replace(NOAA_demography$MORPHOLOGY,"Encrusting \\(mounding\\)","Encrusting")
NOAA_demography$MORPHOLOGY=str_replace(NOAA_demography$MORPHOLOGY,"Encrusting \\(flat\\)","Encrusting")
NOAA_demography$MORPHOLOGY=str_replace(NOAA_demography$MORPHOLOGY,"Disc - free living","Free-living")
NOAA_demography$MORPHOLOGY=tolower(NOAA_demography$MORPHOLOGY)
unique(NOAA_demography$MORPHOLOGY)

#find mean colony lengths for each genus
NOAA_genus_lengths=
  NOAA_demography %>%
  group_by(GENUS) %>%
  summarise(length_mean=quantile(COLONYLENGTH, probs = 0.5, na.rm=TRUE),length_lwr=quantile(COLONYLENGTH, probs = 0.025, na.rm=TRUE),length_upr=quantile(COLONYLENGTH, probs = 0.975, na.rm=TRUE))
NOAA_genus_lengths=NOAA_genus_lengths[!is.na(NOAA_genus_lengths$GENUS),]

#find mean colony lengths for each genus, morphology
NOAA_genus_morphology_lengths=
  NOAA_demography %>%
  group_by(GENUS,MORPHOLOGY) %>%
  summarise(length_mean=quantile(COLONYLENGTH, probs = 0.5, na.rm=TRUE),length_lwr=quantile(COLONYLENGTH, probs = 0.025, na.rm=TRUE),length_upr=quantile(COLONYLENGTH, probs = 0.975, na.rm=TRUE))
NOAA_genus_morphology_lengths=NOAA_genus_morphology_lengths[!is.na(NOAA_genus_morphology_lengths$GENUS),]

#find mean colony lengths for each morphology
NOAA_morphology_lengths=
  NOAA_demography %>%
  group_by(MORPHOLOGY) %>%
  summarise(length_mean=quantile(COLONYLENGTH, probs = 0.5, na.rm=TRUE),length_lwr=quantile(COLONYLENGTH, probs = 0.025, na.rm=TRUE),length_upr=quantile(COLONYLENGTH, probs = 0.975, na.rm=TRUE))
NOAA_morphology_lengths=NOAA_morphology_lengths[!is.na(NOAA_morphology_lengths$MORPHOLOGY),]

#Import NOAA CRED labels
NOAA_CRED_CoralNet_labelset$TAXA = substr(NOAA_CRED_CoralNet_labelset$Name, start = 6, stop=1000)
NOAA_CRED_CoralNet_labelset$GENUS = word(NOAA_CRED_CoralNet_labelset$TAXA,1)
NOAA_CRED_CoralNet_labelset$MORPHOLOGY = substr(word(NOAA_CRED_CoralNet_labelset$TAXA,2),start=5,stop=1000)
colnames(NOAA_CRED_CoralNet_labelset)=c("NAME","CODE","GROUP","TAXA","GENUS","MORPHOLOGY")
NOAA_CRED_morph_lengths=merge(subset(NOAA_CRED_CoralNet_labelset,GROUP=="Hard coral"), NOAA_genus_morphology_lengths, by=c("GENUS","MORPHOLOGY"),all.x=TRUE)

#fill in missing species data with genera means
NOAA_CRED_morph_lengths_all=NOAA_CRED_morph_lengths[!is.na(NOAA_CRED_morph_lengths$length_mean),]
NOAA_CRED_morph_lengths_na=NOAA_CRED_morph_lengths[is.na(NOAA_CRED_morph_lengths$length_mean),]
NOAA_CRED_morph_genus_lengths=merge(NOAA_CRED_morph_lengths_na[,c(1:6)],NOAA_genus_lengths,by="GENUS",all.x=TRUE)

#use consistent labels throughout
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Foliose",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="foliose",]
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Branching",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="branching",]
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Columnar",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="columnar",]
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Encrusting",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="encrusting",]
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Free-living",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="free-living",]
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Massive",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="massive",]
#count whether Goniopora or Alveopora has more records and use most frequent genus
nrow(subset(NOAA_demography,GENUS=="Goniopora")) #199 colonies
nrow(subset(NOAA_demography,GENUS=="Alveopora")) #1 colony
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Goniopora/Alveopora",c(7:9)]=NOAA_genus_lengths[NOAA_genus_lengths$GENUS=="Goniopora",c(2:4)] #Goniopora is more common so use Goniopora
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Diploastrea",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="massive",] #assign massive length to Diploastrea

#find most common morphology for each genus
NOAA_typical_morphology=
  NOAA_demography %>%
  group_by(GENUS) %>%
  summarise(MORPHOLOGY=names(which.max(table(MORPHOLOGY))))

#pair most common morphology for each genera with CRED labels
NOAA_CRED_morph_genus_lengths_filled=NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$MORPHOLOGY != "",]
NOAA_CRED_morph_genus_lengths_na=NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$MORPHOLOGY == "",c(1,3:9)]
NOAA_CRED_morph_genus_lengths_typical=merge(NOAA_CRED_morph_genus_lengths_na,NOAA_typical_morphology,by="GENUS",all.x=TRUE)

#combine all morphology and length data
NOAA_CRED_morph_genus_lengths_all=rbind(NOAA_CRED_morph_lengths_all,NOAA_CRED_morph_genus_lengths_filled,NOAA_CRED_morph_genus_lengths_typical)

#substitute typical morphologies for Goniopora
NOAA_CRED_morph_genus_lengths_all[NOAA_CRED_morph_genus_lengths_all$GENUS=="Goniopora/Alveopora",2]="encrusting" #Goniopora is more common so use Goniopora

#extract rugosity ± sd for each genus and morphology from Husband et al. (2022) data
H_rugosities=bind_cols(
  GENUS=word(Husband_et_al_rugosity$"Genus Morphotype", 1),
  MORPHOLOGY=str_to_lower(word(Husband_et_al_rugosity$"Genus Morphotype", 2)),
  rugosity = Husband_et_al_rugosity$"Mean Rugosity",
  rugosity_sd = Husband_et_al_rugosity$"Std. Dev.")

#add rugosities to dataframe
NOAA_mgr_all=merge(NOAA_CRED_morph_genus_lengths_all, H_rugosities, by=c("GENUS","MORPHOLOGY"),all.x=TRUE)

#replace acropora rugosity data with mean of all Acropora branching categories
NOAA_mgr_all[NOAA_mgr_all$TAXA=="Acropora spp_branching",10]=mean(subset(H_rugosities,GENUS=="Acropora"&MORPHOLOGY!="flat")$rugosity)
NOAA_mgr_all[NOAA_mgr_all$TAXA=="Acropora spp_branching",11]=sd(subset(H_rugosities,GENUS=="Acropora"&MORPHOLOGY!="flat")$rugosity)
NOAA_mgr_all[NOAA_mgr_all$TAXA=="Pocillopora spp",10:11]=subset(H_rugosities,GENUS=="Pocillopora"&MORPHOLOGY=="short,")[,3:4]
NOAA_mgr_all[NOAA_mgr_all$TAXA=="Porites spp_branching",10]=mean(subset(H_rugosities,GENUS=="Porites"&MORPHOLOGY!="massive"&MORPHOLOGY!="unknown")$rugosity)
NOAA_mgr_all[NOAA_mgr_all$TAXA=="Porites spp_branching",11]=sd(subset(H_rugosities,GENUS=="Porites"&MORPHOLOGY!="massive"&MORPHOLOGY!="unknown")$rugosity)

#fill in missing rugosity data
NOAA_mgr_all[NOAA_mgr_all$TAXA=="Acropora spp_tabulate",10:11]=subset(H_rugosities,GENUS=="Acropora"&MORPHOLOGY=="plates")[,3:4]
NOAA_mgr_all[NOAA_mgr_all$TAXA=="Heliopora spp",10:11]=subset(H_rugosities,GENUS=="Heliopora"&MORPHOLOGY=="branching*")[,3:4]
NOAA_mgr_all[NOAA_mgr_all$TAXA=="Pavona spp_foliose",10:11]=subset(H_rugosities,GENUS=="Pavona"&MORPHOLOGY=="")[,3:4]
NOAA_mgr_all[NOAA_mgr_all$TAXA=="Symphyllia spp",10:11]=subset(H_rugosities,GENUS=="Symphyllia")[,3:4]

#Find missing rugosity data to fill in gaps
NOAA_mgr_rugosity=NOAA_mgr_all[!is.na(NOAA_mgr_all$rugosity),]
NOAA_mgr_no_rugosity=NOAA_mgr_all[is.na(NOAA_mgr_all$rugosity),]

NOAA_mgr_no_rugosity[NOAA_mgr_no_rugosity$MORPHOLOGY=="encrusting",10:11]=subset(H_rugosities,GENUS=="Encrusting")[,3:4]
NOAA_mgr_no_rugosity[NOAA_mgr_no_rugosity$MORPHOLOGY=="branching",10:11]=subset(H_rugosities,GENUS=="Branching")[,3:4]
NOAA_mgr_no_rugosity[NOAA_mgr_no_rugosity$MORPHOLOGY=="columnar",10:11]=subset(H_rugosities,GENUS=="Columnar")[,3:4]
NOAA_mgr_no_rugosity[NOAA_mgr_no_rugosity$MORPHOLOGY=="free-living",10:11]=subset(H_rugosities,GENUS=="Free-living")[,3:4]
NOAA_mgr_no_rugosity[NOAA_mgr_no_rugosity$MORPHOLOGY=="foliose",10:11]=subset(H_rugosities,GENUS=="Foliose")[,3:4]
NOAA_mgr_no_rugosity[NOAA_mgr_no_rugosity$MORPHOLOGY=="massive",10:11]=subset(H_rugosities,GENUS=="Massive")[,3:4]
NOAA_mgr_no_rugosity[NOAA_mgr_no_rugosity$MORPHOLOGY=="plating",10:11]=subset(H_rugosities,GENUS=="Plates")[,3:4]

#Merge all rugosity data after gap-filling
NOAA_mgr=bind_rows(NOAA_mgr_rugosity,NOAA_mgr_no_rugosity)

#Import ReefBudget Indo-Pacific data
colnames(ReefBudget_Pacific)=c("code","genera","morphology","ext_mean","ext_sd","dens_mean","dens_sd","cf","coeff_mean","coeff_lwr","coeff_upr","int_mean","int_lwr","int_upr","notes")
ReefBudget_Pacific=ReefBudget_Pacific[!is.na(ReefBudget_Pacific$coeff_mean),]

NOAA_mgr$GENUS_MORPHOLOGY=paste(NOAA_mgr$GENUS,NOAA_mgr$MORPHOLOGY) #combine genus and morphology labels for consistent use
ReefBudget_Pacific$GENUS_MORPHOLOGY=paste(ReefBudget_Pacific$genera,ReefBudget_Pacific$morphology) #combine genus and morphology labels for consistent use
ReefBudget_Pacific_short = as.data.frame(ReefBudget_Pacific[,c(16,9:14)]) #select only columns needed

#Merge NOAA size and rugosity data with ReefBudget rates
NOAA_mgr_rates= merge(NOAA_mgr,ReefBudget_Pacific_short,by=c("GENUS_MORPHOLOGY"),all.x=TRUE)

#use consistent labels throughout
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Acropora tabulate",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Acropora table",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Branching branching",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Hard coral branching",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Columnar columnar",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Hard coral columnar",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Foliose foliose",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Hard coral foliose",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Encrusting encrusting",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Hard coral encrusting",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Massive massive",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Hard coral massive",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Montastraea encrusting",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Montastraea encrusting",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Fungia free-living",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Fungia freeliving",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Goniopora/Alveopora encrusting",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Goniopora encrusting",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Merulina foliose",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Merulina plating",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Symphyllia massive lobate",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Symphyllia massive",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Turbinaria foliose",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Turbinaria plating",2:7]

#create mean rates for all freeliving corals
ReefBudget_Pacific_freeliving_rates =
  subset(ReefBudget_Pacific,morphology=="freeliving") %>% #take all ReefBudget data
  summarize(coeff_mean=mean(coeff_mean),
            coeff_lwr=mean(coeff_lwr),
            coeff_upr=mean(coeff_upr),
            int_mean=mean(int_mean),
            int_lwr=mean(int_lwr), 
            int_upr=mean(int_upr)) 
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Free-living free-living",13:18]=ReefBudget_Pacific_freeliving_rates

#create mean rates for Acropora branching
ReefBudget_Pacific_acropora_branching_rates =
  ReefBudget_Pacific %>% 
  subset(genera=="Acropora") %>% 
  subset(morphology=="arborescent"|morphology=="caespitose"|morphology=="digitate"|morphology=="hispidose") %>% 
  summarize(coeff_mean=mean(coeff_mean),
            coeff_lwr=mean(coeff_lwr),
            coeff_upr=mean(coeff_upr),
            int_mean=mean(int_mean),
            int_lwr=mean(int_lwr), 
            int_upr=mean(int_upr)) 
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Acropora branching",13:18]=ReefBudget_Pacific_acropora_branching_rates

#create mean rates for Plerogyra massive instead of Plerogyra encrusting
ReefBudget_Pacific_plerogyra_branching_rates =
  ReefBudget_Pacific %>% 
  subset(genera=="Plerogyra") %>% 
  subset(morphology==c("massive")) %>% 
  summarize(coeff_mean=mean(coeff_mean),
            coeff_lwr=mean(coeff_lwr),
            coeff_upr=mean(coeff_upr),
            int_mean=mean(int_mean),
            int_lwr=mean(int_lwr), 
            int_upr=mean(int_upr)) 
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Plerogyra encrusting",13:18]=ReefBudget_Pacific_plerogyra_branching_rates

#fill in missing ReefBudget rates with substitutions
Millepora_columnar=c(coeff_mean=0.8174,coeff_lwr=0.4948,coeff_upr=1.1399,int_mean=0,int_lwr=0,int_upr=0) #substitute Millepora extension and density into "Hard Coral Columnar" cells in ReefBudget sheet
Porites_foliose=c(coeff_mean=0.1795,coeff_lwr=0.1112,coeff_upr=0.2554,int_mean=3.9329,int_lwr=2.3863,int_upr=5.7092) #substitute mean±SD of all Porites extension and density into "Hard Coral Foliose" cells in ReefBudget sheet
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Millepora columnar",13:18]=Millepora_columnar
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Porites foliose",13:18]=Porites_foliose

#select only necessary columns
NOAA_coral_rates=NOAA_mgr_rates[,4:18]

#add in CCA, hard substrate, sponges, other growth forms besides CRED labels?
NOAA_CCA_rates=cbind(NOAA_CRED_CoralNet_labelset[str_detect(NOAA_CRED_CoralNet_labelset$NAME,"CCA"),1:4],length_mean=1,length_lwr=1,length_upr=1,NOAA_mgr[NOAA_mgr$TAXA=="Encrusting hard coral",10:11],ReefBudget_Pacific[str_detect(ReefBudget_Pacific$code,"CCA"),9:14])
NOAA_all_rates=rbind(NOAA_coral_rates,NOAA_CCA_rates)

#correct intercept rounding from ReefBudget import
NOAA_all_rates$int_mean=round(NOAA_all_rates$int_mean,digits=4)
NOAA_all_rates$int_lwr=round(NOAA_all_rates$int_lwr,digits=4)
NOAA_all_rates$int_upr=round(NOAA_all_rates$int_upr,digits=4)

#extract macro and microbioerosion rates from data entry sheets
macrobioerosion_pacific=ifelse(include_bioerosion==TRUE,-1,0)*as.numeric(macromicro_bioerosion[2,5])
macrobioerosion_pacific_uncertainty=ifelse(include_bioerosion==TRUE,-1,0)*as.numeric(macromicro_bioerosion[3,5])
microbioerosion_pacific=ifelse(include_bioerosion==TRUE,-1,0)*as.numeric(macromicro_bioerosion[20,5])
microbioerosion_pacific_uncertainty=ifelse(include_bioerosion==TRUE,-1,0)*as.numeric(macromicro_bioerosion[21,5])

#determines mean, lower, and upper rate in kg m^-2 yr^-1 for 100% cover of the respective coral using CARICOMP colony sizes and species-specific rugosity from Husband et al
NOAA_calc_rates=NULL
for (k in 1:nrow(NOAA_all_rates))
{
  nsim=100000
  NAME=NOAA_all_rates$NAME[k]
  CODE=NOAA_all_rates$CODE[k]
  GROUP=NOAA_all_rates$GROUP[k]
  TAXA=NOAA_all_rates$TAXA[k]
  s=rtriangle(n=nsim,c=NOAA_all_rates$length_mean[k],a=NOAA_all_rates$length_lwr[k],b=NOAA_all_rates$length_upr[k])
  c=rtriangle(n=nsim,c=NOAA_all_rates$coeff_mean[k],a=NOAA_all_rates$coeff_lwr[k],b=NOAA_all_rates$coeff_upr[k])
  i=rtriangle(n=nsim,c=NOAA_all_rates$int_mean[k],a=NOAA_all_rates$int_lwr[k],b=NOAA_all_rates$int_upr[k])
  r=rnorm(n=nsim,mean=NOAA_all_rates$rugosity[k],sd=NOAA_all_rates$rugosity_sd[k])
  n=100/s
  g=(n*(c*s*r+i))/10
  calc=round(quantile(g,probs=0.5,na.rm=TRUE),digits=2)
  calc_lower=round(quantile(g,probs=0.25,na.rm=TRUE),digits=2)
  calc_upper=round(quantile(g,probs=0.75,na.rm=TRUE),digits=2)
  NOAA_calc_rates=rbind(NOAA_calc_rates,data.frame(NAME,CODE,GROUP,TAXA,calc,calc_lower,calc_upper))
}

#subtract macrobioerosion from CCA rates
NOAA_calc_rates[!is.na(str_match(NOAA_calc_rates$NAME,"CCA")),]$calc=NOAA_calc_rates[!is.na(str_match(NOAA_calc_rates$NAME,"CCA")),]$calc+macrobioerosion_pacific
NOAA_calc_rates[!is.na(str_match(NOAA_calc_rates$NAME,"CCA")),]$calc_lower=NOAA_calc_rates[!is.na(str_match(NOAA_calc_rates$NAME,"CCA")),]$calc_lower+(macrobioerosion_pacific+macrobioerosion_pacific_uncertainty)
NOAA_calc_rates[!is.na(str_match(NOAA_calc_rates$NAME,"CCA")),]$calc_upper=NOAA_calc_rates[!is.na(str_match(NOAA_calc_rates$NAME,"CCA")),]$calc_upper+(macrobioerosion_pacific-macrobioerosion_pacific_uncertainty)

#find hard substrate that does not include CCA
NOAA_CRED_hard_substrate=NOAA_CRED_CoralNet_labelset[str_detect(NOAA_CRED_CoralNet_labelset$GROUP,"Hard Substrate"),]
NOAA_CRED_hard_nonCCA=NOAA_CRED_hard_substrate[!str_detect(NOAA_CRED_hard_substrate$NAME,"CCA"),]

NOAA_rock_rates=NULL
for (k in 1:nrow(NOAA_CRED_hard_nonCCA))
{
  nsim=100
  NAME=NOAA_CRED_hard_nonCCA$NAME[k]
  CODE=NOAA_CRED_hard_nonCCA$CODE[k]
  GROUP=NOAA_CRED_hard_nonCCA$GROUP[k]
  TAXA=NOAA_CRED_hard_nonCCA$TAXA[k]
  macrob=runif(nsim,macrobioerosion_pacific+macrobioerosion_pacific_uncertainty,macrobioerosion_pacific-macrobioerosion_pacific_uncertainty)
  microb=runif(nsim,microbioerosion_pacific+microbioerosion_pacific_uncertainty,microbioerosion_pacific-microbioerosion_pacific_uncertainty)
  g=macrob+microb
  calc=round(quantile(g,probs=0.5),digits=2)
  calc_lower=round(quantile(g,probs=0.25),digits=2)
  calc_upper=round(quantile(g,probs=0.75),digits=2)
  NOAA_rock_rates=rbind(NOAA_rock_rates,data.frame(NAME,CODE,GROUP,TAXA,calc,calc_lower,calc_upper))
}

#assign NA values to all other labels
NOAA_CRED_other_label_rates=cbind(subset(NOAA_CRED_CoralNet_labelset,GROUP!="Hard coral"&GROUP!="Hard Substrate")[,1:4],calc=NA,calc_lower=NA,calc_upper=NA)

#select only necessary columns
NOAA_CRED_all=rbind(NOAA_calc_rates,NOAA_rock_rates,NOAA_CRED_other_label_rates)
NOAA_CRED_all$region="Indo-Pacific"
colnames(NOAA_CRED_all)=c("name","label","group","taxa","calc","calc_lower","calc_upper","region")
NOAA_CRED_all$label=str_remove(NOAA_CRED_all$label, "[*]")
calc_rates_indopacific=bind_cols('Region'=NOAA_CRED_all$region,'Name'=NOAA_CRED_all$name,'Mean'=NOAA_CRED_all$calc,'Lower bound'=NOAA_CRED_all$calc_lower,'Upper bound'=NOAA_CRED_all$calc_upper)

#import duplicate labels for NOAA CRED to general labels
colnames(duplicates)=c("Name","Duplicate")
calc_rates_indopacific_duplicates=merge(calc_rates_indopacific,duplicates,by="Name")
duplicate_labels_indopacific=calc_rates_indopacific_duplicates%>%select('Region','Duplicate','Mean','Lower bound','Upper bound')
duplicate_labels_indopacific=rename(duplicate_labels_indopacific,Name=Duplicate)

#generate genus level and hard coral means from genus+morphology and hard coral+morphology NOAA CRED labels 
Acropora=calc_rates_indopacific[str_which(calc_rates_indopacific$Name,"Acropora"),]
Acropora_rates=bind_cols(Region="Indo-Pacific",Name="Acropora",'Mean'=round(mean(Acropora$'Mean'),digits=2),'Lower bound'=round(mean(Acropora$'Lower bound'),digits=2),'Upper bound'=round(mean(Acropora$'Upper bound'),digits=2))
Montipora=calc_rates_indopacific[str_which(calc_rates_indopacific$Name,"Montipora"),]
Montipora_rates=bind_cols(Region="Indo-Pacific",Name="Montipora",'Mean'=round(mean(Montipora$Mean),digits=2),'Lower bound'=round(mean(Montipora$'Lower bound'),digits=2),'Upper bound'=round(mean(Montipora$'Upper bound'),digits=2))
Porites=calc_rates_indopacific[str_which(calc_rates_indopacific$Name,"Porites"),]
Porites_rates=bind_cols(Region="Indo-Pacific",Name="Porites",'Mean'=round(mean(Porites$Mean),digits=2),'Lower bound'=round(mean(Porites$'Lower bound'),digits=2),'Upper bound'=round(mean(Porites$'Upper bound'),digits=2))
Hard_coral=calc_rates_indopacific[str_which(calc_rates_indopacific$Name,"hard coral"),]
Hard_coral_rates=bind_cols(Region="Indo-Pacific",Name="Hard coral",'Mean'=round(mean(Hard_coral$Mean),digits=2),'Lower bound'=round(mean(Hard_coral$'Lower bound'),digits=2),'Upper bound'=round(mean(Hard_coral$'Upper bound'),digits=2))

#manually pair up CRED CCA and Rock duplicates with commonly used CoralNet labels
CCA_rates=calc_rates_indopacific[str_which(calc_rates_indopacific$Name,"CRED-CCA growing on hard substrate"),]
CCA_rates$Name="CCA (crustose coralline algae)"
Rock_rates=calc_rates_indopacific[str_which(calc_rates_indopacific$Name,"CRED-Hard substrate"),]
Rock_rates$Name="Rock"
Bare_Rock_rates=calc_rates_indopacific[str_which(calc_rates_indopacific$Name,"CRED-Hard substrate"),]
Bare_Rock_rates$Name="Bare Rock"

#combine NOAA CRED labels with duplicate labels and genus level labels
calc_rates_indopacific_all=bind_rows(calc_rates_indopacific,duplicate_labels_indopacific,Acropora_rates,Montipora_rates,Porites_rates,Hard_coral_rates,CCA_rates,Rock_rates,Bare_Rock_rates) %>% arrange(Name)
calc_rates_indopacific_all=calc_rates_indopacific_all[complete.cases(calc_rates_indopacific_all),]

file_name=ifelse(include_bioerosion==TRUE,"CoralNet_IndoPacific_Calcification_With_Bioerosion_Rates_v2.csv","CoralNet_IndoPacific_Calcification_Without_Bioerosion_Rates_v2.csv")
write.csv(calc_rates_indopacific_all,file_name,row.names=FALSE)
