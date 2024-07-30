# Area-normalized scaling of ReefBudget calcification, macrobioerosion, and microbioerosion rates for use with CoralNet Version 2.0
The following code generates area-normalized calcification, macrobioerosion, and microbioerosion rates from ReefBudget methodologies (http://geography.exeter.ac.uk/reefbudget/; Perry et al., 2018; Perry and Lange 2019) for use with CoralNet image identification labels to support the integration of estimated carbonate production rates with CoralNet, an automated benthic image analysis platform (https://coralnet.ucsd.edu/; Beijbom et al., 2015). Calcification and bioerosion rates have been developed separately for the Indo-Pacific and the Western Atlantic following the below methodologies and data sources. This release serves as an update to the previous Version 1.0 (https://doi.org/10.5281/zenodo.5140477; Courtney et al., 2021).

#### Inclusion of Bioerosion
The default rates include coral calcification, CCA calcification, macrobioerosion, and microbioerosion to be consistent with ReefBudget calcification sheet outputs, but the include_bioerosion argument can be set to FALSE at the beginning of the script to calculate the gross carbonate production without any sources of bioerosion. Because rates have been adapted for benthic image labels, the resulting CoralNet total calcification estimates do not account for mobile sources of parrotfish or urchin bioerosion that are included in the net reef carbonate budget following ReefBudget methodologies. 

#### User defined area-normalized rates
The attached code can be adapted and run for location-specific ReefBudget spreadsheets to allow the user to redefine local calcification rates based on additional data. Alternatively, the included ReefBudget2CoralNet.R function and Shiny app (https://traviscourtney.shinyapps.io/reefbudget2coralnet/) can be used to to generate area-normalized calcification rates for additional labels using linear extension rate and skeletal density (input into ReefBudget calcification sheet to generate calcification coefficient and intercept), average colony size, and colony rugosity.

#### Updates compared to Version 1.0

Version 2.0 includes the following updates relative to Version 1.0 (https://doi.org/10.5281/zenodo.5140477; Courtney et al., 2021):

-	Updated from ReefBudget Indo-Pacific Version v1.2 to v1.3 and ReefBudget Caribbean Version v2.1 to v2.3
-	Removed open-space conversion factor for branching corals because conversion factors are already integrated into the coefficient and intercept estimates derived from ReefBudget
-	Updated colony rugosity estimates for the Indo-Pacific from Gonzalez-Barios & Alvarez-Filip (2019) Caribbean morphology estimates to Husband et al (2022) Indo-Pacific measurements
- Removed microbioerosion and macrobioerosion from coral and CCA substrates when calculating area-normalized calcification rates (note: Indo-Pacific CCA still uses macrobioerosion rates following ReefBudget v1.3)
- Updated Western Atlantic colony sizes from CARICOMP to more recent NCRMP survey data
- Revised Western Atlantic equation according to change in colony sizes from 3D (CARICOMP) to 2D (NCRMP)
- Used Husband et al (2022) to fill in morphologies lacking colony rugosity data for Western Atlantic taxa not included in Gonzalez-Barios & Alvarez-Filip (2019)
-	Replaced ReefBudgetRandCalc with rtriangle for simulating asymmetric uncertainties because ReefBudgetRandCalc was no longer supported by the latest version of R
-	Added seed for generating random numbers and increased number of simulations from 10,000 to 100,000 to increase consistency for generating area-normalized calcification rates
- Matched calcification rate data to additional CoralNet labels

The net effect of these updates on the area-normalized calcification rates from V1.0 to V2.0 varied greatly. Rates for branching increased substantially due to the removal of the open-space conversion factor for branching taxa from the area-normalized rate calculations and also for and columnar taxa due to the removal of the conversion factor from the ReefBudget sheets. Changes in rates for the other morphologies and taxa were typically smaller and resulted from the combined, interacting effects of the above changes. The hard coral mean rates for the Indo-Pacific and Western Atlantic are summarized here to show general changes in the area-normalized calcification rates between v1.0 and v2.0.

Indo-Pacific hard coral mean rates:

- Branching rates increased from 11.64 to 30.31 kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>
- Columnar rates increased from 7.23 to 37.69 kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>
- Free-living rates decreased from 38.61 to 27.68 kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>
- Foliose rates increased slightly from 7.26 to 11.46 kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup> 
- Massive rates remained similar from 22.21 to 20.76 kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>
- Encrusting rates remained similar from 4.09 to 4.21 kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>

Western Atlantic hard coral mean rates:

- Branching rates increased from 1.18 to 7.78 kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>
- Foliose rates decreased from 5.27 to 2.81 kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup> 
- Massive rates remained similar from 11.81 to 11.60 kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>
- Encrusting rates decreased slightly from 2.44 to 1.88 kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>

#### Please cite the following data release for use of the area-normalized calcification and bioerosion rates as follows:

Courtney TA, Lange ID, Sannassy Pilly S, Townsend JE, Chan S, Perry CT, Kriegman DJ, Andersson AJ (2024) Area-normalized scaling of ReefBudget calcification, macrobioerosion, and microbioerosion rates for use with CoralNet Version 2.0

## Estimating CoralNet calcification rates for the Indo-Pacific


Taxa-specific area-normalized calcification rates (G=kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>) were calculated iteratively as G=50th percentile, G<sub>lower</sub>=25th percentile, and G<sub>upper</sub>=75th percentile of a Monte-Carlo simulation (n=100,000) using randomly selected values within the range of uncertainties for each of the taxa-specific equation terms in the below equation:
  
$G = \frac{n(csr+i)}{10}$

n = number of colonies per linear meter (±95%)
	
	Source: 100/s
	
s = median colony diameter (cm) (±95%)
	
	Source: NOAA Pacific Islands coral demography data (see multiple references below)
	
c = calcification rate coefficient (± uncertainties) 
	
	Source: ReefBudget Indo-Pacific v1.3 (Perry et al. 2018)
	
i = calcification rate intercept (± uncertainties) 
	
	Source: ReefBudget Indo-Pacific v1.3 (Perry et al. 2018)
	
r = colony rugosity (± uncertainty)
	
	Source: mean (±SD) colony rugosity for morphology (Husband et al. 2022) 
	        the most common morphology for each genus was selected from NOAA Pacific Islands coral demography data

10 = convert units to kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>

#### Substitutions:

1)	Genus and morphology filled in equation terms for labels at coarser taxa resolution
2)	<i>Porites</i> foliose c and i were determined by substituting extension and density data from ReefBudget <i>Porites</i> into ReefBudget Hard Coral Foliose
3)	<i>Millepora</i> columnar c and i were determined by substituting extension and density data for ReefBudget <i>Millepora</i> into ReefBudget Hard Coral Columnar
4)  Macrobioerosion rates were added to the crustose coralline algae calcification rates following ReefBudget Indo-Pacific v1.3
4)  Macrobioerosion and microbioerosion rates (±uncertainty) from ReefBudget Indo-Pacific v1.3 were applied to rock, bare reef substrate, rubble, and algal covered hard substrate (Note: This assumes all rocks are calcium carbonate. The user should replace the rates for non-carbonate rocks with 0 or NA.)

## Estimating CoralNet calcification rates for the Western Atlantic


Species-specific area-normalized calcification rates (G=kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>) were calculated iteratively as G=50th percentile, G<sub>lower</sub>=25th percentile, and G<sub>upper</sub>=75th percentile of a Monte-Carlo simulation (n=100,000) using randomly selected values within the range of uncertainties for each of the taxa-specific equation terms in the below equation:
  
$G = \frac{n(csr+i)}{10}$
  
n = number of colonies per linear meter (±95%)
	
	Source: 100/s
	
s = median colony diameter (cm) (±95%)
	
	Source: NOAA Western Atlantic coral demography data (see multiple references below)

c = calcification rate coefficient (± uncertainties) 

	Source: ReefBudget Caribbean v2.3 (Perry and Lange 2019)
	
i = calcification rate intercept (± uncertainties) 

	Source: ReefBudget Caribbean v2.3 (Perry and Lange 2019)
	
r = colony rugosity (± uncertainty)

	Source: mean (±SD) rugosity for each species with substitutions by morphology (González-Barrios and Álvarez-Filip 2018)
	        mean (±SD) rugosity for free-living and encrusting morphologies (Husband et al. 2022)
	
10 = convert units to kg CaCO<sub>3</sub> m<sup>-2</sup> yr<sup>-1</sup>

#### Substitutions:

1)	Genus and morphology filled in equation terms for labels at coarser taxa resolution
2)	The mean r of all morphologies was used to fill in encrusting morphologies, assuming that encrustation occurs over mean reef structural complexity
3)  Microbioerosion rates from ReefBudget Caribbean v2.1 were applied to rock, bare reef substrate, rubble, algal covered hard substrate, and parrotfish bite scars (Note: This assumes all rocks are calcium carbonate. The user should replace the rates for non-carbonate rocks with 0 or NA.)
4)  <i>Cliona delitrix</i> bioerosion rates from ReefBudget Caribbean v2.1 were applied to <i>Cliona delitrix</i> and the mean (±SD) rate from ReefBudget Caribbean v2.1 of all Clionid sponges were applied to cover of Clionid sponges

### References:

Beijbom O, Edmunds PJ, Roelfsema C, Smith J, Kline DI, Neal BP, Dunlap MJ, Moriarty V, Fan T-Y, Tan C-J, Chan S, Treibitz T, Gamst A, Mitchell BG, Kriegman D (2015) Towards Automated Annotation of Benthic Survey Images: Variability of Human Experts and Operational Modes of Automation. PLoS ONE 10(7): e0130312.

Courtney, T. A., Chan, S., Lange, I. D., Perry, C. T., Kriegman, D. J., & Andersson, A. J. (2021). Area-normalized scaling of ReefBudget calcification, macrobioerosion, and microbioerosion rates for use with CoralNet Version 1.0 (1.0). Zenodo. https://doi.org/10.5281/zenodo.5140477

González-Barrios FJ, Álvarez-Filip L. A framework for measuring coral species-specific contribution to reef functioning in the Caribbean. Ecological Indicators. 2018 Dec 1;95:877-86.

Husband E, Perry CT, Lange ID. Estimating rates of coral carbonate production from aerial and archive imagery by applying colony scale conversion metrics. 2022 April 14;41:1199-1209.

Perry CT and Lange ID (2019) ReefBudget Caribbean v2: online resource and methodology. Retrieved from http://geography.exeter.ac.uk/reefbudget/

Perry CT, Lange I, Januchowski-Hartley FA (2018) ReefBudget Indo Pacific: online resource and methodology. Retrieved from http://geography.exeter.ac.uk/reefbudget/

#### NOAA Pacific Islands Coral Demography Data:

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across American Samoa in 2015 (NCEI Accession 0159173). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0159173.

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across Jarvis Island from 2016-05-16 to 2016-05-22 (NCEI Accession 0159164). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0159164.

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across the Main Hawaiian Islands from 2013-08-02 to 2013-10-29 (NCEI Accession 0159147). NOAA National Centers for Environmental Information. Dataset. https://accession.nodc.noaa.gov/0159147.

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across the Mariana Archipelago from 2017-05-05 to 2017-06-21 (NCEI Accession 0166383). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0166383.

Ecosystem Sciences Division, Pacific Islands Fisheries Science Center; Papahānaumokuākea Marine National Monument (2019). National Coral Reef Monitoring Program: Stratified random surveys (StRS) of coral demography (adult and juvenile corals) at French Frigate Shoals, Lisianski Island, and Midway Atoll of the Northwestern Hawaiian Islands from 2014-08-14 to 2014-08-26 (NCEI Accession 0184908). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0184908. 

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across the Pacific Remote Island Areas from 2015-01-26 to 2015-04-28 (NCEI Accession 0159161). NOAA National Centers for Environmental Information. Dataset. https://accession.nodc.noaa.gov/0159161.

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across Rose Atoll from 2016-05-01 to 2016-05-04 (NCEI Accession 0159167). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0159167. 

Coral Reef Ecosystem Program; Pacific Islands Fisheries Science Center (2017). National Coral Reef Monitoring Program: Stratified Random Surveys (StRS) of Coral Demography (Adult and Juvenile Corals) across Wake Island from 2014-03-16 to 2014-03-20 (NCEI Accession 0159162). NOAA National Centers for Environmental Information. Dataset. https://www.ncei.noaa.gov/archive/accession/0159162.

#### NOAA Western Atlantic Coral Demography Data:

NOAA National Centers for Coastal Ocean Science (2018). National Coral Reef Monitoring Program: Assessment of coral reef benthic communities in the U.S. Virgin Islands. [indicate subset used]. NOAA National Centers for Environmental Information. Dataset. https://doi.org/10.7289/v5ww7fqk.

NOAA National Centers for Coastal Ocean Science; NOAA Southeast Fisheries Science Center (2018). National Coral Reef Monitoring Program: Assessment of coral reef benthic communities in Puerto Rico. [indicate subset used]. NOAA National Centers for Environmental Information. Dataset. https://doi.org/10.7289/v5pg1q23.

NOAA Southeast Fisheries Science Center; NOAA National Centers for Coastal Ocean Science (2018). National Coral Reef Monitoring Program: Assessment of coral reef benthic communities in the Florida Reef Tract. [indicate subset used]. NOAA National Centers for Environmental Information. Dataset. https://doi.org/10.7289/v5xw4h4z.
