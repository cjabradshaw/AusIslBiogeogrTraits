# Island biogeography predicts vertebrate trait space
<img align="right" src="www/aus islands map.jpg" alt="islands of Australia" width="350" style="margin-bottom: 20px">

<br>
Corey Bradshaw<br>
<a href="http://globalecologyflinders.com">Global Ecology</a>, Flinders University, Australia<br>
<a href="mailto:corey.bradshaw@flinders.edu.au">e-mail</a><br>
<br>
John Llewelyn<br>
<a href="http://globalecologyflinders.com">Global Ecology</a>, Flinders University, Australia<br>
<a href="mailto:john.llewelyn@flinders.edu.au">e-mail</a><br>
<br>
Accompanies paper:<br>
<br>
Bradshaw, CJA, A Naglis, F Saltré, C Mudge, J Llewelyn. Island biogeography similarly predicts species and trait diversity of vertebrate communities on oceanic islands (in preparation)<br>

## Scripts
- <code>islbiogeogrtraits.R</code>: R code for all analyses

## Data
### GIS layers
- <em>ausislavala.zip</em>: Atlas of Living Australia data for bird species (shapefile). Due to Github file-size constraints, we have broken the zip file into similar-sized chunks. Combine chunks aa to ah using the following Terminal (or equivalent) command: <code>cat ausislavalazip_chunk_* > ausislavala.zip</code>. Once chunks are recombined, unzip the corresponding shapefile.
- <em>reptislavala.zip</em>: Atlas of Living Australia data for reptile species (shapefile).
- <em>mamaislavala.zip</em>: Atlas of Living Australia data for mammal species (shapefile).
- <em>amphislavala.zip</em>: Atlas of Living Australia data for amphibians species (shapefile).
- <em>ausislands.zip</em>: Australian islands shapefile (Australian Albers Equal Area projection).
- <em>ausislandsAlbersGeom_corrected.zip: Australian islands shapefile  (Australian Albers Equal Area projection) after removing ~ 30% of islands following criteria outlined in the main text.
- <em>mainlandsLL.zip</em>: mainlands (Australia and surrounding continents) in lat/lon datum.

### Trait data
#### Birds
- <em>AVONET1_Birdlife.xlsx</em>: <a href="https://doi.org/10.1111/ele.13898">AVONET</a> BirdLife traits dataset
- <em>AVONET2_eBird.xlsx</em>: <a href="https://doi.org/10.1111/ele.13898">AVONET</a> eBird traits dataset
- <em>AVONET3_BirdTree.xlsx</em>: <a href="https://doi.org/10.1111/ele.13898">AVONET</a> BirdTree traits dataset

#### Mammals
- <em>SahulTraitsMam.csv</em>: SahulTraits mammal traits database

## Required R libraries
<code>dismo</code>, <code>dplyr</code>, <code>galah</code>, <code>gawdis</code>, <code>gbm</code>, <code>ggplot2</code>, <code>gridExtra</code>, <code>mFD</code>, <code>readxl</code>,<code>sf</code>, <code>terra</code>
<br>
<br>
<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" width="180" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="100" style="margin-top: 20px"></a>  &nbsp; &nbsp; <a href="https://ciehf.au"><img align="bottom-left" src="www/CIEHFlogoHorizTransp.png" alt="CIEHF logo" width="200" style="margin-top: 20px"></a>  &nbsp; &nbsp; &nbsp; <a href="https://www.epicaustralia.org.au"><img align="bottom-left" src="www/CabahFCL.jpg" alt="CABAH logo" width="180" style="margin-top: 20px"></a></p>
