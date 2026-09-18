# CETRACES: A molecular toolkit for enhancing cetacean monitoring in the Bay of Biscay.

This is the repository where the scripts and input files associated to the paper "XXX" by Cristina Claver, Amaia Astarloa, Iñaki Mendibil, Ignacio Molpeceres-Diego, Antonio Fernandez, Maite Louzao, Guillermo Boyra and Naiara Rodríguez-Ezpeleta are locaed. 
*Link to the preprint:  XXXXX*

## Metabarcoding data processing
Scripts used for raw data (available at SRA under Bioproject XXXXX) preparation and analysis.

*Demultiplexing step:*
- 1.script_demultiplexing.sh
- Input files: M*_barcodes.sh
- Input files: Mix_path_to_R1_R2.txt

*Cleaning step:*
- 2.script_prepare_reads.sh
- Input files: Ceto_samples_20251020.txt

*Taxonomic assignments:*
- 3.script_tax_assignment.sh
- Input files are outputted from previous step

*Visualization of results:*
- Plot_mock_eDNA_samples.R

## Sighting and dPCR data processing

*Spatial gridding and zonation of data:*
- Data_processing_visualization.R (cleans visual and molecular (dPCR) datasets. Maps effort and detections across the study area)
- Standardization_visu_molecular.R (builds a spatial grid, assigns visual and molecular observations to each cell, computes encounter and detection rates and maps those metrics)

## Modelling

*Bayesian models*
- edna_only_2.stan (biomass estimation model based on molecular data)
- sightings_only.stan (biomass estimation model based on sighitngs data)
- joint_eDNA_signtings.stan (biomass estimation model based on both molecular and sightings data combined)

*Visualization*

These scripts prepare sightings and molecular data, fit three models for each species, visualize and compare posterior abundance estimates. 
- modeling_dolphins.R
- modeling_fin_whales.R
- modeling_beaked_whales.R
- modeling_pilot_whales.R

## Contact
Cristina Claver (cclaver@azti.es/cristinac@setur.fo)
