This project contains the scripts used to download and process SubX's project data of outgoing longwave radiation (Pegion, K., and Coauthors, 2019: The Subseasonal Experiment (SubX): A Multimodel Subseasonal Prediction Experiment.Bull. Amer. Meteor. Soc., 100, 2043–2060).
(https://journals.ametsoc.org/bams/article/100/10/2043/344809/The-Subseasonal-Experiment-SubX-A-Multimodel)

More on the SubX Project: http://cola.gmu.edu/subx/

Preliminar results are shown in Alvarez M. & Vera C. 2020: Evaluation of a SubX Multimodel Ensemble to predict an intraseasonal index for South America based on Outgoing Long Wave Radiation, presented at the 45th NOAA Climate Diagnostics and Prediction Workshop 
(https://www.researchgate.net/publication/344787078_Evaluation_of_a_SubX_Multimodel_Ensemble_to_predict_an_intraseasonal_index_for_South_America_based_on_Outgoing_Long_Wave_Radiation?channel=doi&linkId=5f90560ba6fdccfd7b72446f&showFulltext=true)

In order to work with the SubX data:

1. Download SubX model data from IRIDL using the "getSUBXIRIXYLMSP.R" (or subdomain version -highly encouraged-) -> netcdf files
2. Download SubX model climatology from IRIDL using "getSUBXIRIXYLMSP_clim_subdomain.R" -> netcdf files
3. Compute anomalies using "compute_rlut_anomalies.R" -> netcdf files
4. Compute ensemble mean using "compute_rlut_ensmean.R" and "compute_rlut_ensmean_NRL.R" -> netcdf files
5. Open ensemble mean from netcdf files and save to single .rdf file to work with R with "opens_ensemblemean_MODEL.R" -> rdf file
6. Compute MME following SubX suggested methodology using "compute_rlut_MME.R" -> .rdf file
