These R scripts compute the model fit parameters for the article "Kinetic modeling of stem cell transcriptome dynamics to identify regulatory modules of normal and disturbed neuroectodermal differentiation". Here, the GPL (>=2) applies.
 
The Snakefile describes the workflow and can be executed using snakemake, downloadable from https://github.com/snakemake/snakemake.

In order to execute the Snakefile, the following additional file is required:

fit_data_in.csv: This file can be created using the code in the repository https://github.com/johannesmg/kinetic_modeling_stem_cell_transcriptome_dynamics.

The configuration of the local directories has to be entered into the config.yml file.