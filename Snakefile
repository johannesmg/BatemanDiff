import pandas as pd
import random

configfile: 'config.yml'

EXPRESSIONDATA=config['expressiondata']
OUTUNTREATED=config['outuntreated']
OUTTREATED=config['outtreated']
PARAMETERFILE=config['parameterfile']
NOCORESTREATEDFIT=config['nocorestreatedfit']

NOREPSTREATEDFITHIGH=config['norepstreatedfithigh']
NOREPSTREATEDFITLOW=config['norepstreatedfitlow']
NOREPSFINAL=config['norepsfinal']

OPTIMMETHOD=config['optimmethod']
GENSATIMEHIGH=config['gensatimehigh']
GENSATIMELOW=config['gensatimelow']
FNSCALEHIGH=config['fnscalehigh']
FACTRHIGH=config['factrhigh']
MAXITHIGH=config['maxithigh']
FNSCALELOW=config['fnscalelow']
FACTRLOW=config['factrlow']
MAXITLOW=config['maxitlow']



genelist = pd.read_csv(EXPRESSIONDATA).set_index("gene", drop=False)
genes=pd.unique(genelist.gene)


rule all:
  input: OUTTREATED+'/parameters_third_round.csv'


rule setup_chunks:
  input: EXPRESSIONDATA
  output: OUTUNTREATED+'/chunks_untreated.csv'
  resources: 
    mem_per_cpu=4000
  message: "--- Generating untreated chunks."
  shell: 'Rscript generate_untreated_chunks.R 1 {input} {output} {PARAMETERFILE}'


rule fit_gene_untreated_first:
  input: [OUTUNTREATED+'/chunks_untreated.csv',EXPRESSIONDATA]
  output: OUTUNTREATED+'/first/fit_out_{gene}.csv'
  resources: 
    mem_per_cpu=4000
  message: "--- Fitting gene untreated first round."
  shell: 'Rscript fit_untreated.R {input[0]} {wildcards.gene} {GENSATIMEHIGH} 1 {input[1]} {output}'


rule fit_gene_treated_first:
  input: [EXPRESSIONDATA,OUTUNTREATED+'/first/fit_out_{gene}.csv']
  output: OUTTREATED+'/first/fit_out_{gene}.csv'
  threads: NOCORESTREATEDFIT
  resources: 
    mem_per_cpu=2000
  message: "--- Fitting gene treated first round."
  shell: 'Rscript fit_treated.R {input[0]} {input[1]} {PARAMETERFILE} {output} {threads} {NOREPSTREATEDFITLOW} {FNSCALELOW} {FACTRLOW} {MAXITLOW}'

rule obtain_first_round_parameters:
  input: expand("{outtreated}/first/fit_out_{gene}.csv",outtreated=OUTTREATED,gene=genes)
  output: OUTTREATED+'/first/parameters_first_round.csv'
  threads: 1
  resources: 
    mem_per_cpu=2000
  message: "--- Obtaining parameters first round."
  shell: 'Rscript obtain_new_parameters.R {OUTTREATED}/first {output}'

rule scale_errors:
  input: [OUTTREATED+'/first/parameters_first_round.csv',EXPRESSIONDATA]
  output: OUTTREATED+'/data/fit_data_scaled_error.csv'
  threads: 1
  resources: 
    mem_per_cpu=2000
  message: "--- Scaling data errors."
  shell: 'Rscript determine_scaling_factor.R {input[0]} {input[1]} {output}'

         
rule fit_gene_untreated_second:
  input: [OUTUNTREATED+'/chunks_untreated.csv',OUTTREATED+'/data/fit_data_scaled_error.csv']
  output: OUTUNTREATED+'/second/fit_out_{gene}.csv'
  resources: 
    mem_per_cpu=4000
  message: "--- Fitting gene untreated second round."
  shell: 'Rscript fit_untreated.R {input[0]} {wildcards.gene} {GENSATIMEHIGH} 1 {input[1]} {output}'


rule fit_gene_treated_second:
  input: [OUTTREATED+'/data/fit_data_scaled_error.csv',OUTUNTREATED+'/second/fit_out_{gene}.csv']
  output: OUTTREATED+'/second/fit_out_{gene}_{fnscalehigh}_{factrhigh}.csv'
  threads: NOCORESTREATEDFIT
  resources: 
    mem_per_cpu=2000
  message: "--- Fitting gene treated second round."
  shell: 'Rscript fit_treated.R {input[0]} {input[1]} {PARAMETERFILE} {output} {threads} {NOREPSTREATEDFITHIGH} {wildcards.fnscalehigh} {wildcards.factrhigh} {MAXITHIGH}'
         
rule fit_gene_treated_third:
  input: [OUTTREATED+'/data/fit_data_scaled_error.csv',OUTTREATED+'/second/fit_out_{gene}_{fnscalehigh}_{factrhigh}.csv']
  output: OUTTREATED+'/third/final_2_out_{gene}_{method}_{fnscalehigh}_{factrhigh}.csv'
  threads: NOCORESTREATEDFIT
  resources: 
    mem_per_cpu=2000
  message: "--- Fitting gene treated third round."
  shell: 'Rscript final_optimization_different_optimizers.R {input[0]} {input[1]} {PARAMETERFILE} {output} {threads} {NOREPSFINAL} {wildcards.method} {wildcards.fnscalehigh} {wildcards.factrhigh} {MAXITHIGH} {GENSATIMELOW}'

rule obtain_third_round_parameters:
  input: expand("{outtreated}/third/final_2_out_{gene}_{method}_{fnscalehigh}_{factrhigh}.csv", gene=genes,outtreated=OUTTREATED,method=OPTIMMETHOD,factrhigh=FACTRHIGH,fnscalehigh=FNSCALEHIGH,norepstreatedfithigh=NOREPSTREATEDFITHIGH)
  output: OUTTREATED+'/parameters_third_round.csv'
  threads: 1
  resources: 
    mem_per_cpu=2000
  message: "--- Obtaining parameters first round."
  shell: 'Rscript obtain_new_parameters.R {OUTTREATED}/third {output}'
