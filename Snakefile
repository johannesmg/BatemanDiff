import pandas as pd

configfile: 'config.yml'

EXPRESSIONDATA=config['expressiondata']
OUTUNTREATED=config['outuntreated']
OUTTREATED=config['outtreated']
PARAMETERFILE=config['parameterfile']
GENSATIME=config['gensatime']
NOCORESTREATEDFIT=config['nocorestreatedfit']
NOREPSTREATEDFIT=config['norepstreatedfit']

genelist = pd.read_csv(EXPRESSIONDATA).set_index("gene", drop=False)
genes=pd.unique(genelist.gene)

rule all:
  input: expand("{outfinal}/fit_out_{gene}.csv", gene=genes,outfinal=OUTTREATED)


rule setup_chunks:
  input: EXPRESSIONDATA
  output: OUTUNTREATED+'/chunks_untreated.csv'
  message: "--- Generating untreated chunks."
  shell: 'Rscript generate_untreated_chunks.R 1 {input} {output} {PARAMETERFILE}'


rule fit_gene:
  input: [OUTUNTREATED+'/chunks_untreated.csv',EXPRESSIONDATA]
  output: OUTUNTREATED+'/fit_out_{gene}.csv'
  message: "--- Fitting gene untreated."
  shell: 'Rscript fit_untreated.R {input[0]} {wildcards.gene} {GENSATIME} 1 {input[1]} {output}'


rule fit_gene_treated:
  input: [EXPRESSIONDATA,OUTUNTREATED+'/fit_out_{gene}.csv']
  output: OUTTREATED+'/fit_out_{gene}.csv'
  message: "--- Fitting gene treated."
  shell: 'Rscript fit_treated.R {input[0]} {input[1]} {PARAMETERFILE} {output} {NOCORESTREATEDFIT} {NOREPSTREATEDFIT}'


