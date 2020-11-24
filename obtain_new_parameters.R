library(dplyr)
library(readr)

args = commandArgs(trailingOnly=TRUE)
in.dir=args[1]
out.file=args[2]


file.list=list.files(path=in.dir,full.names=T)

options(readr.num_columns = 0)
new.pars=bind_rows(lapply(file.list,read_csv,col_names=T))

new.pars %>%
    group_by(gene) %>%
    top_n(1,-cost) ->
    out.pars

write_csv(out.pars,path=out.file)
