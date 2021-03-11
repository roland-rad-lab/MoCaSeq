library(data.table)
library(splitstackshape)
args = commandArgs(trailingOnly=TRUE)

tool <- args[1] # e.g. theta2
file <- args[2] # e.g. theta2/SAMPLE.n2.results
  
if (length(args)!=2) {
  print("Exactly two argument must be supplied, all others are ignored")
}

if(tool == "theta2"){
  dt <- fread(file, select = "mu")
  dt <- cSplit(dt, "mu", ",", direction = "long")
  dt[1, population := "normal"]
  
  if(nrow(dt[is.na(population)]) == 1){
    dt[is.na(population), population := "tumor"]
  } else {
    dt[is.na(population), population := paste0("tumor_sub", .I)]
  }
  
  dt[, mu := round(mu, digits = 2)]
  
  dt <- dt[, .(population, purity=mu)]
  print(dt)
  
  outfile <- gsub(".results$",".results.table",file)
  fwrite(dt, outfile, sep="\t")
} else {
  stop(paste0("Tool (",tool,") unknown"))
}