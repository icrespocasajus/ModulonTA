## code to prepare `DATASET` dataset goes here
library(usethis)

network.LCMV = readRDS(file="./network.LCMV.Rds")
network.LCMV = unique(network.LCMV)
modulons.LCMV=readRDS(file="./TF.AUC.clusters.LCMV.subset.Rds")
usethis::use_data(network.LCMV, overwrite = TRUE)
usethis::use_data(modulons.LCMV, overwrite = TRUE)

network.TILs = readRDS(file="./network.TILs.Rds")
network.TILs = unique(network.TILs)
modulons.TILs=readRDS(file="./TF.AUC.clusters.TILs.Rds")
usethis::use_data(network.TILs, overwrite = TRUE)
usethis::use_data(modulons.TILs, overwrite = TRUE)

cc.TILs=readRDS(file="./cc.TILs.Rds")
usethis::use_data(cc.TILs, overwrite = TRUE)
Modulon.Cores.TILs=readRDS(file="./Modulon.Cores.TILs.Rds")
usethis::use_data(Modulon.Cores.TILs, overwrite = TRUE)

cc.LCMV=readRDS(file="./cc.LCMV.Rds")
usethis::use_data(cc.LCMV, overwrite = TRUE)
Modulon.Cores.LCMV=readRDS(file="./Modulon.Cores.LCMV.Rds")
usethis::use_data(Modulon.Cores.LCMV, overwrite = TRUE)

discriminancy.LCMV = readRDS(file = "./REGULON&RNA.LCMV.metaAUC.all_subpopulations_Vs._background.Results.OPLS.Rds")
DA.LCMV = lapply(discriminancy.LCMV,function(x){tv=grep('_REGULON',rownames(x));tmp.df=as.data.frame(x[tv,'Weights']);rownames(tmp.df)=gsub('_REGULON','',rownames(x)[tv]);colnames(tmp.df)='Weights';return(tmp.df)})
usethis::use_data(DA.LCMV, overwrite = TRUE)

discriminancy.TILs = readRDS(file = "./REGULON&RNA.TILs.metaAUC.all_subpopulations_Vs._background.Results.OPLS.Rds")
DA.TILs = lapply(discriminancy.TILs,function(x){tv=grep('_REGULON',rownames(x));tmp.df=as.data.frame(x[tv,'Weights']);rownames(tmp.df)=gsub('_REGULON','',rownames(x)[tv]);colnames(tmp.df)='Weights';return(tmp.df)})
usethis::use_data(DA.TILs, overwrite = TRUE)

RegAUC.TILs=readRDS(file="./RegAUC.TILs.Rds")
usethis::use_data(RegAUC.TILs, overwrite = TRUE)

RegAUC.LCMV=readRDS(file="./RegAUC.LCMV.Rds")
usethis::use_data(RegAUC.LCMV, overwrite = TRUE)
