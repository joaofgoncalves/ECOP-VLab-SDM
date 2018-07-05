
### ------------------------------------------------------------------------------ ###
### 
### ECOPOTENTIAL VirtualLab r-biomod2 model porting
### ICETA / CIBIO / InBIO team
###
### Salvador Arenas-Castro, Joao Goncalves
### Porto, 2018
###
### This script contains a generalized version of the SDM correlative model developed 
### using R/biomod2 and aiming to predict the current distribution of Iris boissieri 
### based on satellite-derived Ecosystem Functional Attributes presented in:
###
### Arenas-Castro, S. et (2018). Assessing the multi-scale predictive ability of 
### ecosystem functional attributes for species distribution modelling. PLoS One.
### 
### http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0199292
### https://doi.org/10.1371/journal.pone.0199292
###
### ------------------------------------------------------------------------------ ###


library(raster)
library(biomod2)

# Set working directory
setwd("/usr/local/src")

# Get arguments from the command line
# [1] a zip file containing the raster covariates
# [2] a csv file containing two columns [X, Y] with presence record coordinates
# [3] a csv files with parameters

# Get command line values
args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied!", call.=FALSE)
}


### ------------------------------------------------------------------------------ ###
### 1 | LOAD PARAMETERS AND DATA FROM FILES                                        ###
### ------------------------------------------------------------------------------ ###

# Read parameter csv file
pars <- read.csv(args[3], stringsAsFactors = FALSE, comment.char = "#")

# if(nrow(pars) != 3)
#   stop("The csv parameter file must have four rows!")
if(ncol(pars) != 2)
  stop("The csv parameter file must have two columns!")

rownames(pars) <- pars[,1]
print(pars)

resp.name <- as.character(pars["resp.name",2])
PA.nb.absences <- as.integer(pars["PA.nb.absences",2])
PA.nb.rep  <- as.integer(pars["PA.nb.rep",2])
models <- eval(parse(text=as.character(pars["models",2])))
NbRunEval <- as.integer(pars["NbRunEval",2])
DataSplit <- as.integer(pars["DataSplit",2])
Yweights <- as.character(pars["Yweights",2])
Prevalence <- as.numeric(pars["Prevalence",2])
VarImport <- as.integer(pars["VarImport",2])
models.eval.meth <- eval(parse(text=as.character(pars["models.eval.meth",2])))
SaveObj <- as.logical(as.character(pars["SaveObj",2]))
rescal.all.models <- as.logical(as.character(pars["rescal.all.models",2]))
do.full.models <- as.logical(as.character(pars["do.full.models",2]))

print(resp.name)
print(PA.nb.rep)
print(PA.nb.absences)
print(models.eval.meth)

# Unzip the covariates data
unzip(zipfile = args[1], overwrite = TRUE)

# List all GeoTIFF files to be used for 
fl <- list.files(pattern=".tif$", full.names = TRUE)

# Environmental predictors
myExpl <- stack(fl)
print(myExpl)

# Load species data
DataSpecies <- read.csv(args[2], sep=',')
# Check species data
head(DataSpecies)

# Create spatial points object / Presence/absence data for the target species
# Set CRS as the same for myExpl raster stack
myResp <- SpatialPoints(DataSpecies[,c("X","Y")], proj4string = crs(myExpl))
print(myResp)


### ------------------------------------------------------------------------------ ###
### 2 | BIOMOD DATA FORMATING/PREPARATION                                          ###
### ------------------------------------------------------------------------------ ###

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.name = resp.name,
                                     eval.resp.var = NULL,
                                     eval.expl.var = NULL,
                                     eval.resp.xy = NULL,
                                     PA.nb.rep = PA.nb.rep,
                                     PA.nb.absences = PA.nb.absences,
                                     PA.strategy = 'random',
                                     PA.dist.min = 0,
                                     PA.dist.max = NULL,
                                     PA.sre.quant = 0.025,
                                     PA.table = NULL,
                                     na.rm = TRUE)

# Check biomod formatted data
print(myBiomodData)


### ------------------------------------------------------------------------------ ###
### 3 | DEFINE MODEL OPTIONS                                                       ###
### ------------------------------------------------------------------------------ ###

opts <- BIOMOD_ModelingOptions(
  GLM = eval(parse(text=as.character(pars["opts.GLM",2]))),
  GBM = eval(parse(text=as.character(pars["opts.GBM",2]))),
  GAM = eval(parse(text=as.character(pars["opts.GAM",2]))),
  CTA = eval(parse(text=as.character(pars["opts.CTA",2]))),
  ANN = eval(parse(text=as.character(pars["opts.ANN",2]))),
  SRE = eval(parse(text=as.character(pars["opts.SRE",2]))),
  FDA = eval(parse(text=as.character(pars["opts.FDA",2]))),
  MARS = eval(parse(text=as.character(pars["opts.MARS",2]))),
  RF = eval(parse(text=as.character(pars["opts.RF",2]))),
  MAXENT.Phillips = eval(parse(text=as.character(pars["opts.MAXENT.Phillips",2]))),
  MAXENT.Tsuruoka = eval(parse(text=as.character(pars["opts.MAXENT.Tsuruoka",2])))
)

print(opts)

myBiomodOption <- BIOMOD_ModelingOptions(GAM = list(k = 4))


### ------------------------------------------------------------------------------ ###
### 4 | MODEL TRAINING/CALIBRATION                                                 ###
### ------------------------------------------------------------------------------ ###

myBiomodModelOut <- BIOMOD_Modeling( 
  myBiomodData, # data including species records, species name, and variables
  models = c('GLM','GAM'), # models to run
  models.options = myBiomodOption, # options for modelling
  NbRunEval = 3, # 30,# number of evaluations runs
  DataSplit = 80, # % of data used as training
  Yweights = NULL,
  Prevalence = 0.5,
  VarImport = 10, #number of permutation to estimate variable importance
  models.eval.meth = c('TSS','ROC','KAPPA'), #evaluation metrics
  SaveObj = TRUE,            # keep all results and outputs on hard drive
  rescal.all.models = TRUE,  # all models scaled with a binomial GLM
  do.full.models = TRUE)     # models calibrated and evaluated with whole dataset
  # modeling.id = "allmodels") # ID (=NAME) of modelling procedure

# MODELING SUMMARY
print(myBiomodModelOut)
print(list.files(paste("./",resp.name,sep="")))


### ------------------------------------------------------------------------------ ###
### 5 | MODEL EVALUATION STATS                                                     ###
### ------------------------------------------------------------------------------ ###

# GET ALL MODEL EVALUATIONs
modEvalStats <- as.data.frame(get_evaluations(myBiomodModelOut))
write.csv(modEvalStats, paste("./",resp.name,"/modEvalStats.csv",sep=""))

# PRINT VARIABLE IMPORTANCES                                    
varImpObj <- as.data.frame(get_variables_importance(myBiomodModelOut))
write.csv(varImpObj, paste("./",resp.name,"/variableImportances.csv",sep=""))


### ------------------------------------------------------------------------------ ###
### 6 | ENSEMBLE MODELING                                                          ###
### ------------------------------------------------------------------------------ ###

myBiomodEM <- BIOMOD_EnsembleModeling( 
  modeling.output = myBiomodModelOut, # model results
  chosen.models = 'all', # models to be included when ensembling
  em.by='all', # flag defining the way the models will be combined to build the ensemble models: 'PA_dataset+repet' (default), 'PA_dataset+algo', 'PA_dataset', 'algo', 'all'
  eval.metric = c('ROC'), # evaluation metric used to build ensemble models
  eval.metric.quality.threshold = c(0.7), # If not NULL, the minimum scores below which models will be excluded of the ensemble-models building
  prob.mean = TRUE, # estimate the mean probabilities across predictions
  prob.cv = FALSE, # estimate the coefficient of variation across predictions
  prob.ci = FALSE, # estimate the confidence interval around the prob.mean
  prob.ci.alpha = 0.05, # significance level for estimating the confidence interval. Default = 0.05
  prob.median = FALSE, # estimate the mediane of probabilities
  committee.averaging = FALSE, # estimate the committee averaging across predictions
  prob.mean.weight = FALSE, # estimate the weighted sum of probabilities
  prob.mean.weight.decay = 'proportional') # define the relative importance of the weights

#PRINT SUMMARY                     
print(myBiomodEM)


### ---------------------------------------------------------------------------------- ###
### 7 | PROJECTION OF MODELS TO THE FULL GEO/ENV SPACE USING THE TRAINING RASTER DATA  ###
### ---------------------------------------------------------------------------------- ###

#PROJECTION OVER THE STUDY AREA UNDER CURRENT CONDITIONS
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut, # modelling results
  new.env = myExpl, # environmental variables
  proj.name = 'current', # name of projections
  selected.models = get_kept_models(myBiomodEM, model = 1), ## Changed this to include only selected models from BIOMOD_EnsembleModeling criteria
  binary.meth = 'TSS', # a vector of a subset of models evaluation method computed before
  build.clamping.mask = TRUE, # if TRUE, a clamping mask will be saved on hard drive different
  output.format = '.grd', # the format of the GIS files
  do.stack = TRUE)

mod_projPres <- get_predictions(myBiomodProj)
writeRaster(mod_projPres, filename = paste("./",resp.name,"/modelProjectionsByModel.tif",sep=""), overwrite = TRUE)


### ---------------------------------------------------------------------------------- ###
### 8 | ENSEMBLE FORECASTING                                                           ###
### ---------------------------------------------------------------------------------- ###

myBiomodEF <- BIOMOD_EnsembleForecasting( 
  EM.output = myBiomodEM,
  projection.output = myBiomodProj,
  binary.meth = c('ROC', 'TSS'),
  filtered.meth = c('ROC', 'TSS'),
  compress = TRUE,
  output.format = '.grd')

mod_projBiomodEF <- get_predictions(myBiomodEF)
writeRaster(mod_projBiomodEF, filename = paste("./",resp.name,"/modelProjEnsembleForecasting.tif",sep=""), overwrite = TRUE)

## Make a zip archive with all outputs created in the resp.name folder
## This file should be exported from VLab
zip("output.zip", paste("./", resp.name, sep=""))


