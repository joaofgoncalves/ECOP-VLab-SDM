FROM jfgoncalves/r-biomod2:latest

MAINTAINER "Joao F Goncalves" "joaofgo@gmail.com"

## Copy biomod2 script
COPY run_biomod2.r /usr/local/src

## Input covariates data 
COPY vars.zip /usr/local/src

## Input PA records for the target species
COPY input_records.csv /usr/local/src

## Parameters for running biomod
COPY parameters.csv /usr/local/src

## Maxent java files
COPY maxent.jar /usr/local/src
COPY maxent.sh /usr/local/src

## Run R script
WORKDIR /usr/local/src
CMD ["Rscript", "--vanilla run_biomod2.r vars.zip input_records.csv parameters.csv"]