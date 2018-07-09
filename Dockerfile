FROM jfgoncalves/r-biomod2:latest

MAINTAINER "Joao F Goncalves" "joaofgo@gmail.com"

## Copy biomod2 script
COPY ./scripts/run_biomod2.r /
COPY ./scripts/run_biomod2.sh /

RUN mkdir -p ./input
RUN mkdir -p ./output

## Input covariates data 
COPY ./input/vars.zip ./input

## Input PA records for the target species
COPY ./input/input_records.csv ./input

## Parameters for running biomod
COPY ./input/parameters.csv ./input

## Maxent java files
COPY ./scripts/maxent.jar ./output

## Run R script -- don't use this
WORKDIR ./output
#CMD ["Rscript", "--vanilla ../run_biomod2.r ../input/vars.zip ../input/input_records.csv ../input/parameters.csv"]