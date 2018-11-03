# Create the main output directory
mkdir ./output

# Run R script file
Rscript --vanilla run_biomod2.r "./input/vars.zip" "./input/input_records.csv" "./input/parameters.csv"