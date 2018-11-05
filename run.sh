# Create the main output directory
mkdir -p ./output

# Run R script file
Rscript --vanilla run_biomod2.r "vars.zip" "input_records.csv" "parameters.csv"
