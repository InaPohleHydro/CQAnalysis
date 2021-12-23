source("Functions_cQ_Analysis.R") #loads the R functions to carry out the c-Q analysis
example_data<- read.table("example_data.txt", header = TRUE) #reads the example data
cQ_example <- analyse_type_cQ_and_chemostatic_chemodynamic(cQ_input_data = example_data) #carries out complete c-Q analysis of example data. if interim results are of interest, please carry out the other functions independently


