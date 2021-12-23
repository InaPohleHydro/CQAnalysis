#citation:  Pohle, I., Baggaley, N., Palarea-Albaladejo, J., Stutter, M., & Glendell, M. (2021). A framework for assessing concentration-discharge catchment behavior from low-frequency water quality data. Water Resources Research, 57, e2021WR029692. https://doi.org/10.1029/2021WR029692 
source("Functions_cQ_Analysis.R") #loads the R functions to carry out the c-Q analysis
example_data<- read.table("example_data.txt", header = TRUE) #reads the example data
cQ_example <- analyse_type_cQ_and_chemostatic_chemodynamic(cQ_input_data = example_data) #carries out complete c-Q analysis of example data. if interim results are of interest, please carry out the other functions independently


