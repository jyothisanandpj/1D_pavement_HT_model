# 1D_pavement_HT_model
The model will estimate the pavement surface temperature and convectional heat transfer between the pavement surface and the surrounding air



								RUN THE MODEL: Double click "EXECUTE.EXE" file in the folder
								
			NOTE: OUTPUT.CSV file should NOT BE OPENED at the time of execution
				
==============================================================================================================================================				
				
				MODIFY INPUTS: Both material and radiative properties of the pavement can be modified
				
	MAT_PROP.CSV: In this file, you can edit: denisty (kg/cu.m), Specific heat capacity (J/kg-K), conductivity(W/m-K), and length (m).
				  The model is developed for 4 layer pavement structure, therefore while executing the simulation there should be properties 
				  for all the 4 layers, as shown in the DEFUALT VALUES. In case, if your model has only 3 layers, please split the 3rd layer 
				  into 2 layers as shown in the DEFUALT VALUES.
				  NOTE: Length of each layer can be modified, however, total length of all the 4 layers (sum) should be 1 m to meet the 
				  ground boundary condition set in the model.				  
	DEFUALT VALUES:
	
denisty (kg/cu.m),Specific heat capacity (J/kg-K),conductivity (W/m-K),length (m)
2238,921,1.21,0.1
2238,921,1.21,0.15
1500,1900,1,0.25
1500,1900,1,0.5	
	
	RAD_PROP.CSV: The radiative properties (albedo and emissivity) of the pavement surface can be modified in the RAD_PROP.CSV file.
				  NOTE: Both values should be between 0 and 1, otherwise leads to incorrect results.		  
	DEFUALT VALUES:
	
albedo,emissivity
0.17,0.85
	
==============================================================================================================================================

								RESULTS OF THE SIMULATION can be seen in the OUTPUT.CSV file 
	
	Simulation ran for five days based on Phoenix, AZ, USA weather data from 6/22/2004 to 6/26/2004. The first 3 days simulation is used for 
	spinning the model and remaining two says results can be used for analaysis. The hourly averaged results (pavement surface 
	tempaerature and sensible heat flux) can be seen in the OUTPUT.CSV file and the weather profile for 6/25/2004 and 6/26/2004 can be seen in 
	the	WEATHERFILE.CSV.


==============================================================================================================================================

								Model Validation and Citation: https://doi.org/10.1016/j.solener.2021.10.056
				
