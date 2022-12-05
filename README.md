# KRE (kernel regression estimation)
These are the data and code mentioned in the manuscript (Family-specific training improves linear B cell epitope prediction for emerging viruses).

`KRE` utilizes sequence homology to predict the epitopes of a novel virus. Users can define their own training dataset, while we recommend a training dataset only including the data from the other viruses in the same family. In the manuscript, we compared the family-specific and general predictions in detail. Family-specific prediction is a better choice.

The file description is as follows:  
+---code  
|　　　　cross_virus_prediction_example.R # an example for seropositive rate data (four viruses)  
|　　　　family_dist_show.R # plot the virus distance map.  
|　　　　family_pred.R # family and general prediction on IEDB data.  
|　　　　find_threshold.R # try to find the optimal evolution threshold for each virus.  
|　　　　kre.R # the main function  
|　　　　kre_parallel.R # the main function in a parallel version  
|　　　　phylosim_generate.R # generate the simulated dataset  
|　　　　plot_family_pred.R # plot the bar plot for predictions on IEDB
|  
+---data # seropositive rate data.  
|　　　　DENV_B_cell_epitope.csv   
|　　　　SARS-CoV+-+B+cell+epitope.csv    
|　　　　SARS-CoV-2+-+B+cell+epitope.csv   
|　　　　ZIKV_B_cell_epitope.csv   
|  
+---patients   
|　　+---code  
|　　|　　　ana.R # calculate the inconsistency of seropositive rate data.  
|　　|　　　dist_show.R # plot two distributions to imitate the positive and negative distributions.   
|　　|       
|　　\---data  
|　　　　　　SARS-CoV-2+Health.csv # healthy people data for SARS-CoV-2.  
|　　　　　　SARS-CoV-2+Patient.csv # patient data for SARS-CoV-2.  
|  
+---Virus # epitope data from IEDB.  

I hope you can find something useful here. :)
