# final_masters_code
 Repo to store the final code that is attached to Master's mini-dissertation

 # Folder structure
 All the files for the Master's mini-dissertation was generated using R. You do need to state where you stored the files in the file_path variable.
 
There are three main files where the sections' plots and values used in tables can be generated.
1. one for the univariate simulation `final_simulation_code.R`
2. one for the multivariate simulation `final_multivariate_simulation_code.R`
3. and one for the application section `final_application_code.R`.
4. For the log-likelihood graph in the univariate simulation section, you can generate the graph by running `final_likelihood_graph.R`.

The models' logic is stored in two files (univariate and bivariate)
1. Incremental Model
- Univariate: `final_IEM_model_univariate.R`
- Bivariate: `final_IEM_model_multivariate.R`
2. Parallel Model
- Univariate: `final_parallel_model_univariate.R`
- Bivariate: `final_parallel_multivariate.R`
3. Centralised Model
- Univariate: `final_centralised_model_univariate.R`
- Multivariate: `final_centralised_multivariate.R`

Stored logic that is used for general purpose or by all the models in two helper scripts:
1. Script to store logic that all/most models user: `final_all_models_use.R`
2. Script to store logic for general data generation and data plotting: `final_data_generation.R`

Note: There is another script `Checking_of_model_results_final.R` that does not necessarily generate any output that was used by the other scripts or in the paper but is a nice script to run to check the prediction results of the different multivariate models. 
