# Adelie_forecasts

Forecasts of Adélie penguin colony growth rates using space-for-time substitution.

To replicate the results first run "run_adpe_models.R" then "forecast_growth.R".

**scripts/**

  &ensp; **run_adpe_models.R**: Fits the relationship between sea ice concentration and average colony growth with quadratic relationship in log scale.

  **wrangle_data_forecast_adpe.R**: Shapes the climate projection data into sutaible format for population projections with necessary transformations and bias corrections.

  **forecast_growth.R**: Using the model fit above project Adélie penguin colony growth rates.

  **get_starting_pop_size.R**: Calculate starting population size to be later used in meta-population models.

**data/**

  **data_pop_adpe.rds**: General population data later to be shaped into proper format to use in population models to estimate colony growth.

  **data_stan_null_adpe.rds**: Data shaped into stan format to be used in population models to estimate colony growth 

  **data_forced_finn_500km_adpe.rds**: Observed sea ice concentration data used in model fitting.

  **data_aice_coupled_adpe.rds**: Raw climate projection data of sea ice concentration across 50 ensemble members.

  **data_coupled_transformed_adpe.rds**: Bias corrected and shaped climate data to be used in colony growth projections.

  **data_coupled_normalized_adpe.rds**: Bias corrected, shaped and standardized (with observed mean and variance) climate data to be used in colony growth projections.
  
  **colonies_adelie.csv**: Colony coordinates.

  **pop_size_starting.csv**: Starting population size for Adélie colonies reported with mean, median, sd and credible intervals.

  **colony_info_adelie.csv**: Colony coordinates combined with starting pop size information.


