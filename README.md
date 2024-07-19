# Adelie_forecasts

Forecasts of Adélie penguin colony growth rates using space-for-time substitution.

To replicate the results first run "run_adpe_models.R" then "forecast_growth.R".

**scripts/**

&emsp; **run_adpe_models.R**: Fits the relationship between sea ice concentration and average colony growth with quadratic relationship in log scale.
&emsp;
&emsp; **wrangle_data_forecast_adpe.R**: Shapes the climate projection data into sutaible format for population projections with necessary transformations and bias corrections.
&emsp;
&emsp; **forecast_growth.R**: Using the model fit above project Adélie penguin colony growth rates.
&emsp;
&emsp; **get_starting_pop_size.R**: Calculate starting population size to be later used in meta-population models.

**data/**

&emsp; **data_pop_adpe.rds**: General population data later to be shaped into proper format to use in population models to estimate colony growth.
&emsp;
&emsp; **data_stan_null_adpe.rds**: Data shaped into stan format to be used in population models to estimate colony growth 
&emsp;
&emsp; **data_forced_finn_500km_adpe.rds**: Observed sea ice concentration data used in model fitting.
&emsp;
&emsp; **data_aice_coupled_adpe.rds**: Raw climate projection data of sea ice concentration across 50 ensemble members.
&emsp;
&emsp; **data_coupled_transformed_adpe.rds**: Bias corrected and shaped climate data to be used in colony growth projections.
&emsp;
&emsp; **data_coupled_normalized_adpe.rds**: Bias corrected, shaped and standardized (with observed mean and variance) climate data to be used in colony growth projections.
&emsp;
&emsp; **colonies_adelie.csv**: Colony coordinates.
&emsp;
&emsp; **pop_size_starting.csv**: Starting population size for Adélie colonies reported with mean, median, sd and credible intervals.
&emsp;
&emsp; **colony_info_adelie.csv**: Colony coordinates combined with starting pop size information.


