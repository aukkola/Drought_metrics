Python codes for calculating common drought metrics from gridded data

SPI-based metrics are calculated as per:
Ukkola et al., Evaluating CMIP5 model agreement for multiple drought metrics,
Journal of Hydrometeorology, 2018
https://journals.ametsoc.org/doi/abs/10.1175/JHM-D-17-0099.1

see: Example for_calculating_drought_metrics_SPI_Precip_CMIP5.py


Percentile-based metrics use a set percentile for determining droughts. Percentile codes updates to codes used in:
Ukkola et al., Robust future changes in meteorological drought in CMIP6 projections despite uncertainty in precipitation, Geophysical Research Letters, 46, e2020GL087820. 
https://doi.org/10.1029/2020GL087820

see: Example_for_calculating_historical_drought_metrics_from_CMIP6_using_percentiles.py


For more examples on the use of percentile codes, including calculating future drought indices using a historical baseline, see:
https://bitbucket.org/aukkola/cmip6_drought_projections/src/master/Python/


The above percentile codes use the same threshold (default: 15th percentile) for drought onset and termination. This file adapts these codes to use different onset and termination thresholds (10th and 50th percentile by default, respectively). With this method, drought commences when the time series goes below the onset threshold and ends when the termination threshold is next exceeded:

Example_for_calculating_historical_drought_metrics_from_CMIP6_using_percentiles_two_threshold.py

