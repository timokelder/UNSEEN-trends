# UNSEEN-Trends
In this repository, we provide code for the paper UNSEEN-trends: Towards detection of changes in 100-year precipitation events over the last 35 years.


## Short intro
In this study, we use Norway and Svalbard as a case study to explore the UNSEEN-trends approach: using [SEAS5, Johnson et al., 2019](https://www.geosci-model-dev.net/12/1087/2019/) to analyze trends in extreme precipitation. These regions have recently faced severe events, raising the question whether these kind of events have occurred by chance or a the new norm. From observations, it is impossible to analyze the changes in these severe events: How have the 100-year precipitation events changed over last 35 years? By pooling ensemble members and lead times into an UNprecedented Simulated Extreme ENsemble (UNSEEN, [Thompson et al., 2017](https://www.nature.com/articles/s41467-017-00275-3)), we can try to answer this question. 

## Data
SEAS5 hindcasts are provided every month between the years 1981-2016, with 25 ensemble members. The hindcasts are downloaded from the ECMWF MARS server onto the server of the Norwegian Meteorological Institute (MetNor). We then mined the hindcasts, i.e. we extracted the target variable (the maximum three-day precipitation event in autumn: SON-3DP) from each forecast to generate the UNSEEN ensemble. This process is performed on the cluster of MetNor [generating the UNSEEN ensemble](Mining/Mine_regional_average.sh) but unfortunately cannot be reproduced without access to the server of MetNor. For verification over Norway, we use the daily gridded observed record [SeNorge](https://www.earth-syst-sci-data.net/11/1531/2019/) and upscale this to the same resolution as SEAS5 and extract the SON-3DP events. We provide the resulting UNSEEN ensembles for Norway and Svalbard, along with the processed observed record for Norway, in the [data folder](Data/Extremes). 

## Evaluation of the UNSEEN ensemble

The framework to evaluate the UNSEEN ensemble presented here consists of testing the [ensemble member independence](Notebooks/Independence.md), [model stability](Notebooks/Model_stability.md) and [model fidelity](Notebooks/Model_fidelity.md).

The independence of ensemble members is an important requirement for the UNSEEN approach, as dependent members would merely duplicate the data and artificially inflate the sample size, without adding new information. A second issue when using seasonal forecasts for pooled evaluation of extremes is the presence of inconsistencies in climate models, which can introduce drift in the simulated climatology and alter precipitation extremes over longer lead times. Therefore, model stability is a requirement for pooling lead times. Trusting the simulated `unprecedented extremes' in large ensembles is further complicated by the inability to validate extremes, given the limited sample sizes of observations. We follow previous UNSEEN studies ( [Thompson et al., 2017](https://www.nature.com/articles/s41467-017-00275-3)) in evaluating the modelled extremes by bootstrapping the ensemble into datasets of 35 years and assessing whether the observations fell within the range of the bootstrapped distribution. We perform this analysis for the SEAS5-100 SON-3DP ensemble over Norway. Norway's dense station network enables large-scale precipitation comparisons  (in contrast with Svalbard). For a comprehensive global model validation of SEAS5, see [Johnson et al., 2019](https://www.geosci-model-dev.net/12/1087/2019/).

## The UNSEEN-trends approach

How have the 100-year precipitation events changed over last 35 years? The large UNSEEN sample can be used to answer this question! We illustrate the [UNSEEN-Trends approach](Notebooks/Trend_analysis.md) for Norway and Svalbard and show that the 100-year event in 1981 over Svalbard can now be expecter with a return period of 41 years. For Norway, we find a small increase of 2%, with confidence intervals from -3%/7%. 

## Outlook

We envision that further applications can 1) help estimateing design values, especially relevant for data scarce regions; 2) improve risk estimation of natural hazards by coupling UNSEEN to impact models; 3) detect trends in rare climate extremes, including variables other than precipitation; and 4) increase our physical understanding of the drivers of non-stationarity of climate extremes, through the possible attribution of detected trends. 

We hope this approach may see many applications across a range of scientific fields!

Please don't hesitate to ask any questions. 
