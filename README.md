# UNSEEN-Trends
In this repo, we provide code for the paper UNSEEN-trends: Towards detection of changes in 100-year precipitation events over the last 35 years.


## Short intro
In this study, we use Norway and Svalbard as a case study to explore the UNSEEN-trends approach: using [SEAS5, Johnson et al., 2019](https://www.geosci-model-dev.net/12/1087/2019/gmd-12-1087-2019.pdf) to analyze trends in extreme precipitation. These regions have recently faced severe events, raising the question whether these kind of events have occurred by chance or a the new norm. From observations, it is impossible to analyze the changes in these severe events: How have the 100-year precipitation events changed over last 35 years? By pooling ensemble members and lead times into an UNprecedented Simulated Extreme ENsemble (UNSEEN, [Thompson et al., 2017](https://www.nature.com/articles/s41467-017-00275-3)), we can try to answer this question. 

In notebooks, We discuss the requirements of [ensemble member independence](Notebooks/Predictive_skill.md)) and [model stability](Notebooks/Predictive_skill.md) and we perform a [fidelity analysis](Notebooks/Predictive_skill.md)). We then illustrate the [UNSEEN-Trends approach](Notebooks/Predictive_skill.md). We hope this approach may see many applications across a range of scientific fields!
