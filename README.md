# bmitb
BMI distributions &amp; TB

## Software used in this analysis

This analysis used R veresion 4.5.0 and the following packages:


| Package | Version |
|---------|---------|
| data.table | 1.17.8 |
| ggplot2 | 3.5.2 |
| ggpubr | 0.6.1 |
| ggrepel | 0.9.6 |
| glue | 1.8.0 |
| googlesheets4 | 1.1.1 |
| here | 1.0.1 |
| officer | 0.6.10 |
| paletteer | 1.6.0 |
| readxl | 1.4.5 |
| rvg | 0.3.5 |
| scales | 1.4.0 |
| sf | 1.0-21 |
| sp | 2.2-0 |
| wbmapdata | 0.0.0.9000 |
| MASS | 7.3-65 |

These are available from CRAN excep `wbmapdata` which is available from https://github.com/petedodd/wbmapdata/


## Data used in this analysis

A single archive of the public input data to reproduce this analysis has been posted on Zenodo at:
https://zenodo.org/records/16900137

which includes data from the following sources.

Estimates for 2022 of BMI by country, sex, and age from the NCD-RisC consortium,[^1] which were downloaded from https://ncdrisc.org/data-downloads-adiposity.html on 14/Feb/2024. 
In particular, we used the files:

- `NCD_RisC_Lancet_2024_BMI_child_adolescent_country.csv`
- `NCD_RisC_Lancet_2024_BMI_female_age_specific_country.csv`
- `NCD_RisC_Lancet_2024_BMI_male_age_specific_country.csv`

[^1]: NCD Risk Factor Collaboration (NCD-RisC). Worldwide trends in underweight and obesity from 1990 to 2022: a pooled analysis of 3663 population-representative studies with 222 million children, adolescents, and adults. Lancet 2024; 403: 1027–50.


For the adolescent group aged 15-19 years, we also made use of WHO reference tables to convert NCD-RisC estimates of z-scores into BMIs. These were downloaded from https://www.who.int/tools/growth-reference-data-for-5to19-years/indicators/bmi-for-age on 27/Nov/2023. In particular, we used the files:

- `bmi-boys-z-who-2007-exp.xlsx`
- `bmi-girls-z-who-2007-exp.xlsx`

We used WHO estimates of TB incidence, which were downloaded from https://www.who.int/teams/global-programme-on-tuberculosis-and-lung-health/data on 30/Oct/2024. In particular, we used the files:

- `TB_burden_age_sex_2024-10-30.csv`
- `TB_notifications_2024-10-30.csv`

We used World Population Prospects 2024 demographic estimates from the United Nations Population Division, downloaded from  https://population.un.org/wpp/ on 19/August/2025. These were aggregated into relevant age categories for 2023 and included as the file:

- `N8523.Rdata`

We also used coefficients and variance covariance matrices from the regression in Saunders et al.[^2]:

- `general_population_vcov_matrix.csv`
- `general_population_piecewise_parameters.csv`

[^2]: Saunders MJ, Cegielski JP, Clark R, Houben RMGJ, McQuaid CF. Body mass index and tuberculosis risk – an updated systematic literature review and dose-response meta-analysis. Int J Epidemiol 2025; 54. DOI:10.1093/ije/dyaf154.

