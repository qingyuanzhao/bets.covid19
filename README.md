# Data for the 2019-nCoV outbreak

On 30th of January, 2020, The World Health Organization has declared the [2019-nCoV outbreak a Public Health Emergency of International Concern](https://www.who.int/news-room/detail/30-01-2020-statement-on-the-second-meeting-of-the-international-health-regulations-(2005)-emergency-committee-regarding-the-outbreak-of-novel-coronavirus-(2019-ncov)). **The outbreak originated from my hometown, Wuhan.**

This dataset is collected from public agencies or news media, containing information about the 2019-nCov cases confirmed in and outside China. This dataset is free to use and share given that appropriate credits are given (see the [license](./LICENSE.md)). It can be loaded in R as a package:
```r
devtools::install_github("qingyuanzhao/2019-nCov-Data")
library(nCoV2019.data)
```

A preliminary report is avaiable as a [on medRxiv](https://www.medrxiv.org/content/10.1101/2020.02.06.20020941v1). This report is currently being peer-reviewed. It is based on the analysis in [Feb6.R](https://github.com/qingyuanzhao/2019-nCov-Data/blob/master/1st-Report/Feb6.R) and [Feb6.Rmd](https://github.com/qingyuanzhao/2019-nCov-Data/blob/master/1st-Report/Feb6.Rmd) (results available [here](https://htmlpreview.github.io/?https://github.com/qingyuanzhao/2019-nCov-Data/blob/master/1st-Report/Feb6.html)).

## To Contribute

I am hoping this can become a collaborative and transparent project by people across the world. You can contribute by either updating the dataset or by forming a team to analyze the dataset.

### Dataset

There are two ways to contribute to building this dataset (please only use publicly available information that I can confirm):

1. You can suggest comments in [this Google Spreadsheet](https://docs.google.com/spreadsheets/d/1H4MzVxkug2txyzkiDJsGVKB04YveYcsHg9ijuer8clE/edit?usp=sharing). The easiest way to help is to pick a random row and verify the information is correct by reading the link source. Record your result by commenting on the "Verified" cell in that row. Currently I am having a hard time to obtain detailed information for cases in *Australia*, *France*, *Thailand*, *United Kingdom*, *United States*, and *Vietnam*.

2. You can also use the [Issues](https://github.com/qingyuanzhao/2019-nCov-Data/issues) to record information for new cases. Make sure you read the lessons below before posting.

**Lessons I learned when building this dataset:**
1. News articles don't always report the cases in the same order. It's useful to record the nationality/residence, gender and age of the cases to distinguish them.
2. The most useful columns for data analysis are
  - *Outside* (if the case is infected outside Wuhan). "Y" means yes, "L" means likely, empty means (almost certainly) no.
  - *Infected* (when the case was initially infected). This is rarely available, but anything (for example an interval) can help.
  - *Arrive* (when the case first arrived in the country/region). This is helpful to narrow down the infection time.
  - *Symptom* (when the case first showed symptom). This is useful because we can impute the infected time if we know the distribution of the incubation period.
3. Make sure to record the URL to your source so everyone can verify.

### Analysis

If you would like to form a team to analyze this dataset, please register your interest [#1](https://github.com/qingyuanzhao/2019-nCov-Data/issues/1). You are also welcome to use it for your own research.

Please use the GitHub Issues to make any long suggestion or discussion. You may want to first read the report of a [preliminary analysis](https://htmlpreview.github.io/?https://github.com/qingyuanzhao/2019-nCov-Data/blob/master/Feb1.html). Some problems that we can all think about include:
- How to impute the missing values [#3](https://github.com/qingyuanzhao/2019-nCov-Data/issues/3).
- How to better model the dynamics (from infection, international arrival, symptom onset, initial medical visit to case confimation) recorded in the dataset [#4](https://github.com/qingyuanzhao/2019-nCov-Data/issues/4).
- How to incorporate the lockdown on 23rd of January in the model [#5](https://github.com/qingyuanzhao/2019-nCov-Data/issues/5).

## Update: Febraury 11th

### Major update to the GitHub repository

- The project has been restructured as a R package.

### Major update to the dataset

- Now more than 300 cases in China.

## Update: February 3rd

### Major update to the dataset

- Included a new dataset for cases in China (suggested by Cindy Chen). Elsa Yang have recorded 76 cases in Hefei that are confirmed by February 2nd.
- Updated most countries to February 3rd.

### New report

Please click [here](https://htmlpreview.github.io/?https://github.com/qingyuanzhao/2019-nCov-Data/blob/master/1st-Report/Feb3.html) to view the report. This analysis fits an exponential growth model to infection time imputed using symptom onset and reported incubation interval. This report **has NOT been peer-reviewed** and extra caution is required to interpret the results.


## Update: February 1st

### Major update to the dataset

- Added cases in Hong Kong and Macau (suggested by Cindy Chen).
- The "Hospital" column has been splitted to "Initial" and "Hospital". The "Initial" column record when the patient first went (or was taken) to an outpatient clinic or emergency room after developing symptoms. The "Hospital" column records if the patient was not admitted immediately during the first visit, when he/she was eventually admitted to an hospital. **This split has only been done for Japan, Singapore, Taiwan, Hong Kong, Macau, South Korea.**

### New preliminary analysis

Please click [here](https://htmlpreview.github.io/?https://github.com/qingyuanzhao/2019-nCov-Data/blob/master/1st-Report/Feb1.html) to view the report. This report **has NOT been peer-reviewed** and extra caution is required to interpret the results.
