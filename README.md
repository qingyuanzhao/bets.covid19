# Data for the 2019-nCoV outbreak

On 30th of January, 2020, The World Health Organization has declared the [2019-nCoV outbreak a Public Health Emergency of International Concern](https://www.who.int/news-room/detail/30-01-2020-statement-on-the-second-meeting-of-the-international-health-regulations-(2005)-emergency-committee-regarding-the-outbreak-of-novel-coronavirus-(2019-ncov)). **The outbreak originated from my hometown, Wuhan.**

I am sharing this dataset I collected from public agencies or news media, containing information about the 2019-nCov cases confirmed outside China. This dataset is free to use and share given that appropriate credits are given (see the [license](./LICENSE.md)).

## To Contribute

I am hoping this can become a collaborative project by people across the world. You can contribute by either updating the dataset or by forming a team to analyze the dataset.

### Dataset
The easiest way to contribute is to suggest comments in [this Google Spreadsheet](https://docs.google.com/spreadsheets/d/1H4MzVxkug2txyzkiDJsGVKB04YveYcsHg9ijuer8clE/edit?usp=sharing). I will monitor the suggestions and update this GitHub every day.

**Currently I need help with:**
1. Japanese-speaking people to verify the information recorded for the Japanese cases.
2. Korean-speaking people to verify the information recorded for the Korean cases.
3. I am having a hard time to obtain detailed information for cases in 
- Australia;
- France;
- Russian;
- Thailand;
- United Kingdom;
- United States;
- Vietnam.

Please only use publicly available information that I can confirm when making suggestions.

If you like to be acknowledged for your contribution, please tell me your (real or nick) name in the comments.

### Analysis

If you would like to form a team to analyze this dataset, please register your interest [here](https://github.com/qingyuanzhao/2019-nCov-Data/issues/1). You are also welcome to use it for your own research.

Please use [Issues](https://github.com/qingyuanzhao/Coronavirus-Data/issues) to make any long suggestion or discussion. You may want to first read the report of a [preliminary analysis](https://htmlpreview.github.io/?https://github.com/qingyuanzhao/2019-nCov-Data/blob/master/Feb1.html). Some problems that we can all think about include:
- How to impute the missing values [#3](https://github.com/qingyuanzhao/2019-nCov-Data/issues/3).
- How to better model the dynamics (from infection, international arrival, symptom onset, initial medical visit to case confimation) recorded in the dataset [#4](https://github.com/qingyuanzhao/2019-nCov-Data/issues/4).
- How to incorporate the lockdown on 23rd of January in the model [#5](https://github.com/qingyuanzhao/2019-nCov-Data/issues/5).

## Update: February 1st

### Major update to the dataset

- Added cases in Hong Kong and Macau (suggested by Cindy Chen).
- The "Hospital" column has been splitted to "Initial" and "Hospital". The "Initial" column record when the patient first went (or was taken) to an outpatient clinic or emergency room after developing symptoms. The "Hospital" column records if the patient was not admitted immediately during the first visit, when he/she was eventually admitted to an hospital. **This split has only been done for Japan, Singapore, Taiwan, Hong Kong, Macau, South Korea.**

### New preliminary analysis

Please click [here](https://htmlpreview.github.io/?https://github.com/qingyuanzhao/2019-nCov-Data/blob/master/Feb1.html) to view the report. This report **has NOT been peer-reviewed** and extra caution is required to interpret the results.
