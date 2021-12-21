# FMI
On the conventional definition of path-specific effects - fully mediated interaction with multiple ordered mediators.

To install this package, please copy and paste following codes into your R session:

1. install.packages("devtools")
2. library(devtools)
3. install_github("AshTai/FMI")



## Example
```R
library(FMI)
FMI(data=sim.data,outcome="feeltherm",treat = "spifev1per",mediators=c("m1_crqdysmnP","m2_hadsanxN"),
baseline.conf = c("countryB","female","ageT0","bmi","smopackyrs"),
med.model=c("gaussian","gaussian"),
outcome.model = "gaussian",nboost=100)
```

## Contact information
An-Shun Tai ([anshuntai@nctu.edu.tw](mailto:anshuntai@nctu.edu.tw))  https://anshuntai.weebly.com
