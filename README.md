![alt text](https://github.com/BenedettaBanzi/SpaceTrooper_BB/blob/devel/SpaceTrooper_logo.png)
# SpaceTrooper
An R package for the preprocessing and quality control of imaging-based spatial transcriptomics data

### Disclaimer
This is still a development version. A lot of fixing is required, please be patient for any
bug or anomaly.

### References
Manuscript and quoting coming soon!

### Requirements
To install the package, R must be at least 4.3.2.

### Installation
First try to run the following either in R or in a virgin environment:
```
remotes::install_github(repo="https://github.com/BenedettaBanzi/SpaceTrooper_BB", ref = "devel")
```

It is highly probable that it won't finish up because some dependencies' installation is a bit cumbersome.
PLEASE, DON'T GET DISCOURAGED! Just look at the error message and jot down which packages were not successfully
installed. Then install manually the packages either using BioCManager for Bioconductor packages, the standard 
"install.packages()" in R, or through conda or bioconda. For example:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scater")
```
Then try to run again SpaceTrooper installation, now it should work fine and take a couple of minutes.

### Demo
If you were able to install the package correctly, CONGRATULATIONS!!! 🎉🎉🎉
Now try to follow along the "11_15_SpaceTrooper_demo_DBKero" vignette in "vignette" folder.
Otherwise, you can download the html file from the same folder and get a look at that.

Thanks for checking out! 🌸

