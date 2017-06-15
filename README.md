This README file introduces *HRitools*, a package created by Héctor Puigdomènech Gómez. *HRitools* is a tool created in order to analyse adaptation and recombination data with the purpose of **quantifying Hill-Robertson interference** (HRi) with a curvilinear model, as suggested by [Castellano et al. (2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4794616/).

# Installation

In order to use *HRitools* in your R session, it must be installed. Dependencies must be installed previously.

```R
install.packages(c("nls2","curl","rtacklayer","ggplot2"))
```

Once dependencies have been installed, *devtools* is the package required to install *HRitools*, as well as all the packages which are installed from Github.

```R
install.packages("devtools")
library(devtools)
install_github("hectorpuigdo/HRitools")
```

Then, the library can simply be loaded with the comand `library(HRitools) `.

# Usage

Three functions are included within the package: **HRi**, **LVNLtest** and **rhokbPopFly**. Each function has its own documentation, which can be found by, for example, `?HRi`.

## HRi

It is the primordial function of the package. It allows to calculate HRi curvilinear functions over an adapatation-recombination dataset ponderating with 0-fold sites. It simply needs a data frame with recombination, adaptation and 0-fold sites data points. Its usage is as simple as:

```R
dataset <- data.frame(recombination,adapatation,fold0sites)
HRi(dataset)

#results can be stored in an object
results <- HRi(dataset)

#the function results contain several results, such as graphs, forumulas, and vectors
results$Graph
```

## LVNLtest

This function is designed to calculate HRi curvilinear and linear models in order to compare them using Akaike information criterion. The same dataframe used with the function `HRi` can be reused with this function such as:

```R
LVNLtest(dataset)
```

This function provides three objects as a result: a data frame with AIC results, and the results of generating both linear and curvilinear models.

## rhokbPopFly

This function is very useful for researchers who do not have the possibility to get exprimental *Drosophila melanogaster* population recombination data and need to find some population-scaled recombination data. [Hervas et al. (2017)](https://doi.org/10.1093/bioinformatics/btx301) created [PopFly, the *Drosophila* population genomics browser](popfly.uab.cat), and it contains a lot of recombination data from different populations, so `rhokbPopFly` is a function to download and put recombination data into a data frame.

Population and window size must be specified, otherwise the Raleigh —RAL— population recombination dataset with a window size of 100kb is downloaded.

```R
recombination_data <- rhokbPopFly("ZI","10kb")
```
----

Happy R coding!

![](https://media.giphy.com/media/3oEdvdmg0utG4LVBrW/giphy.gif)
