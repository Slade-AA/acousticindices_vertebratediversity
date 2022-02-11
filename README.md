
# acousticindices\_vertebratediversity

Data and analysis scripts for paper looking at correlation between
various acoustic indices and vertebrate diversity

Six different acoustic indices were used: ACI, ADI, AEI, BI, NDSI and
EVN?

## To do

  - Rerun analyses with recleaned dataset (waiting on Seb and Don)
  - Check whether I should be running glmer for some analyses -
    fitdistrplus (zero-inflated beta distribution for frogs?)

# Table of Contents

[Species richness](#species-richness)  
[Shannon diversity index](#shannon-diversity-index)

## Results

### Species richness

Correlation between bird species richness and various acoustic indices
![Alt text](outputs/figures/birds_richness.png)

Correlation between frog species richness and various acoustic indices
![Alt text](outputs/figures/frogs_richness.png)

### Shannon diversity index

Correlation between shannon diversity index for birds and various
acoustic indices ![Alt text](outputs/figures/birds_shannon.png)

Correlation between shannon diversity index for frogs and various
acoustic indices ![Alt text](outputs/figures/frogs_shannon.png)

### Total abundance

Correlation between total bird abundance and various acoustic indices
![Alt text](outputs/figures/birds_count.png)

Correlation between total frog abundance and various acoustic indices
![Alt text](outputs/figures/frogs_count.png)

### Bird diversity vs acoustic indices (morning, day, all)

There was little difference in the predictive performance of random
forest models on bird biodiversity when using acoustic indices
calculated for all recordings, morning recordings (6am-9am), or day
recordings (6am-6pm). ![Alt
text](outputs/figures/randomforestperformance/Birds_MorningDayAll.png)
