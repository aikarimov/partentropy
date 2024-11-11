# partentropy
Calculate partitioned entropy for a signal. The main idea is taken from the publication [1]. 
Typical use is as follows:

```matlab
load t, Y %load time and Y of the  experimental data
k = 64; %number of clusters
taumax = 3; %max time for the experiment
ntimes = 20; %number of max times of leaving cluster for calculating the slope hKS
hKS = partentropy(t,Y,k,taumax,ntimes)
```

The code estimated the partitioned entropy over `k` clusters, and then estimates the Kolmogorov-Sinai entropy as the slope of the linear fit over first `ntimes*max(tauleave)` where `tauleave` is the time for a signal point to leave the initial cluster.

![Illustration](https://github.com/aikarimov/partentropy/blob/main/hP.jpg)

## Literature
1. Shiozawa, Kota, and Isao T. Tokuda. "Estimating Kolmogorov–Sinai entropy from time series of high-dimensional complex systems." Physics Letters A 510 (2024): 129531.
