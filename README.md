# A Homebrewed Concotion For Spatial Cluster Dynamics in Vegetation.


**The code in this repository is published under the open-access BSD 3 license. To learn more about what this entails, please refer to [this](https://opensource.org/licenses/BSD-3-Clause).**

This project --- and the code listed here --- has two primary functions. The first is to determine the critical exponents associated with macroscopic observables in the various percolation systems (which we shall elaborate on to an extent in the following sections), while the second is to capture cluster dynamics between two temporal frames of the percolation system, which are spaced a variable time-steps apart.

*For a deeper look into the theory and the particulars of the algorithms, look at the Methods Section of this [document](https://drive.google.com/file/d/17A0xSbOQrUHNOF3O6vys0rPSrMBqjLdn/view?usp=sharing).*

## Null Percolation (NP) & Directed Percolation (DP)

We make use of the **stochastic cellular automata (SCA)**, which are lattice-based discrete-time random processes.

In this model we have two possible states - vegetated or barren - for a given site on an square lattice of side length L, the former represented by 1 and the latter by 0. We consider non-periodic boundaries, although the code also contains functions for periodic boundary simulations.

The transition rules for NP are:

<div align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=0&space;\overset{p}{\rightarrow}&space;1:&space;Birth" target="_blank"><img src="https://latex.codecogs.com/svg.latex?0&space;\overset{p}{\rightarrow}&space;1:&space;Birth" title="0 \overset{p}{\rightarrow} 1: Birth" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=1&space;\overset{1-p}{\rightarrow}&space;0:&space;Death" target="_blank"><img src="https://latex.codecogs.com/svg.latex?1&space;\overset{1-p}{\rightarrow}&space;0:&space;Death" title="1 \overset{1-p}{\rightarrow} 0: Death" /></a>
</div>

A single parameter p taking values from 0 to 1 governs the model. Both birth and death are frequency-dependent.

In DP, we encounter the simplest model for encoding local positive feedback (of Van-Neumann radius 1). Here the transition rules are given by:

<div align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=01&space;\overset{p}{\rightarrow}&space;11:&space;Birth" target="_blank"><img src="https://latex.codecogs.com/svg.latex?01&space;\overset{p}{\rightarrow}&space;11:&space;Birth" title="01 \overset{p}{\rightarrow} 11: Birth" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=1&space;\overset{1-p}{\rightarrow}&space;0:&space;Death" target="_blank"><img src="https://latex.codecogs.com/svg.latex?1&space;\overset{1-p}{\rightarrow}&space;0:&space;Death" title="1 \overset{1-p}{\rightarrow} 0: Death" /></a>
</div>

One may observe for the NP & DP processes in the following figure.

![NP-DP-Schema](/MSThesisFigures/NPDPSchema.png "NP-DP-Schema")

## Tricritical Percolation (TCP)

 The TCP process introduces an additional higher-order parameter (q), that modulates the strength of the local positive feedback.  TCP is a generalisation of DP, being exactly identical in the case q=0. The transition rules are:

<div align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=01&space;\overset{p}{\rightarrow}&space;11:&space;Birth" target="_blank"><img src="https://latex.codecogs.com/svg.latex?01&space;\overset{p}{\rightarrow}&space;11:&space;Birth" title="01 \overset{p}{\rightarrow} 11: Birth" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=1&space;\overset{1-p}{\rightarrow}&space;0:&space;Death" target="_blank"><img src="https://latex.codecogs.com/svg.latex?1&space;\overset{1-p}{\rightarrow}&space;0:&space;Death" title="1 \overset{1-p}{\rightarrow} 0: Death" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=011&space;\overset{q}{\rightarrow}&space;111:&space;Facilitation" target="_blank"><img src="https://latex.codecogs.com/svg.latex?011&space;\overset{q}{\rightarrow}&space;111:&space;Facilitation" title="011 \overset{q}{\rightarrow} 111: Facilitation" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=11&space;\overset{(1-q)(1-p)}{\rightarrow}&space;01:&space;Facilitation" target="_blank"><img src="https://latex.codecogs.com/svg.latex?11&space;\overset{(1-q)(1-p)}{\rightarrow}&space;01:&space;Facilitation" title="11 \overset{(1-q)(1-p)}{\rightarrow} 01: Facilitation" /></a>
</div>

One may observe for the TCP process in the following figure.

![TCP-Schema](/MSThesisFigures/TCPSchema.png "TCP-Schema")


## Critical Exponents

Certain macroscopic observables, such as the average cluster size (S[p]) & percolation strength (P[p]) display a power-law (a power law with size of the system (L)), close to certain value of p, known as the geometric percolation threshold, pc. In the case of NP, we have:
<div align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=P\[p\]&space;\propto&space;L^{-(\beta/\nu)}" target="_blank"><img src="https://latex.codecogs.com/png.latex?P\[p\]&space;\propto&space;L^{-(\beta/\nu)}" title="P\[p\] \propto L^{-(\beta/\nu)}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=S\[p\]&space;\propto&space;L^{(\gamma/\nu)}" target="_blank"><img src="https://latex.codecogs.com/png.latex?S\[p\]&space;\propto&space;L^{(\gamma/\nu)}" title="S\[p\] \propto L^{(\gamma/\nu)}" /></a>
</div>

We set out to examine if the same power law regimes hold for DP & TCP processes. For further details refer [here](https://drive.google.com/file/d/17A0xSbOQrUHNOF3O6vys0rPSrMBqjLdn/view?usp=sharing).

## The Code

All the relevant code for generating raw simulation data (stored as CSVs under `/simulations/.../.csv`) are present in the `/simulations` sub-directory. Further analysis of this raw data can be done by using the Python scripts present in the various sub-directories of `/analysis`. Final plots of this analysed data is then stored under the appropriate sub-directories of `/figures`.

### Finding The Geometric Percolation Threshold (NP/DP/TCP)

*NOTE: THE CODE FOR THIS SECTION MAY BE FOUND IN THE BRANCH "MAIN".*

For finding the geometric percolation threshold for the NP/DP processes, make use of the following command:

```
g++ finding_pc_dp.cpp cluster_dynamics.cpp -fopenmp
```

*NOTE: Much of these scripts make use of OpenMP for parallelising simulations over the individual cores of a CPU.*

A sample input would resemble:

```
Enter grid size: 512
Enter number of random trials for each value of p: 8
Enter number of census (measurements taken in a given random trial for a given p): 15
Enter lag in terms of frames: 0
```

For the TCP process, run the following to determine the percolation threshold.

```
g++ finding_pc_tcp.cpp cluster_dynamics.cpp -fopenmp
```

The input is similar to the above, but with the addition of an input q parameter.

```
Enter grid size: 512
Enter number of random trials for each value of p: 8
Enter number of census (measurements taken in a given random trial for a given p): 15
Enter lag in terms of frames: 0
Enter q-value: 0.1
```

### Finite Scaling For DP/NP Process:

*NOTE: THE CODE FOR THIS SECTION MAY BE FOUND IN THE BRANCH "MAIN".*

Run the following command:

```
g++ finite_scaling_sp.cpp cluster_dynamics.cpp -fopenmp
```

An input parameter for this code is an accurate approximation of the geometric percolation threshold for the NP/DP process. If such an approximation is unavailable, one may obtain the same using the code provided in the previous section. A sample input would resemble:

```
Enter p (must be very, very close to p_c): 0.7285
Enter starting grid size (as a power of 2, such as 5 ~ 32):  4.5
Enter ending grid size (as a power of 2, such as 9 ~ 512): 8.5
Enter number of points at which to sample Grid Length: 33
Enter number of random trials to be inititated per value of the Grid Length: 18
Enter number of censuses per random trial: 14
Enter type of scaling experiment to be performed (choose b/w 'Beta', 'Nu' or 'Gam' (default is Beta)): Gam
```

Note for the DP process, only *'Beta'* & *'Gam'* constitute valid inputs, representing the composite finite-scaling exponents Beta/Nu & Gamma/Nu respectively.

### Finite Scaling For TCP/DP Process:

*NOTE: THE CODE FOR THIS SECTION MAY BE FOUND IN THE BRANCH "MAIN".*

Run the following command:

```
g++ finite_scaling_TCP.cpp cluster_dynamics.cpp -fopenmp
```

An input parameter for this code is an accurate approximation of the geometric percolation threshold for the TCP/DP process. If such an approximation is unavailable, one may obtain the same using the code provided in the previous section. A sample input would resemble:

```
Enter q: 0.03
Enter p (must be very, very close to p_c(q)): 0.72338
Enter starting grid size (as a power of 2, such as 5 ~ 32):  5
Enter ending grid size (as a power of 2, such as 9 ~ 512): 5.875
Enter number of points at which to sample Grid Length: 8
Enter number of random trials to be inititated per value of the Grid Length: 2
Enter number of censuses per random trial: 5
Enter lag between temporal frames: 5000
```

This gives the finite scaling exponents for the DP process in the marginal case q=0. Note that there is no selection of which particular exponent to estimate, by default the raw data that is outputted can be used to estimate both Beta/Nu & Gamma/Nu (this estimation from Stauffer-Aharony Pg. 70-75).


### Multi-shot Transformations For The DP Process:

*NOTE: THE CODE FOR THIS SECTION MAY BE FOUND IN THE BRANCH "ARBITRER".*

The dynamics of change in cluster sizes between (correlated) temporal frames that are a variable time-steps apart. For this, the spatial cross-correlation function between frames had to be determined first for individual values of p (each corresponding to an unique average active density of sites, as shown in the table below), after which these functional forms were plugged in to measure the dynamics of cluster size changes between frames that are variably spaced (with the cluster thus bearing 'multi-shot' transformations between the frames).

To calculate the ij-th index of the Cross-Correlation Matrix, we use the following formula:

<div align=center>

<a href="https://www.codecogs.com/eqnedit.php?latex=r_{i&space;j}=\frac{\sum_{m}&space;\sum_{n}[f(m&plus;i,&space;n&plus;j)-\bar{f}][g(m,&space;n)-\bar{g}]}{\sqrt{\sum_{m}&space;\sum_{n}[f(m,&space;n)-\bar{f}]^{2}&space;\sum_{m}&space;\sum_{n}[g(m,&space;n)-\bar{g}]^{2}}}" target="_blank"><img src="https://latex.codecogs.com/png.latex?r_{i&space;j}=\frac{\sum_{m}&space;\sum_{n}[f(m&plus;i,&space;n&plus;j)-\bar{f}][g(m,&space;n)-\bar{g}]}{\sqrt{\sum_{m}&space;\sum_{n}[f(m,&space;n)-\bar{f}]^{2}&space;\sum_{m}&space;\sum_{n}[g(m,&space;n)-\bar{g}]^{2}}}" title="r_{i j}=\frac{\sum_{m} \sum_{n}[f(m+i, n+j)-\bar{f}][g(m, n)-\bar{g}]}{\sqrt{\sum_{m} \sum_{n}[f(m, n)-\bar{f}]^{2} \sum_{m} \sum_{n}[g(m, n)-\bar{g}]^{2}}}" /></a>

</div>

Here f() & g() represent two temporally distinct frames, with average density of sites given by f-bar and g-bar. We are primarily interested in the max element of this cross-correlation matrix, which is almost always r(0,0) for correlated temporal frames (i.e. max cross-correlation is achieved when the two frames are superimposed on each other).

The max cross-correlation function, dependent on the time difference seperating two frames is calculated seperately for p values corresponding to unit increments of the active density of sites (by 0.01) in the region 0.3-0.7 as shown in the table below (courtsey Ayan Das, TEE Lab, IISc):

![Density Equivalence](/MSThesisFigures/DensityTable.png "Density Equivalence")

To calculate the max cross correlation between temporal frames as a function of the time difference between two frames at a particular value of p, run the following:

```
g++ CrosCol_DP.cpp cluster_dynamics.cpp -fopenmp
```

A typical input would resemble:

```
Enter q: 0.03
Enter p (must be very, very close to p_c(q)): 0.72338
Enter grid size:   256
Enter p value of contention: 0.707
Enter number of points at which to sample Grid Length: 8
Enter number of random trials to be inititated per value of the Grid Length: 2
Enter number of census: 12500
Enter lag: 1
Enter number of random trials: 8
```

Here "number of census" is proportional to the maximum time-step difference upto which the Cross-Correlation function is calculated.

By building a wrapper linking each p-value in the above table to a unique Cross-Correlation functional form (the exponents that describe this form, basically), we can obtain the relevant function by simply inputting the relevant occupancy rate(p). In our case, we replicated this binder through a CSV, where the values on each row corresponded to the functional form at a given occupancy rate (15_16_KungF---U.csv).

Unpacking the rows of this binder CSV, we are ready to compute the multi-shot transformation statistics of cluster dynamics, i.e. examine the probability distribution function of the magnitude of change in the sizes of clusters between two temporal (correlated) frames that are seperated by a significant number of time steps.

The exploration of multi-shot transformations was carried out for temporal frames that were seperated by time-steps that corresponded to Max-Cross-Correlations of 0.95, 0.9, 0.8, 0.75, 0.7 & 0.6.

These transformations can be run by using:

```
g++ MultiTransformersDP.cpp cluster_dynamics.cpp -fopenmp
```
