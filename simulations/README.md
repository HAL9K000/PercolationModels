## Vegetation Cover 

*What is the vegetation cover (the fraction of 1s) for some given parameter(s)?* 

For null percolation, the answer is simple. There are no interactions and a trivial mean-field analysis will reveal that the vegetation cover for a given value of p is simply p itself. For DP and TP, the script `densities_tp.cpp` may be used to determine the vegetation cover. From the current directory (`cwd`) compile as:

```shell
g++ densities_tp.cpp cluster_dynamics.cpp -fopenmp 
```
A executable `a.out` will be generated in the `cwd` which may be run by typing `./a.out`. The script will ask you for a number of inputs. The prompts along with a typical set of responses are mentioned below. 

```
Enter the number of slaves you want: 4
Enter grid size: 128
Enter 4 birth probabilities: 0.724
0.728
0.731
0.735
Enter feedback strength: 0.0 
Enter updates per site: 25000
Enter number of census: 1000
Enter lag in terms of frames: 10
Do you want to collect frames? Answer with 1 and 0: 1
```
The above specification generates four slave processes each of which performs a simulation of thick percolation in parallel for the specified grid size, birth probability and feedback strength. Each birth probability must be followed by an return stroke. `updates per site` specifies the duration for the simulation to reach its steady state. For DP, a value of 25000 works well. For TP, much larger values may be required depending on the parameters. `Number of census` specifies the number of frames to be sampled from the steady state to calculate the average vegetation cover. `Lag` species the average updates per site between consecutive samples. If you want to collect the frames from which the average cover was calculated, you may specify `1` for the last prompt. 

After the calculations are done, the terminal will print the following output:

```
p: 0.728 q: 0 Mean: 0.56 Standard Deviation: 0.0059

p: 0.724 q: 0 Mean: 0.55 Standard Deviation: 0.0058

p: 0.731 q: 0 Mean: 0.57 Standard Deviation: 0.0055 

p: 0.735 q: 0 Mean: 0.58 Standard Deviation: 0.0053

CPU Time: 90 seconds 
```

If you opted for collecting the frames, they may be found in the `dump` folder in `cwd`. They are labelled as follows (only shown for p=0.725):

```
tp_frame_128_0_0.725_1000_999.txt
tp_frame_128_0_0.725_1000_998.txt
...
tp_frame_128_0_0.725_1000_0.txt
```
Keep in mind that for large grid sizes these files can become very large. Because DP is a special case of TP when q=0, this script can calculate vegetation cover for both TP and DP: the above specification actually gives the densities for DP. Below is a table giving the vegetation cover and associated value of p for the DP model. These were calculated by averaging over 10000 frames at the steady state. The uncertainty in cover is 0.01 for all values.  

| Cover | p     | Cover | p     | Cover | p     | Cover | p     | Cover | p     | Cover | p     |
|:-----:|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| 0.20  | 0.638 | 0.30  | 0.656 | 0.40  | 0.678 | 0.50  | 0.707 | 0.60  | 0.743 | 0.70  | 0.788 |
| 0.21  | 0.640 | 0.31  | 0.658 | 0.41  | 0.681 | 0.51  | 0.710 | 0.61  | 0.747 | 0.71  | 0.793 |
| 0.22  | 0.642 | 0.32  | 0.660 | 0.42  | 0.684 | 0.52  | 0.714 | 0.62  | 0.751 | 0.72  | 0.798 |
| 0.23  | 0.643 | 0.33  | 0.662 | 0.43  | 0.686 | 0.53  | 0.717 | 0.63  | 0.755 | 0.73  | 0.803 |
| 0.24  | 0.645 | 0.34  | 0.664 | 0.44  | 0.689 | 0.54  | 0.721 | 0.64  | 0.760 | 0.74  | 0.808 |
| 0.25  | 0.647 | 0.35  | 0.666 | 0.45  | 0.692 | 0.55  | 0.724 | 0.65  | 0.764 | 0.75  | 0.814 |
| 0.26  | 0.648 | 0.36  | 0.668 | 0.46  | 0.695 | 0.56  | 0.728 | 0.66  | 0.769 | 0.76  | 0.820 |
| 0.27  | 0.650 | 0.37  | 0.671 | 0.47  | 0.698 | 0.57  | 0.731 | 0.67  | 0.773 | 0.77  | 0.825 |
| 0.28  | 0.652 | 0.38  | 0.673 | 0.48  | 0.701 | 0.58  | 0.735 | 0.68  | 0.778 | 0.78  | 0.831 |
| 0.29  | 0.654 | 0.39  | 0.676 | 0.49  | 0.704 | 0.59  | 0.739 | 0.69  | 0.783 | 0.79  | 0.837 |

## Probability of Percolation 

If we tune the parameters of these models in a way that the number of vegetated sites in the model slowly increases, a transition is observed from a phase
with disconnected clusters of various sizes to a phase that has a single lattice-spanning cluster - the so-called spanning cluster. Such a transition occurs at the percolation point - a point in parameter space where the probability that an arbitrary vegetated site belongs to the spanning cluster (also called the percolation probability) shows a transition from zero to one. The percolation transition is essentially a density-dependent one.

After a model's simulation reaches the steady state, the lattice may or may not percolate depending on the parameters. Further, for the same parameter, it may percolate for some frames and not for others. From a sample of frames from the steady state, the fraction of frames that percolate gives an estimate of the probability of percolation P. The scripts `percolation_probabilities_np.cpp` and `percolation_probabilities_dp.cpp` allow us to calculate this quantity. For DP, compile as

```
g++ percolation_probabilities_dp.cpp cluster_dynamics.cpp -fopenmp 
```
and execute as `./a.out`. The prompts and some typical responses are:

```
Enter grid size: 32
Enter starting p: 0.717
Enter ending p: 0.727
Enter divisions: 11
Enter number of census: 10000
Enter lag in terms of frames: 10 
```

Reliable estimates for smaller grid sizes require averaging over a larger number of frames (`number of census`). A value of 25000 usually works well. After the calculations are done, the terminal will print:

```
p: 0.717 P: 0.144
p: 0.718 P: 0.149
p: 0.719 P: 0.176
p: 0.720 P: 0.187
p: 0.721 P: 0.21
p: 0.722 P: 0.228
p: 0.723 P: 0.239
p: 0.724 P: 0.264
p: 0.725 P: 0.275
p: 0.726 P: 0.296
p: 0.727 P: 0.309

CPU Time: 216 seconds 
```

The procedure is exactly the same for NP for which the script is `percolation_probabilities_np.cpp`. We don't have a script for TP yet but it is straightforward to implement should the user require it. 

## Cluster Size Distribution 

A group of vegetated sites that share nearest neighbours form a cluster. To obtain the size of all the clusters in a given number frames from the steady state, the user may use the script `cluster_sizes_np.cpp` and `cluster_sizes_dp.cpp` Compile and execute as

```
g++ cluster_sizes_dp.cpp cluster_dynamics.cpp 
./a.out 
```
The prompt and some typical responses are given next. 

```
Enter grid size: 512
Enter birth probability: 0.724
What is the corresponding density? 0.55
Enter updates per site: 25000
Enter number of census: 1000
Enter lag between frames: 10
Which replicate is this?: 1
Do you want to collect frames? Answer with 1 or 0: 0 
```
After completion, the output - a text file of whitespace separated integers - can be found in `dump`. 

```
dp_cluster_sizes_512_0.55_1000_1.txt
```
## Average Cluster Size 

Use scripts `average_cluster_size_np.cpp` and `average_cluster_size_dp.cpp` to find the average cluster size S. 

```
g++ average_cluster_size_np.cpp cluster_dynamics.cpp -fopenmp
./a.out 
```

```
Enter grid size: 64
Enter starting p: 0.55
Enter ending p: 0.65
Enter divisions: 11
Enter number of census: 100
Enter lag in terms of frames: 10 
```
As output, the terminal prints

```
p: 0.55 S: 11.2
p: 0.56 S: 12.3
p: 0.57 S: 14
p: 0.58 S: 14.2
p: 0.59 S: 14.1
p: 0.6 S: 12.2
p: 0.61 S: 10.8
p: 0.62 S: 9.12 
p: 0.63 S: 6.85 
p: 0.64 S: 4.95
p: 0.65 S: 3.69 

CPU Time: 74 seconds 
```


## Cluster Dynamics 

The clusters can grow, shrink, merge with other clusters, and split. We use a method suggested in Seri et al. (2012) to track the changes in size of clusters. This is done by `transformations_dp.cpp` and `transformations_np.cpp`. 

```
g++ transformations_dp.cpp cluster_dynamics.cpp 
./a.out 
```

```
Enter grid size: 512
Enter birth probability: 0.724
What is the corresponding density? : 0.55
Enter updates per site: 25000
How many transformations do you want? Enter a number: 2000000
Which replicate is this? : 1
```
After completion, the output is found in `dump`. 

```
dp_transformations_512_0.55_2000000_1.txt
```

The first and second column of the file give the size of the cluster before and after the update respectively. 

```
34922 34928
22286 22285
34928 35063
35063 35062
3964 3963
22285 22290
22290 22291
...
```
