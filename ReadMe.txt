Data:

The data files are in the MATLAB .dat format.
We include two datasets, i.e., 
(1) the LED data from 
Yang, G. (2007), Life Cycle Reliability Engineering, John Wiley & Sons
and 
(2) the HDD data from
Ye, Z.-S., Wang, Y., Tsui, K.-L., and Pecht, M. (2013), Degradation data analysis using Wiener
processes with measurement errors, IEEE Transactions on Reliability, 62(4), 772-780.

Code:

(1) Real data analysis.
You can open the data file using MATLAB, and a variable named "obs" will appear in the workspace. 
Then, you can run ¡°WienerNIG_L(obs)¡± to obtain the estimate for linear model,
and ¡°WienerNIG(obs)¡± for the power law model.
(2) Simulation.
The code "simu_WienerNIG" is used to replicate the simulation result in Table 1 of Section 4.
The code "simu_compareMCMC" is used to replicate the result in Table 2 and 3 of Section 4, on the comparison of EM and MCEM.
