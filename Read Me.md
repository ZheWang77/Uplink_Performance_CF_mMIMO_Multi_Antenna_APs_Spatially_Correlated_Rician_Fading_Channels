# Uplink Performance of Cell-Free Massive MIMO Over Spatially Correlated Rician Fading Channels

This is a code package is related to the following scientific article:

Zhe Wang, Jiayi Zhang, Emil Björnson, and Bo Ai, "[Uplink Performance of Cell-Free Massive MIMO Over Spatially Correlated Rician Fading Channels](https://ieeexplore.ieee.org/document/9276421)," *IEEE Communications Letters*, vol. 25, no. 4, pp. 1348-1352, April 2021.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*

## Abstract of Article
We consider a practical cell-free massive multiple-input-multiple-output (MIMO) system with multiantenna access points (APs) and spatially correlated Rician fading channels. 
The significant phase-shift of the line-of-sight component induced by the user equipment movement is modeled randomly. 
Furthermore, we investigate the uplink spectral efficiency (SE) with maximum ratio (MR)/local minimum mean squared error (L-MMSE) combining and optimal large-scale fading decoding based on the phase-aware MMSE, phase-aware element-wise MMSE and linear MMSE (LMMSE) estimators. 
Then new closed-form SE expressions with MR combining are derived. Numerical results validate our derived expressions and show that the SE benefits from the spatial correlation. 
It is important to observe that the performance gap between L-MMSE and MR combining increases with the number of antennas per AP and the SE of the LMMSE estimator is lower than that of other estimators due to the lack of phase-shifts knowledge.

## Content of Code Package

The package generates the simulation SE results which are used in Figure 1,  and Figure 3. To be specific:

%--functionChannelEstimates_EW_MMSE: EW-MMSE Estimator

%--functionChannelEstimates_LMMSE: LMMSE Estimator

%--functionChannelEstimates_MMSE: MMSE Estimator

%--functionChannelGeneration: Generate the channel in this paper

%--functionComputeMonteCarloSE_UL: Compute SE with MR combining by the Monte-Carlo method

%--functionComputeSE_AP_MMSE_Combining_UL: Compute SE with L-MMSE combining by the Monte-Carlo method

%--functionMatrixGeneration: Generate matrices used below

%--functionPilotAllocation
%--functionMMSE_interferenceLevels
Pilot allocation

%--functionTheoreticalCellFreeULSE_EW_MMSE
%--functionTheoreticalCellFreeULSE_LMMSE
%--functionTheoreticalCellFreeULSE_MMSE
Compute SE by analytical results with MMSE/EW-MMSE/LMMSE estimators

%--Letter_Fig1: Fig.1 in this paper

%--Letter_Fig3: Fig.3 in this paper

%--RandomAP_generateSetup_Rician_Multi_Antenna: System generation


See each file for further documentation.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
