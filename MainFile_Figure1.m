%Fig.1 in this paper

%%=============================================================
%This function was developed as a part of the paper:
%
%Zhe Wang, Jiayi Zhang, Emil Bjornson, and Bo Ai, "Uplink Performance of Cell-Free Massive MIMO Over Spatially Correlated Rician Fading Channels,"
%IEEE Communications Letters, vol. 25, no. 4, pp. 1348-1352, April 2021, %doi: 10.1109/LCOMM.2020.3041899.
%
%Download article: https://ieeexplore.ieee.org/document/9276421 or https://arxiv.org/abs/2110.05796
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.
%This is version 1.0 (Last edited: 2020-05-12)
%%=============================================================

%Empty workspace and close figures
clear
close all
clc

%Define the range of number of Access Points (APs)
M = 40:10:100;

 

%Number of UEs
K = 20;

%Number of antennas per AP
N = [1,2,4];

%Pilot length
tau_p = 10;


%Select length of coherence block
tau_c = 200;


%Uplink transmit power per UE (W)
p = 0.2; %200 mW
%Create the power vector for all UEs (The uplink power is the same
%(p)at each UE)
pv = p*ones(1,K );


%Select the number of setups with random AP/UE locations
nbrOfSetups = 50;
%Select the number of channel realizations per setup
nbrOfRealizations = 1000;


%Prepare to save simulation results 
userSE_MMSE_LSFD_MR_Combining= zeros(length(M),nbrOfSetups,length(N));
userSE_MMSE_LSFD_MMSE_Combining= zeros(length(M),nbrOfSetups,length(N));



%Prepare to save simulation results 
SE_MMSE_LSFD_MR_Combining= zeros(K,nbrOfSetups,length(M),length(N));
SE_MMSE_LSFD_MMSE_Combining= zeros(K,nbrOfSetups,length(M),length(N));


       
%Go through the number of APs
for i = 1:length(N)
for m = 1:length(M)

    A_singleLayer = reshape(repmat(eye(M(m)),1,K),M(m),M(m),K);
    for n=1:nbrOfSetups
       %Deploy UEs and generate the covariance and mean matrices
       [R_AP,HMean_Withoutphase,~] = RandomAP_generateSetup_Rician_Multi_Antenna(M(m),K,N(i),1,1);
  
       %Create channel generations for each UE-AP pair
        [H,HMean] = functionChannelGeneration( R_AP,HMean_Withoutphase,M(m),K,N(i),nbrOfRealizations);
        disp(['ChannelGeneration of setup ' num2str(n)]);

       [Pset] = functionPilotAllocation( R_AP,HMean_Withoutphase,A_singleLayer,M(m),K,N(i),tau_p,pv);
       disp(['PilotAllocation of setup ' num2str(n)]);

       [Rp,Dk,Lk,Phi,Phip,Phi_EW,Omega,Omegap,Sigma,C_MMSE,C_EW_MMSE,C_LMMSE] = functionMatrixGeneration(R_AP,HMean_Withoutphase,M(m),K,N(i),tau_p,pv,Pset);
       disp(['MatrixGeneration of setup ' num2str(n)]);

       [Hhat_MMSE] = functionChannelEstimates_MMSE(R_AP,HMean,H,nbrOfRealizations,M(m),K,N(i),tau_p,pv,Pset);
       disp(['MMSE Estimate of setup ' num2str(n)]);

       %Second layer decoding Spectral Efficiency (SE) computation
       %SE with MMSE estimator

       [SE_MR_MMSE] = functionComputeMonteCarloSE_UL(Hhat_MMSE,H,tau_c,tau_p,nbrOfRealizations,N(i),K,M(m),pv);
       disp(['MR Combining Monte Compution of MMSE Estimate of setup ' num2str(n)]);
       [SE_MMSE_Combining_MMSE] = functionComputeSE_AP_MMSE_Combining_UL(Hhat_MMSE,H,C_MMSE,tau_c,tau_p,nbrOfRealizations,N(i),K,M(m),pv);
       disp(['MMSE Combining Monte Compution of MMSE Estimate of setup ' num2str(n)]);


       
       %--Average SE
       userSE_MMSE_LSFD_MR_Combining(m,n,i) = mean(SE_MR_MMSE,1);
       userSE_MMSE_LSFD_MMSE_Combining(m,n,i) = mean(SE_MMSE_Combining_MMSE,1);

      
       %--SE
       SE_MMSE_LSFD_MR_Combining(:,n,m,i) = SE_MR_MMSE;
       SE_MMSE_LSFD_MMSE_Combining(:,n,m,i)= SE_MMSE_Combining_MMSE;

       
       disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    end 
    clear R_AP H HMean Hhat_MMSE Hhat_LMMSE 

    disp([num2str(M(m)) ' APs out of ' num2str(M(end))]);
end
end


