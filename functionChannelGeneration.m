function [H,HMean] = functionChannelGeneration( R_AP,HMean_Withoutphase,M,K,N,nbrOfRealizations)

%---This function is used to generate the channel realizations for given covariance and
%mean matrices. The outputs are channel realizations and channel means with
%random phase shifts at a coherence block.
%And each AP is equipped with N antennas.

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

%INPUT:
%R_AP                  = Matrix with dimension N x N x M x K x nbrOfRealizations where (:,:,m,i,k) is 
%                        the spatial correlation matrix between AP m and UE k in i^th channel realization, normalized by noise
%                        power
%M                     = Number of APs
%nbrOfRealizations     = Number of realizations
%K                     = Number of UEs 
%N                     = Number of antennas per AP
%
%OUTPUT:
%
%H                     = Matrix with dimension MN x K x nbrOfRealzations
%                        where (mn,k,i) is the i^th channel realization
%                        between the n^th antenna of AP m and UE k
%                     
%HMean                 = Matrix with dimension MN x K x nbrOfRealzations 
%                        where (mn,k,i) is the i^th realization of the channel mean
%                        between the n^th antenna of AP m and UE k


%Prepare to store the results   
R = zeros(M*N,M*N,K);
W = (randn(M*N,nbrOfRealizations,K)+1i*randn(M*N,nbrOfRealizations,K));
H = zeros(M*N,nbrOfRealizations,K);

%Same channelGain_LoS and channelGain_NLoS for all realizations (at each setup) but phases shift of LoS are different at each
%coherence block
HMean = zeros(M*N,nbrOfRealizations,K); 
HMeanx = reshape(repmat(HMean_Withoutphase,nbrOfRealizations,1),M*N,nbrOfRealizations,K); 


%--Phase shift same for all antennas
angles = -pi + 2*pi*rand(M,nbrOfRealizations,K);
phaseMatrix = exp(1i*angles);
v_kron = ones(N,1);
PhaseMatrix = zeros(M*N,nbrOfRealizations,K);

for k = 1:K
    PhaseMatrix(:,:,k) = kron(phaseMatrix(:,:,k),v_kron);
end
    
    


for m = 1:M
    for k = 1:K
        
        R((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) = R_AP(:,:,m,k);
        
    end
end


%Go through all UEs and apply the channel gains to the spatial
%correlation and mean matrices and introduce the phase shifts 
for k = 1:K
    
    HMean(:,:,k)= PhaseMatrix(:,:,k).*HMeanx(:,:,k);
    Rsqrt = sqrtm(R(:,:,k));
    H(:,:,k) = sqrt(0.5)*Rsqrt*W(:,:,k) + HMean(:,:,k);
       
end








