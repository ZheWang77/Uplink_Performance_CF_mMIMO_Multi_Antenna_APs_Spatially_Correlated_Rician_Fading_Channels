function [Hhat_EW_MMSE] = functionChannelEstimates_EW_MMSE(R_AP,HMean,H,nbrOfRealizations,M,K,N,tau_p,p,Pset)

%EW_MMSE channel estimator for Cell-Free setup. The estimation is locally
%performed at the APs.
%And each AP is equipped with N antennas.

%%=============================================================
%This function was developed as a part of the paper:
%
%Zhe Wang, Jiayi Zhang, Emil Bj?rnson, and Bo Ai, "Uplink Performance of %Cell-Free Massive MIMO Over Spatially Correlated Rician Fading Channels,"
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
%R_AP                 = Matrix with dimension N x N x M x K x nbrOfRealzations where (:,:,m,k,i) is
%                       the spatial correlation matrix between AP l and UE k 
%                       in i^th realization, normalized by the noise power
%HMean                = Matrix with dimension MN  x K x nbrOfRealzations
%                       where (mn,i,k) is the i^th realization of the channel mean
%                       between the n^th antenna of AP m and UE k (with phase shift)                  
%H                    = Matrix with dimension MN x K x nbrOfRealzations 
%                       where (mn,i,k) is the i^th channel realization
%                       between the n^th antenna of AP m and UE k
%nbrOfRealizations    = Number of realizations
%M                    = Number of APs
%K                    = Number of UEs 
%N                    = Number of AP antennas
%p                    = 1xK vector, uplink power at each UE
%tau_p                = Pilot length
%Pset                 = Pilot allocation set
%
%OUTPUT:
%Hhat_EW_MMSE            = Matrix with dimension MN x nbrOfRealzations x K
%                       where (mn,k,i) is the i^th  realization of EW-MMSE
%                       channel estimate between the n^th antenna of AP m and UE k
%                     


%Prepare to store EW-MMSE channel estimates
Hhat_EW_MMSE = zeros(M*N,nbrOfRealizations,K);

%Store identity matrix of size N x N
eyeN = eye(N);

%Generate realizations of normalized noise 
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,M,K) + 1i*randn(N,nbrOfRealizations,M,K));


for m = 1:M
    for k = 1:K
        
        yp = zeros(N,nbrOfRealizations);
        yMean = zeros(N,nbrOfRealizations);
        PsiInv = zeros(N,N);
        inds = Pset(:,k); 
        
        for z = 1:length(inds)
            
            yp = yp + sqrt(p(inds(z)))*tau_p*H((m-1)*N+1:m*N,:,inds(z));
            yMean = yMean + sqrt(p(inds(z)))*tau_p*HMean((m-1)*N+1:m*N,:,inds(z));
            PsiInv = PsiInv + p(inds(z))*tau_p*diag(diag(R_AP(:,:,m,inds(z))));
            
        end
        yp = yp + sqrt(tau_p)*Np(:,:,m,k);
        PsiInv = PsiInv + eyeN;
        
      
        for z = 1:length(inds)
            
            RPsi = diag(diag(R_AP(:,:,m,inds(z))))/PsiInv;
            Hhat_EW_MMSE((m-1)*N+1:m*N,:,inds(z)) = HMean((m-1)*N+1:m*N,:,inds(z)) + sqrt(p(inds(z)))*RPsi*(yp-yMean);
            
        end
    end
end
