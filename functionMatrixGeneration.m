function [Rp,Dk,Lk,Phi,Phip,Phi_EW,Omega,Omegap,Sigma,C_MMSE,C_EW_MMSE,C_LMMSE] = functionMatrixGeneration(R_AP,HMean_Withoutphase,M,K,N,tau_p,p,Pset)

%---This function is used to generate the matrix used in the next
%subsequent calculation 

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
%R_AP                 = Matrix with dimension N x N x M x K where (:,:,m,k) is
%                       the spatial correlation matrix between AP l and UE k 
%                       in setup n, normalized by the noise power
%HMeanWithoutPhase    = Matrix with dimension MN x K ,where (mn,k) is the
%                      channel mean between the n^th antenna of AP m and UE k, normalized by
%                      noise power and without random phase shifts                
%H                    = Matrix with dimension MN x nbrOfRealzations x K
%                       where (mn,i,k) is the i^th channel realization
%                       between the n^th antenna of AP m and UE k
%M                    = Number of APs
%K                    = Number of UEs 
%N                    = Number of antennas per AP
%p                    = 1xK vector, uplink power at each UE
%tau_p                = Pilot length
%Pset                 = Pilot allocation set
%
%OUTPUT:
%R                    = Matrix with dimension MN x MN x K where ((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) is
%                       the spatial correlation matrix between AP l and UE k 
%                       normalized by the noise power
%Lk                   = Matrix with dimension MN x MN x K 
%                      
%Rp                   = Matrix with dimension MN x MN x K
%
%Phi                  = Matrix with dimension MN x MN x K
%
%Phip                 = Matrix with dimension MN x MN x K 
%                      
%Omega                = Matrix with dimension MN x MN x K
%
%Omegap               = Matrix with dimension MN x MN x K
%
%C_MMSE               = Matrix with dimension N x N x M x K where (:,:,m,k)
%                       is the covariance matrix of channel error for AP m - UE k based on the
%                       MMSE estimator
%
%C_EW_MMSE            = Matrix with dimension N x N x M x K where (:,:,m,k)
%                       is the covariance matrix of channel error for AP m - UE k based on the
%                       EW-MMSE estimator
%
%C_EW_MMSE            = Matrix with dimension N x N x M x K where (:,:,m,k)
%                       is the covariance matrix of channel error for AP m - UE k based on the
%                       EW-MMSE estimator
%
%C_LMMSE              = Matrix with dimension N x N x M x K where (:,:,m,k)
%                       is the covariance matrix of channel error for AP m - UE k based on the
%                       LMMSE estimator



%Prepare to store the result
Lk = zeros(N,N,M,K);
Phi = zeros(N,N,M,K);
Phi_EW = zeros(N,N,M,K);
Phip = zeros(N,N,M,K);
Omega = zeros(N,N,M,K);
Omegap = zeros(N,N,M,K);
Sigma = zeros(N,N,M,K);
Dk = zeros(N,N,M,K);
C_MMSE = zeros(N,N,M,K);
C_EW_MMSE = zeros(N,N,M,K);
C_LMMSE = zeros(N,N,M,K);


for m = 1:M
    for k = 1:K
        

         Lk(:,:,m,k) = HMean_Withoutphase((m-1)*N+1:m*N,k)*(HMean_Withoutphase((m-1)*N+1:m*N,k))';
         Dk(:,:,m,k) = diag(diag(R_AP(:,:,m,k)));
         
    end
end
 Rp = R_AP + Lk;

 
%Go through all APs          
for m = 1:M
    
    %Go through all UEs
    for k = 1:K
        
        %Compute the UEs indexes that use the same pilot as UE k
        inds = Pset(:,k);
        PsiInv = zeros(N,N);
        PsiInv_EW_MMSE = zeros(N,N);
        PsiInv_LMMSE = zeros(N,N);
        
        %Go through all UEs that use the same pilot as UE k 
        for z = 1:length(inds)   
            
            PsiInv = PsiInv + p(inds(z))*tau_p*R_AP(:,:,m,inds(z));
            PsiInv_EW_MMSE = PsiInv_EW_MMSE + p(inds(z))*tau_p*Dk(:,:,m,inds(z));
            PsiInv_LMMSE = PsiInv_LMMSE + p(inds(z))*tau_p*Rp(:,:,m,inds(z));

        
        end
            PsiInv = PsiInv + eye(N);
            PsiInv_EW_MMSE = PsiInv_EW_MMSE + eye(N);
            PsiInv_LMMSE = PsiInv_LMMSE + eye(N);
            
            for z = 1:length(inds)
                
                Phi(:,:,m,inds(z)) = PsiInv;
                Phi_EW(:,:,m,inds(z)) = PsiInv_EW_MMSE;
                Phip(:,:,m,inds(z)) = PsiInv_LMMSE;
            
            end
            
            Omega(:,:,m,k) = R_AP(:,:,m,k)/PsiInv*R_AP(:,:,m,k);
            Sigma(:,:,m,k) = p(k)*tau_p*Dk(:,:,m,k)/PsiInv_EW_MMSE*PsiInv/PsiInv_EW_MMSE*Dk(:,:,m,k);
            Omegap(:,:,m,k) = Rp(:,:,m,k)/PsiInv_LMMSE*Rp(:,:,m,k);
            
    end
end

%Generate estimation error correlation matrices
for k = 1:K
    
    C_MMSE(:,:,:,k) = R_AP(:,:,:,k) - p(k)*tau_p*Omega(:,:,:,k);
    C_LMMSE(:,:,:,k) = Rp(:,:,:,k) - p(k)*tau_p*Omegap(:,:,:,k);
    
    for m = 1:M
        
        C_EW_MMSE(:,:,m,k) = R_AP(:,:,m,k) - p(k)*tau_p*R_AP(:,:,m,k)/Phi_EW(:,:,m,k)*Dk(:,:,m,k) - p(k)*tau_p*Dk(:,:,m,k)/Phi_EW(:,:,m,k)*R_AP(:,:,m,k) + Sigma(:,:,m,k);
        
    end
    
end
    
    
