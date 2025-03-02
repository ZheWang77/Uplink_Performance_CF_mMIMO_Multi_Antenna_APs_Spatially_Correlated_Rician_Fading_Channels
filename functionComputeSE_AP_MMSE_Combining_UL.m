function [SE_MMSE] = functionComputeSE_AP_MMSE_Combining_UL(Hhat,H,C,tau_c,tau_p,nbrOfRealizations,N,K,M,p)
%Compute SE with L-MMSE combining by the Monte-Carlo method

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
%Hhat              = Matrix with dimension M*N x nbrOfRealizations x K
%                    where (:,i,k) is the estimated collective channel from
%                    all APs to UE k at channel realization i.
%H                 = Matrix with dimension M*N x nbrOfRealizations x K
%                    where (:,i,k) is the true collective channel from all
%                    APs to UE k at channel realization i.
%H                 = Matrix with dimension M*N x nbrOfRealizations x K
%                    where (:,i,k) is the true collective channel from all
%                    APs to UE k at channel realization i.
%C                 = Covariance matrix of channel error with N x N x M x K
%                    where (:,:,m,k) is the covariance matrix of channel error for AP m - UE k
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Total number of UEs
%M                 = Number of APs
%p                 = Matrix K x 1 where element k is the uplink transmit
%                    power of UE k (If it is a scalar, the same value is
%                    used for all users)

%
%OUTPUT:
%SE_MMSE     = K x 1 vector where the (k,1):th element is the uplink SE of 
%            UE k achieved with L-MMSE combining



%If only one transmit power is provided, use the same for all the UEs
if length(p) == 1
   p = p*ones(K,1);
end

%If no specific Level 1 transmit powers are provided, use the same as for
%the other levels
if nargin<12
    p1 = p;
end


%Store identity matrices of different sizes
eyeN = eye(N);

%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_MMSE = zeros(K,1);

%Compute sum of all estimation error correlation matrices at every BS
% C_tot = zeros(N,N,M);
C_tot = sum(C,4);
% for k = 1:K
%     
%     C_tot = C_tot + p(k)*(R(:,:,:,k)-B(:,:,:,k));
% 
% end


%Diagonal matrix with transmit powers and its square root
Dp = diag(p);
Dp12 = diag(sqrt(p));


%Prepare to save simulation results

signal_MMSE = zeros(M,K);
scaling_MMSE = zeros(M,K);
G_MMSE = zeros(M,M,K);
A = zeros(M,M,K);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    
    
    %-----------------Levels 1-3
    gp_MMSE = zeros(M,K,K);
    
    
    %Go through all APs
    for m = 1:M
        
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H(1+(m-1)*N:m*N,n,:),[N K]);
        
        %Extract channel estimate realizations from all UEs to AP l
        Hhatallj = reshape(Hhat(1+(m-1)*N:m*N,n,:),[N K]);
        
        
        %Compute MR combining
        V_MMSE = ((Hhatallj*Dp*Hhatallj')+C_tot(:,:,m)+eyeN)\(Hhatallj*Dp);
        

        %Go through all UEs
        for k = 1:K
            
           
            
            
            %%MMSE combining
            v = V_MMSE(:,k); %Extract combining vector

            signal_MMSE(m,k) = signal_MMSE(m,k) + (v'*Hallj(:,k))/nbrOfRealizations;
            gp_MMSE(m,:,k) = gp_MMSE(m,:,k) + (v'*Hallj)*Dp12;
            scaling_MMSE(m,k) = scaling_MMSE(m,k) + norm(v).^2/nbrOfRealizations;
            
            
   
           
        end
        
    end
    
    for k = 1:K
        
        G_MMSE(:,:,k) = G_MMSE(:,:,k) + gp_MMSE(:,:,k)*gp_MMSE(:,:,k)'/nbrOfRealizations;
        
    end
    
end



for k = 1:K
    

    %With L-MMSE combining
    b = signal_MMSE(:,k);
    A(:,:,k) = G_MMSE(:,:,k) + diag(scaling_MMSE(:,k)) - p(k)*(b*b');
    SE_MMSE(k) = prelogFactor*real(log2(1+p(k)*b'*(A(:,:,k)\b)));  
    


end
