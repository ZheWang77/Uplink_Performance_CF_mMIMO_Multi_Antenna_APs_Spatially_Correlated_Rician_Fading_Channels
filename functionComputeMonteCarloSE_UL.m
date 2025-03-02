function [SE_MR] = functionComputeMonteCarloSE_UL(Hhat,H,tau_c,tau_p,nbrOfRealizations,N,K,M,p)
%Compute SE with MR combining by the Monte-Carlo method

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

%
%INPUT:
%Hhat              = Matrix with dimension M*N x nbrOfRealizations x K
%                    where (:,i,k) is the estimated collective channel from
%                    all APs to UE k at channel realization i.
%H                 = Matrix with dimension M*N x nbrOfRealizations x K
%                    where (:,i,k) is the true collective channel from all
%                    APs to UE k at channel realization i.
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
%SE_MR     = K x 1 vector where the (k,1):th element is the uplink SE of 
%            UE k achieved with MR combining




%If only one transmit power is provided, use the same for all the UEs
if length(p) == 1
   p = p*ones(K,1);
end


%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_MR = zeros(K,1);



%Diagonal matrix with transmit powers and its square root
Dp12 = diag(sqrt(p));



signal_MR = zeros(M,K);
scaling_MR = zeros(M,K);
Gp_MR = zeros(M,M,K);
A = zeros(M,M,K);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    
    
    %-----------------Levels 1-3
    gp_MR = zeros(M,K,K);
    
    
    %Go through all APs
    for m = 1:M
        
        %Extract channel realizations from all UEs to AP m
        Hallj = reshape(H(1+(m-1)*N:m*N,n,:),[N K]);
        
        %Extract channel estimate realizations from all UEs to AP m
        Hhatallj = reshape(Hhat(1+(m-1)*N:m*N,n,:),[N K]);
        
        
        %Compute MR combining
        V_MR = Hhatallj;

        
        %Go through all UEs
        for k = 1:K
            
            %%MR combining
            v = V_MR(:,k); %Extract combining vector
            signal_MR(m,k) = signal_MR(m,k) + (v'*Hallj(:,k))/nbrOfRealizations; % (v_mk)'h_mk
            gp_MR(m,:,k) = gp_MR(m,:,k) + (v'*Hallj)*Dp12;                       % 
            scaling_MR(m,k) = scaling_MR(m,k) + norm(v).^2/nbrOfRealizations;    % ||v_mk||^2
            
            
        end
        
    end
    
    %Compute averaging of terms for Level 2 and Level 3
    for k = 1:K
        
        Gp_MR(:,:,k) = Gp_MR(:,:,k) + gp_MR(:,:,k)*gp_MR(:,:,k)'/nbrOfRealizations;
        
    end
    
end


%Compute the SE
for k = 1:K
    
    %With MR combining
    b = signal_MR(:,k);
    A(:,:,k) = Gp_MR(:,:,k) + diag(scaling_MR(:,k)) - p(k)*(b*b');
    SE_MR(k) = prelogFactor*real(log2(1+p(k)*b'*(A(:,:,k)\b)));   
 

end