% This Matlab script generates Figure 2 in the paper:
%
% Emil Björnson, Özgecan Özdogan, Erik G. Larsson, “Reconfigurable
% Intelligent Surfaces: Three Myths and Two Critical Questions,” IEEE
% Communications Magazine, vol. 58, no. 12, pp. 90-96, December 2020.
%
% Download article: https://arxiv.org/pdf/2006.03377.pdf
%
% This is version 1.0 (Last edited: 2021-01-02)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.


close all;
clear;

%Wavelength (3 GHz)
lambda = 0.1;

%Propagation distances and angles for source and destination
d = 300;
delta = 10;
eta = pi/3;
omega = 0;


%Transmit powers (computed over 20 MHz, to get representative numbers)
Ptx = 10000; %mW
Prelay = 100; %mW

%Bandwidth
B = 20e6;

%Noise figure (in dB)
noiseFiguredB = 10;

%Compute the noise power in dBm
sigma2dBm = -174 + 10*log10(B) + noiseFiguredB;
sigma2 = db2pow(sigma2dBm);

%Transmit SNR values
Ptxsigma2 = Ptx/sigma2;
Prelaysigma2 = Prelay/sigma2;

%Area of isotropic antenna (used for comparison)
A = lambda^2/(4*pi);

%Number of RIS elements (of size equivalent to isotropic antenna)
NvaluesNonInteger = logspace(0,6.5,200);
Nsqrt = ceil(sqrt(NvaluesNonInteger));
Nvalues = Nsqrt.^2;

%Side length of each element
a = sqrt(A);

%Set the antenna gain at the source (10 dBi)
antennaGainTx = 10;

%Set the extra propagation loss through the window (20 dB)
penetrationLoss = 100;

%Define locations of the source and destination
p_t = [d*sin(eta); 0; d*cos(eta)];
p_r = [delta*sin(omega); 0; delta*cos(omega)];


%Prepare to save simulation results
SNR_RIS = zeros(length(Nvalues),1);
SNR_DFrelay = zeros(length(Nvalues),1);


%% Go through the different number elements
for j = 1:length(Nvalues)
    
    %Extract the number of elements/antennas
    N = Nvalues(j);
    
    %Prepare to store channel gains for individual elements/antennas
    betaHn = zeros(N,1);
    betaGn = zeros(N,1);
    
    %Go through each element/antenna and compute pathlosses
    for n = 1:N
        
        %Compute location using Eqs. (22)-(23) in [11]
        x = -a*(sqrt(N)-1)/2 + a*mod(n-1,sqrt(N));
        y = a*(sqrt(N)-1)/2 - a*floor((n-1)/sqrt(N));
        
        %Compute channel gain for the n:th element
        betaHn(n) = channelgainGeneral(p_t,[x; y; 0],a);
        betaGn(n) = channelgainGeneral(p_r,[x; y; 0],a);
        
    end
    
    %Compute the SNR with the RIS using Eq. (42) in [11]
    SNR_RIS(j) = (Ptxsigma2*antennaGainTx/penetrationLoss) * sum(sqrt(betaHn.*betaGn)).^2;
    
    %Compute the SNR with the DF relay using Eq. (37) in [11]
    SNR_DFrelay(j) = min([(Ptxsigma2*antennaGainTx/penetrationLoss)*channelGainArray(d,eta,N,A) Prelaysigma2*channelGainArray(delta,omega,N,A)]);
    
end


%Compute the spectral efficiencies in the different setups
rate_RIS = log2(1+SNR_RIS);
rate_DFrelay = (1/2)*log2(1+SNR_DFrelay);


%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(Nvalues*A,rate_RIS,'b-','LineWidth',2);
plot(Nvalues*A,rate_DFrelay,'k--','LineWidth',2);
set(gca,'XScale','log');
xlabel('Surface area [m$^2$]','Interpreter','Latex');
ylabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
legend({'RIS','DF relay'},'Interpreter','Latex','Location','SouthEast');
set(gca,'fontsize',18);
xlim([1e-3 1e3]);
ylim([0 14]);
xticks([1e-3 1e-2 1e-1 1 1e1 1e2 1e3]);
