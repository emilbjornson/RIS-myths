% This Matlab script generates Figure 4 in the paper:
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
eta = pi/3;
deltaValues = linspace(0.5,50,200);

%Area of each element
A = (lambda/5)^2;

%Number of RIS elements
Nsqrt = 20*5;
N = Nsqrt.^2;

%Side length of each element
a = sqrt(A);

%Set the antenna gain at the source (10 dBi)
antennaGainTx = 10;

%Set the extra propagation loss through the window (20 dB)
penetrationLoss = 100;

%Define location of the source
p_t = [d*sin(eta); 0; d*cos(eta)];


%Prepare to save simulation results
channelGain_RIS = zeros(length(deltaValues),1);
channelGain_RIS_approx = zeros(length(deltaValues),1);
channelGain_mirror = zeros(length(deltaValues),1);


%% Go through the different number elements
for j = 1:length(deltaValues)
    
    %Position of the destination
    p_r = [0; 0; deltaValues(j)];
    
    %Prepare to store channel gains for individual elements/antennas
    betaHn = zeros(N,1);
    betaGn = zeros(N,1);
    phaseHn = zeros(N,1);
    phaseGn = zeros(N,1);
    
    %Go through each element/antenna and compute pathlosses
    for n = 1:N
        
        %Compute location using Eqs. (22)-(23) in [11]
        x = -a*(sqrt(N)-1)/2 + a*mod(n-1,sqrt(N));
        y = a*(sqrt(N)-1)/2 - a*floor((n-1)/sqrt(N));
        
        %Compute channel gain for the n:th element
        betaHn(n) = channelgainGeneral(p_t,[x; y; 0],a);
        betaGn(n) = channelgainGeneral(p_r,[x; y; 0],a);
        
        %Compute phase-shift for the n:th element
        phaseHn(n) = mod(norm(p_t-[x; y; 0])/lambda,1)*2*pi;
        phaseGn(n) = mod(norm(p_r-[x; y; 0])/lambda,1)*2*pi;
        
    end
    
    
    %Compute the total channel gain with the RIS using Eq. (42) in [11],
    %by removing the P/sigma^2 term
    channelGain_RIS(j) = sum(sqrt(betaHn.*betaGn)).^2/penetrationLoss;
    
    %Compute the  total channel gain with the RIS using Eq. (21) in [11]
    %when mimicking a mirror using theta_n=0 and mu_n=1 (by removing the
    %P/sigma^2 term)
    channelGain_RIS_approx(j) = abs(sum(sqrt(betaHn.*betaGn).*exp(-1i*(phaseGn)))).^2/penetrationLoss;
    
    %Compute the mirror limit in Eq. (54) in [11]
    channelGain_mirror(j) = lambda^2/((4*pi)^2*(deltaValues(j)+d)^2)/penetrationLoss;
    
end


%% Plot the simulation results
figure;
hold on; box on; grid on;
plot(deltaValues,pow2db(channelGain_RIS),'b-','LineWidth',2);
plot(deltaValues,pow2db(channelGain_RIS_approx),'r-.','LineWidth',2);
plot(deltaValues,pow2db(channelGain_mirror),'k--','LineWidth',2);
xlabel('Distance from surface [m]','Interpreter','Latex');
ylabel('Pathloss [dB]','Interpreter','Latex');
legend({'RIS (Optimal)','RIS (Mirror-mimicking)','Mirror approximation'},'Interpreter','Latex','Location','NorthEast');
set(gca,'fontsize',18);
