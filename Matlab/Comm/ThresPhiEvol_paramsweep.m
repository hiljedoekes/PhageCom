%%% Find distribution of phimax- and thres-values after many transfers for long time between transfers (Tt = 24)
%%% Consider many different parameter settings %%%
%%% Save results to file, this script takes a while (hours/days) to run %%%
%%% Script was run twice (with different random seed) in parallel %%%

clear all;
savefolder = '20200210_ThresPhiEvol_paramsweep1';

%%%%% PARAMETERS - not varied %%%%%

r = 1.0;                % Growth rate of bacteria (per hour)
K = 1.0;                % Scaled carrying capacity of bacteria

mu = 0.001;             % Mutation rate (prob of mutation in strain j)
mut_type = "step";      % Mutation type ("step" or "rand")

cL = 0;                 % Scaled arbitrium production rate by lysogens (per hour)
deltaA = 0.0;           % Decay rate of free arbitrium (per hour)

Tt = 24;                % Time between transfers
ntausteps = 1e2;        % Number of time steps returned in calculation of a single transfer
nt = 1000;               % Number of transfers simulated 

dtau = Tt/ntausteps;    % Time vector for single transfer cycle
tau = 0:dtau:Tt;
ltau = length(tau);

% Characteristics of phages
phimax = 0:0.1:1;
thres = 0:0.1:1;

nphi = length(phimax);  % Number of different phi values.
nthres = length(thres); % Number of different thres values.
ns = nphi*nthres;

InitS = 1.0;               % Initial density of S
InitL = 0.0*ones(1,ns);    % Initial density of lysogens
InitP = (1e-4/ns)*ones(1,ns);   % Initial density of phages

%%%%% Construct vectors that hold phi and thres values for all strains %%%%
phimax_all = zeros(1,ns);
thres_all = zeros(1,ns);
for k=1:ns
    phimax_ind = mod(k,nphi);
    if phimax_ind == 0
        phimax_all(k) = phimax(nphi);
    else
        phimax_all(k) = phimax(phimax_ind);
    end
    thres_ind = 1 + floor((k-1)/nphi);
    thres_all(k) = thres(thres_ind);
end

%%%%% PARAMETERS - varied %%%%%
nsamples = 250;         % number of runs (each with different parameter settings)
B = zeros(1,nsamples);  % Effective burst size
a = zeros(1,nsamples);  % Scaled adsorption rate phages to cells
alpha = zeros(1,nsamples);  % Lysogen reactivation rate
deltaP = zeros(1,nsamples); % Spontaneous degradation rate phages
u = zeros(1,nsamples);      % Scaled arbitrium uptake and degradation rate
dilution_factor = zeros(1,nsamples);    % Dilution factor of phages at transfer

% Seed of random number generator
rng(11);
% rng(12); % for paramsweep2

for i=1:nsamples
    B(i) = 10^(3*rand);     % log-uniform on 1 - 1000
    a(i) = 10^(2*rand);     % log-uniform on 1 - 100
    deltaP(i) = 10^(-(1+2*rand));   % log-uniform on 1e-3 - 0.1
    alpha(i) = 10^(-(2+2*rand));    % log-uniform on 1e-4 - 1e-2
    u(i) = 10^(-(3*rand));          % log-uniform on 1e-3 - 1
    dilution_factor(i) = 10^(-(1+2*rand));  % log-uniform on 0.001 - 0.1
end

% Save parameters
dlmwrite([savefolder,'/B.txt'],B,'precision',15);
dlmwrite([savefolder,'/alpha.txt'],alpha,'precision',15);
dlmwrite([savefolder,'/deltaP.txt'],deltaP,'precision',15);
dlmwrite([savefolder,'/a.txt'],a,'precision',15);
dlmwrite([savefolder,'/dilution_factor.txt'],dilution_factor,'precision',15);


%%%%% Matrices to store final values in %%%%%
S = zeros(1,nsamples);
L = zeros(ns,nsamples);
P = zeros(ns,nsamples);
l = zeros(ns,nsamples);      % Relative frequencies
p = zeros(ns,nsamples);      % Relative frequencies

%%%% For each parameter setting, find S, L and P after nt transfers
for i=1:nsamples
    disp(['Run ', num2str(i), ' of ', num2str(nsamples)]);
    Init = [InitS,0,InitL,InitP];
    exitflag = 0;
    % Solve ODE system with ode-solver to simulate single transfer, and restart after transfer
    for j = 1:nt
        [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phimax_all,thres_all,nphi,r,K,B(i),alpha(i),deltaP(i),a(i),cL,deltaA,u(i),mu,mut_type);
        if length(S_temp) ~= ltau
            exitflag = -1;
            break
        else
            % New initial condition
            Init = [InitS,0,zeros(1,ns),dilution_factor(i)*(P_temp(:,end))'];
        end
    end
    if exitflag == 0        % ODE-integration successful
        % Store final outcome
        S(1,i) = S_temp(1,ltau);
        L(:,i) = L_temp(:,ltau);
        P(:,i) = P_temp(:,ltau);
        % Store relative frequencies
        l(:,i) = L(:,i) / sum(L(:,i));
        p(:,i) = P(:,i) / sum(P(:,i));
    else
        S(1,i) = -1;
        L(:,i) = -1*ones(1,ns);
        l(:,i) = -1*ones(1,ns);
        P(:,i) = -1*ones(1,ns);
        p(:,i) = -1*ones(1,ns);
    end
end


%%%% Write output to file
dlmwrite([savefolder,'/S.txt'],S,'precision',15);
dlmwrite([savefolder,'/P.txt'],P,'precision',15);
dlmwrite([savefolder,'/p.txt'],p,'precision',15);
dlmwrite([savefolder,'/L.txt'],L,'precision',15);
dlmwrite([savefolder,'/l.txt'],l,'precision',15);


%%% Exit matlab
quit;

