%%% Look at evolution of thres for long time between transfers (Tt = 24) %%%
%%% Set phimax = 1, consider evolution of thres only %%%
%%% Save data, because script takes some time to run

clear all;
savefolder = '20200214_ThresEvol_varyB';

%%%%% PARAMETERS - default %%%%%

r = 1.0;                % Growth rate of bacteria (per hour)
K = 1.0;                % Scaled carrying capacity of bacteria
alpha = 0.001;          % Reactivation rate of lysogens (per hour)
deltaP = 0.01;          % Decay rate of free phage particles (per hour)
a = 10;                 % Scaled adsorption rate free phages to cells (per hour)

mu = 0.001;              % Mutation rate (prob of mutation in strain j)
mut_type = "step";      % Mutation type ("step" or "rand")

cL = 0;                 % Scaled arbitrium production rate by lysogens (per hour)
deltaA = 0.0;           % Decay rate of free arbitrium (per hour)
u = 0.1;                % Scaled arbitrium uptake + break down rate by cells (per hour)

dilution_factor = 0.01;

Tt = 24;                % Time between transfers
ntausteps = 1e2;        % Number of time steps returned in calculation of a single transfer
nt = 1000;               % Number of transfers simulated

dtau = Tt/ntausteps;    % Time vector for single transfer cycle
tau = 0:dtau:Tt;
ltau = length(tau);

% Characteristics of phages
thres = 0:0.02:1;      % Vector of threshold values of different strains
ns = length(thres);    % Number of strains
phimax = ones(1,ns);        % Vector of phimax values of different strains
nphi = 1;

InitS = 1.0;               % Initial density of S
InitL = 0.0*ones(1,ns);    % Initial density of lysogens
InitP = (1e-4/ns)*ones(1,ns);   % Initial density of phages

%%%%% Parameter varied %%%%%
logB = 0.03:0.03:3;         % logB = 0 cannot be included because phages with B=1 are not viable
B = 10.^(logB);
nsamples = length(B);


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
        [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phimax,thres,nphi,r,K,B(i),alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
        if length(S_temp) ~= ltau   % try-catch for integration errors
            exitflag = -1;
            break
        end
        % New initial condition
        Init = [InitS,0,zeros(1,ns),dilution_factor*(P_temp(:,end))'];
    end
    if exitflag == -1
        S(1,i) = j;
        L(:,i) = -1*ones(ns,1);
        P(:,i) = -1*ones(ns,1);
        l(:,i) = -1*ones(ns,1);
        p(:,i) = -1*ones(ns,1);
    else    % ODE-integration successful
        % Store final outcome
        S(1,i) = S_temp(1,ltau);
        L(:,i) = L_temp(:,ltau);
        P(:,i) = P_temp(:,ltau);
        % Store relative frequencies
        l(:,i) = L(:,i) / sum(L(:,i));
        p(:,i) = P(:,i) / sum(P(:,i));
    end
end

%%%% Write output to file
dlmwrite([savefolder,'/S.txt'],S,'precision',15);
dlmwrite([savefolder,'/P.txt'],P,'precision',15);
dlmwrite([savefolder,'/p.txt'],p,'precision',15);
dlmwrite([savefolder,'/L.txt'],L,'precision',15);
dlmwrite([savefolder,'/l.txt'],l,'precision',15);


