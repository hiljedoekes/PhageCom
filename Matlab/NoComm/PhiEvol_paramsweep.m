%%% Simulate evolution of phi for long time between transfers (Tt = 24) %%%
%%% Consider many different parameter settings %%%
%%% Save results to file (script takes hours - days to run) %%%

clear all;
savefolder = '20200210_PhiEvol_paramsweep';

%%%%% PARAMETERS - not varied %%%%%

r = 1.0;                % Growth rate of bacteria (per hour)
K = 1.0;                % Scaled carrying capacity of bacteria

mu = 0.001;             % Mutation rate (prob of mutation in strain j)
mut_type = "step";      % Mutation type ("step" or "rand")

cL = 0;                 % Scaled arbitrium production rate by lysogens (per hour)
deltaA = 0.0;           % Decay rate of free arbitrium (per hour)
u = 0.1;                % Scaled arbitrium uptake + break down rate by cells (per hour)
                        % Not varied because not relevant in absence of communication

Tt = 24;                % Time between transfers
ntausteps = 1e2;        % Number of time steps returned in calculation of a single transfer
nt = 100;               % Number of transfers simulated before equilibrium finder is used
ntmax = 1e4;            % Number of transfers maximally performed per Ttransfer (if equilibrium finder doesn't converge)

dtau = Tt/ntausteps;    % Time vectors
tau = 0:dtau:Tt;
ltau = length(tau);

phi = 0:0.01:1;         % Vector of phi values of different strains
ns = length(phi);       % Number of strains

InitS = 1.0;               % Initial density of S
InitL = 0.0*ones(1,ns);    % Initial density of lysogens
InitP = (1e-4/ns)*ones(1,ns);  % Initial density of phages

%%%%% PARAMETERS - varied %%%%%

% Vectors to store parameter values in
nsamples = 500;
B = zeros(1,nsamples);
a = zeros(1,nsamples);
alpha = zeros(1,nsamples);
deltaP = zeros(1,nsamples);
dilution_factor = zeros(1,nsamples);

% Seed of random number generator (vary if running sets in parallel)
rng(11);

% Randomly sample parameter values
for i=1:nsamples
    B(i) = 10^(3*rand);     % log-uniform on 1 - 1000
    a(i) = 10^(2*rand);     % log-uniform on 1 - 100
    deltaP(i) = 10^(-(1+2*rand));   % log-uniform on 1e-3 - 0.1
    alpha(i) = 10^(-(2+2*rand));    % log-uniform on 1e-4 - 1e-2
    dilution_factor(i) = 10^(-(1+2*rand));  % log-uniform on 1e-3 - 0.1
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
phisel = zeros(1,nsamples);  % Phi-value selected = most-prevalent strain
solve_exitflag = zeros(1,nsamples);

%%%% For each parameter setting, determine stable S, L and P distribution
%%%% given that transfers occur
for i=1:nsamples
    disp(['Run ', num2str(i), ' of ', num2str(nsamples)]);
    Init = [InitS,0,InitL,InitP];
    
    % First, solve ODE system with for nt transfers
    for j = 1:nt
        [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phi,r,K,B(i),alpha(i),deltaP(i),a(i),cL,deltaA,u,mu,mut_type);
        % New initial condition
        Init = [InitS,0,zeros(1,ns),dilution_factor(i)*(P_temp(:,end))'];
    end
    
    % Then, use S, L and P after nt transfers as input to solve for the
    % 'transfer equilibrium' ('transfer equilibrium' = stable distribution
    % given that the system is exposed to regular transfers).
    
    % Define a function that looks at the difference between input and
    % output-vector over one transfer cycle.
    z0 = [S_temp(end), A_temp(end), L_temp(:,end)',P_temp(:,end)'];
    fun = @(z) (one_transfer_results(tau,[InitS,0,zeros(1,ns),dilution_factor(i)*z((ns+3):end)],phi,r,K,B(i),alpha(i),deltaP(i),a(i),cL,deltaA,u,mu,mut_type) - z);
    % Find root of this function using non-linear least squares
    lb = zeros(1,length(z0));
    options = optimoptions('lsqnonlin','Display','off');
    [z,resnorm,residual,exitflag,output] = lsqnonlin(fun,z0,lb,[],options);
    solve_exitflag(i) = exitflag;
    % Check if non-linear least squares converged
    if exitflag > 0
        % Convergence successful, store final outcome
        S(1,i) = z(1);
        L(:,i) = z(3:(ns+2));
        P(:,i) = z((ns+3):end);
    else
        % Convergence not successful, find equilibrium by simulating many
        % more transfers
        for j = (nt+1):ntmax
            [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phi,r,K,B(i),alpha(i),deltaP(i),a(i),cL,deltaA,u,mu,mut_type);
            % New initial condition
            Init = [InitS,0,zeros(1,ns),dilution_factor(i)*(P_temp(:,end))'];          
        end
        % Store final outcome
        S(1,i) = S_temp(1,ltau);
        L(:,i) = L_temp(:,ltau);
        P(:,i) = P_temp(:,ltau);
    end    
    % Store relative frequencies
    l(:,i) = L(:,i) / sum(L(:,i));
    p(:,i) = P(:,i) / sum(P(:,i));
    % Store most prevalent phi value
    [max_p, max_p_ind] = max(p(:,i));
    modus_phi = (max_p_ind - 1)/100;          % Map back to phi -- SHOULD BE GENERALIZED if different number of strains is included
    phisel(i) = modus_phi;
end


%%%% Write output to file
dlmwrite([savefolder,'/S.txt'],S,'precision',15);
dlmwrite([savefolder,'/P.txt'],P,'precision',15);
%dlmwrite([savefolder,'/p.txt'],p,'precision',15);
dlmwrite([savefolder,'/L.txt'],L,'precision',15);
%dlmwrite([savefolder,'/l.txt'],l,'precision',15);
dlmwrite([savefolder,'/phisel.txt'],phisel,'precision',15);
dlmwrite([savefolder,'/solve_exitflag.txt'],solve_exitflag);

%%% Exit matlab
exit;
