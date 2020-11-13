%%% Find distribution of phi-values after many transfers for varying values of between transfer time T %%%
%%% Save data to use for plotting later (script takes some time to run) %%%

clear all;
savefolder = '20191101_PhiEvol_varyT';

%%%%% PARAMETERS %%%%%

r = 1.0;                % Growth rate of bacteria (per hour)
K = 1.0;                % Scaled carrying capacity of bacteria
a = 10;                 % Scaled adsorption rate free phages to cells (per hour)
B = 2;                  % Effective burst size (per adsorption to susceptible cell)
alpha = 0.001;          % Reactivation rate of lysogens (per hour)
deltaP = 0.01;          % Decay rate of free phage particles (per hour)

mu = 0.001;              % Mutation rate (prob of mutation in strain j)
mut_type = "step";       % Mutation type ("step" or "rand")

cL = 0;                  % Scaled arbitrium production rate by lysogens (per hour)
deltaA = 0.0;           % Decay rate of free arbitrium (per hour)
u = 0.1;                % Scaled arbitrium uptake + break down rate by cells (per hour)

ntausteps = 1e2;        % Number of time steps returned in calculation of a single transfer
nt = 100;               % Number of transfers simulated before equilibrium finder is used
ntmax = 1e4;            % Number of transfers maximally performed per Ttransfer (if equilibrium finder doesn't converge)


% Parameters varied
Ttransfer = [0.5:0.5:12, 24];         % Vector of times between serial transfers (hours)
nTt = length(Ttransfer);

phi = 0:0.005:0.2;      % Vector of phi values of different strains
ns = length(phi);       % Number of strains

% Initial conditions
Init_phage_dens = 1e-5;      % Initial total phage density
dilute = 0.01;               % Dilution factor of supernatant (phages) after each transfer

InitS = 1.0;                % Initial density of S
InitL = 0.0*ones(1,ns);    % Initial vector of lysogen density
InitP = (Init_phage_dens/ns)*ones(1,ns);    % Initial vector of free phage density

% Check if InitL and InitP are correct length
if ( length(InitL) ~= ns || length(InitP) ~= ns)
    disp('Error: InitL and InitP should have length ns')
    return
end

%%%%% Matrices to store final values in %%%%%
S = zeros(1,nTt);
L = zeros(ns,nTt);
P = zeros(ns,nTt);
l = zeros(ns,nTt);      % Relative frequencies
p = zeros(ns,nTt);      % Relative frequencies
solve_exitflag = zeros(1,nTt);  % indicates if equilibrium could be solved with lsqnonlin

%%%% For each Ttransfer, calculate stable S, L and P densities at end of transfer %%%%
for i=1:nTt
    Tt = Ttransfer(i);
    disp(['Tt = ', num2str(Tt)]);
    Init = [InitS,0,InitL,InitP];   % Initial values vector
    % Define time vectors
    dtau = Tt/ntausteps;        
    tau = 0:dtau:Tt;            % Time vector for single transfer cycle
    ltau = length(tau);
    
    % First, solve the ODE system for nt transfers
    for j = 1:nt
        [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
        % New initial condition
        Init = [InitS,0,zeros(1,ns),dilute*(P_temp(:,end))'];
    end
    
    % Then use S, L and P after nt transfers as input to solve for the
    % 'transfer' equilibrium. ('transfer equilibrium' = stable distribution
    % given that the system is exposed to regular transfers)
    
    % Define a function that looks at the difference between input and
    % output-vector over one transfer cycle.
    z0 = [S_temp(end), A_temp(end), L_temp(:,end)',P_temp(:,end)'];
    fun = @(z) (one_transfer_results(tau,[InitS,0,zeros(1,ns),dilute*z((ns+3):end)],phi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type) - z);
    % Find root of this function using non-linear least squares
    lb = zeros(1,length(z0));
    options = optimoptions('lsqnonlin','Display','off');
    [z,resnorm,residual,exitflag,output] = lsqnonlin(fun,z0,lb,[],options);
    solve_exitflag(i) = exitflag;
    % Check if non-linear least squares was successful
    if exitflag > 0
        % Convergence successful, save final outcome
        S(1,i) = z(1);
        L(:,i) = z(3:(ns+2));
        P(:,i) = z((ns+3):end);
    else
        % Convergence not successful, find equilibrium by simulating many
        % more transfers
        for j = (nt+1):ntmax
            [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
            % New initial condition
            Init = [InitS,0,zeros(1,ns),dilute*(P_temp(:,end))'];          
        end
        % Store final outcome
        S(1,i) = S_temp(1,ltau);
        L(:,i) = L_temp(:,ltau);
        P(:,i) = P_temp(:,ltau);
    end    
    % Store relative frequencies
    l(:,i) = L(:,i) / sum(L(:,i));
    p(:,i) = P(:,i) / sum(P(:,i));
end

%%%% Write output to file
dlmwrite([savefolder,'/P.txt'],P,'precision',15);
dlmwrite([savefolder,'/p.txt'],p,'precision',15);
dlmwrite([savefolder,'/L.txt'],L,'precision',15);
dlmwrite([savefolder,'/l.txt'],l,'precision',15);
dlmwrite([savefolder,'/solve_exitflag.txt'],solve_exitflag);


