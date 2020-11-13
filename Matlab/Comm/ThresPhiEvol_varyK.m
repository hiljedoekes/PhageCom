%%% Find distribution of phimax- and thres-values after many transfers for varying values of between transfer time T %%%
%%% Save data to use for plotting later (script takes some time to run) %%%
%%% Vary the carrying capacity from transfer cycle to transfer cycle %%%

clear all;
savefolder = '20200819_ThresPhiEvol_varyK';

%%%%% PARAMETERS %%%%%

r = 1.0;                % Growth rate of bacteria (per hour)
Kmean = 1.0;            % Mean of scaled carrying capacity of bacteria
Kvar = [0.025:0.025:0.2,0.25:0.05:1];   % Vector of variances of K
nruns = length(Kvar);

a = 10;                 % Scalde adsorption rate free phages to cells (per hour)
B = 2;                  % Effective burst size (per adsorption to susceptible cell)
alpha = 0.001;          % Reactivation rate of lysogens (per hour)
deltaP = 0.01;          % Decay rate of free phage particles (per hour)

mu = 0.001;             % Mutation rate (prob of mutation in strain j)
mut_type = "step";      % Mutation type ("step" or "rand")

cL = 0;                 % Scaled arbitrium production rate by lysogens (per hour)
deltaA = 0.0;           % Decay rate of free arbitrium (per hour)
u = 0.1;                % Scaled arbitrium uptake + break down rate by cells (per hour)

ntausteps = 1e2;        % Number of time steps returned in calculation of a single transfer
nt = 1000;              % Number of transfers simulated

Tt = 24;                % Consider a relatively long term between transfers

Init_phage_dens = 1e-5; % Initial total phage density

dilute = 0.01;          % Dilution factor at transfer

% Characteristics of phages
phimax = 0:0.05:1;           % Vector of phimax values of different strains
thres = 0:0.05:1;            % Vector of threshold values of different strains

nphi = length(phimax);  % Number of different phi values.
nthres = length(thres); % Number of different thres values.
ns = nphi*nthres;

%%%%% Construct phi and thres vector for all strains %%%%
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

% Initial conditions
InitS = Kmean;                % Initial density of S
InitL = 0.0*ones(1,ns);    % Initial density of lysogens
InitP = (Init_phage_dens/ns)*ones(1,ns);    % Initial density of free phages

if ( length(InitL) ~= ns || length(InitP) ~= ns)
    disp('Error: InitL and InitP should have length ns')
    return
end

% Define time vector for single transfers
dtau = Tt/ntausteps;
tau = 0:dtau:Tt;
ltau = length(tau);

%%%%% Matrices to store final values in %%%%%
L = zeros(ns,nt,nruns);
P = zeros(ns,nt,nruns);
Kvec = zeros(nt,nruns);
solve_exitflag = zeros(1,nruns);  % indicates if ODE-solving was ok (no integration errors)

%%%% For each Kvar, find S, L and P after nt transfers %%%%
for i=1:nruns
    Kv = Kvar(i);
    disp(['Variance of K = ', num2str(Kv)]);
    Init = [InitS,0,InitL,InitP];
    K = Kmean;      % Start with the mean carrying capacity
    % Solve ODE system with for nt transfers
    exitflag = 0;
    for j = 1:nt
        [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phimax_all,thres_all,nphi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
        if length(S_temp) ~= ltau
            exitflag = -1;
            break;
        else
            if j > (nt-100)     % Store K, and P and L distributions for final 100 transfers
                Kvec(j,i) = K;
                L(:,j,i) = L_temp(:,end);
                P(:,j,i) = P_temp(:,end);
            end
            % Pick a random carrying capacity for the next transfer cycle
            K = gamrnd((Kmean^2/Kv), (Kv/Kmean));
            % New initial condition
            Init = [K,0,zeros(1,ns),dilute*(P_temp(:,end))'];
        end
    end
    solve_exitflag(i) = exitflag;   

end

%%%% Write output to file
dlmwrite([savefolder,'/P.txt'],P(:,:,:),'precision',15);
dlmwrite([savefolder,'/L.txt'],L(:,:,:),'precision',15);
dlmwrite([savefolder,'/K.txt'],Kvec(:,:),'precision',15);
dlmwrite([savefolder,'/solve_exitflag.txt'],solve_exitflag);

% Exit matlab
exit;
