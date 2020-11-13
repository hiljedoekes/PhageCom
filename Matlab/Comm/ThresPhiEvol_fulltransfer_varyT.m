%%% Find distribution of phimax- and thres-values after many transfers for varying values of between transfer time T %%%
%%% Save data to use for plotting later (script takes some time to run) %%%
%%% Transfer arbitrium, lysogens, and phages, instead of phages only %%%

clear all;
savefolder = '20200819_ThresPhiEvol_varyT';

%%%%% PARAMETERS %%%%%

r = 1.0;                % Growth rate of bacteria (per hour)
K = 1.0;                % Scaled carrying capacity of bacteria
a = 10;                 % Scaled adsorption rate free phages to cells (per hour)
B = 2;                  % Effective burst size (per adsorption to susceptible cell)
alpha = 0.001;          % Reactivation rate of lysogens (per hour)
deltaP = 0.01;          % Decay rate of free phage particles (per hour)

mu = 0.001;             % Mutation rate (prob of mutation in strain j)
mut_type = "step";      % Mutation type ("step" or "rand")

cL = 0;                 % Scaled arbitrium production rate by lysogens (per hour)
deltaA = 0.0;           % Decay rate of free arbitrium (per hour)
u = 0.1;                % Scaled arbitrium uptake + break down rate by cells (per hour)

ntausteps = 1e2;        % Number of time steps returned in calculation of a single transfer
nt = 1000;               % Number of transfers simulated

Ttransfer = [0.5:0.5:12,24];         % Vector of periods between serial transfers (hours)
nTt = length(Ttransfer);

Init_phage_dens = 1e-5;   % Initial total phage density

dilute = 0.01;          % Dilution factor at transfer


% Characteristics of phages
phimax = 0:0.05:1;           % Vector of phimax values of different strains
thres = 0:0.05:1;            % Vector of threshold values of different strains

nphi = length(phimax);  % Number of different phi values.
nthres = length(thres); % Number of different thres values.
ns = nphi*nthres;

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

% Initial conditions
InitS = 1.0;                % Initial density of S
InitL = 0.0*ones(1,ns);    % Initial density of lysogens
InitP = (Init_phage_dens/ns)*ones(1,ns);    % Initial density of free phages

if ( length(InitL) ~= ns || length(InitP) ~= ns)
    disp('Error: InitL and InitP should have length ns')
    return
end

%%%%% Matrices to store final values in %%%%%
L = zeros(ns,nTt);
P = zeros(ns,nTt);
l = zeros(ns,nTt);      % Relative frequencies
p = zeros(ns,nTt);      % Relative frequencies
solve_exitflag = zeros(1,nTt);  % indicates if ODE-solving was ok (no integration errors)

%%%% For each Ttransfer, find S, L and P after nt transfers %%%%
for i=1:nTt
    Tt = Ttransfer(i);
    disp(['Run 1, Tt = ', num2str(Tt)]);
    % Reset initial values
    Init = [InitS,0,InitL,InitP];
    % Define single transfer period time vector 
    dtau = Tt/ntausteps;
    tau = 0:dtau:Tt;
    ltau = length(tau);
    % Solve ODE system with ode-solver to simulate single transfer, and restart after transfer
    exitflag = 0;
    for j = 1:nt
        [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phimax_all,thres_all,nphi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
        if length(S_temp) ~= ltau
            exitflag = -1;
            break;
        else
            % New initial condition - transfer arbitrium, lysogens, and
            % phages
            Init = [InitS,dilute*A_temp(end),dilute*(L_temp(:,end))',dilute*(P_temp(:,end))'];
        end
    end
    solve_exitflag(i) = exitflag;
    % Store the final L and P distributions after nt transfers
    if exitflag == 0        % ODE solving successful
        L(:,i) = L_temp(:,ltau);
        P(:,i) = P_temp(:,ltau);
        % Store relative frequencies
        l(:,i) = L(:,i) / sum(L(:,i));
        p(:,i) = P(:,i) / sum(P(:,i));
    else
        L(:,i) = -1*ones(1,ns);
        P(:,i) = -1*ones(1,ns);
        l(:,i) = -1*ones(1,ns);
        p(:,i) = -1*ones(1,ns);
    end    

end

%%%% Write output to file
dlmwrite([savefolder,'/P.txt'],P(:,:),'precision',15);
dlmwrite([savefolder,'/p.txt'],p(:,:),'precision',15);
dlmwrite([savefolder,'/L.txt'],L(:,:),'precision',15);
dlmwrite([savefolder,'/l.txt'],l(:,:),'precision',15);
dlmwrite([savefolder,'/solve_exitflag.txt'],solve_exitflag);

% Exit matlab
exit;