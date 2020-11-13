%%% Plot results on long time scale, only plot last points of transfers %%%

clear all;

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

Ttransfer = 24;         % Period between serial transfers (hours)
ntausteps = 1e2;        % Number of time steps returned in calculation of a single transfer
nt = 100;               % Number of transfers simulated
Tmax = nt*Ttransfer;

phimax = 1;             % Vector of phimax values of different strains
thres = 0:0.2:1;        % Vector of threshold values of different strains

nphi = length(phimax);  % Number of different phi values.
nthres = length(thres); % Number of different thres values.
ns = nphi*nthres;

Init_phage_dens = 1e-5;    % Total initial phage density
dilute = 0.01;             % Dilution factor of supernatant (phages) at transfer

InitS = 1.0;                % Initial density of S
InitL = 0.0*ones(1,ns);    % Initial density of lysogens
InitP = (Init_phage_dens/ns)*ones(1,ns);    % Initial density of free phages

% Check length of initial vectors
if ( length(InitL) ~= ns || length(InitP) ~= ns)
    disp('Error: InitL and InitP should have length ns')
    return
end

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

%%%%% Time vectors %%%%%
dtau = Ttransfer / ntausteps;
tau = 0:dtau:Ttransfer;                 % Within-transfer time vec
ltau = length(tau);

t = 0:Ttransfer:(nt*Ttransfer);         % Long-term time vec
lt = length(t);


%%%%% MODEL %%%%%
tic;
% Matrices to store full dynamics over time
S = zeros(1,lt);
A = zeros(1,lt);
L = zeros(ns,lt);
P = zeros(ns,lt);

% Define intial conditions: Susceptibles, Arbitrium, Lysogens, Phages
Init = [InitS,0,InitL,InitP];
S(1,1) = InitS;
A(1,1) = 0;
L(:,1) = InitL;
P(:,1) = InitP;

% Solve ODE system with ode-solver to simulate single transfer, store final time point, and restart after transfer
for i = 2:lt
    [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phimax_all,thres_all,nphi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
    S(1,i) = S_temp(1,ltau);
    A(1,i) = A_temp(1,ltau);
    L(:,i) = L_temp(:,ltau);
    P(:,i) = P_temp(:,ltau);
    % New initial condition
    Init = [InitS,0,zeros(1,ns),dilute*(P_temp(:,end))'];
end

simultime = toc;
disp(['Simulation run in: ',num2str(simultime,'%.2f'),' seconds'])

%%%%% PLOT %%%%%
figure;
set(gcf,'DefaultAxesColorOrder',parula(ns)) % colour scheme
timeticks = 0:floor(Tmax/10):Tmax;      % location of time tickmarks

% Bacteria
subplot(2,1,1)
plot(t,S,'LineWidth',2,'Color',[0 0 0]);
hold on;
ax = gca;
ax.ColorOrderIndex = 1;     % reset colour numbering
plot(t,L,'LineWidth',2);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Bacteria','Fontsize',14);
%Legend
legendentries = cell(ns+1,1);
legendentries(1,1) = {'Susceptibles'};
for j = 2:(ns+1)
    legendentries(j,1) = { ['phi = ',num2str(phimax_all(j-1)), ' thres = ',num2str(thres_all(j-1))] };
end
leg = legend(legendentries,'Location','eastoutside');

% Free phages
subplot(2,1,2)
plot(t,P,'LineWidth',2);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Free phages','Fontsize',14);
% Legend
leg = legend(legendentries(2:end,:),'Location','eastoutside');

