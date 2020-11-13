%%% Plot short-term results of competition where the selected communicating
%%% strain invades on the selected bet-hedging strain

clear all;

%%%%% PARAMETERS %%%%%

r = 1.0;                % Growth rate of bacteria (per hour)
K = 1.0;                % Scaled carrying capacity of bacteria
a = 10;                 % Scaled adsorption rate free phages to cells (per hour)
B = 2;                  % Effective burst size (per adsorption to susceptible cell)
alpha = 0.001;          % Reactivation rate of lysogens (per hour)
deltaP = 0.01;          % Decay rate of free phage particles (per hour)

mu = 0.0;               % Competition experiment, so no mutations
mut_type = "step";      % Mutation type ("step" or "rand")

cL = 0;                 % Scaled arbitrium production rate by lysogens (per hour)
deltaA = 0.0;           % Decay rate of free arbitrium (per hour)
u = 0.1;                % Scaled arbitrium uptake + break down rate by cells (per hour)

Ttransfer = 12;          % Period between serial transfers (hours)
Tmax = 30;              % Maximum number of hours simulated
ntsteps = 1e4;          % (approximate) Number of time steps plotted

Init_phage_dens = 1e-5;      % Initial density of the bet-hedging phage
eps = 0.01;                  % Initial density of communicating invader compared to the bet-hedging resident
dilute = 0.01;               % Dilution factor of supernatant (phages) at each transfer

phimax_all = [0.04 1];      % phimax values of bet-hedging and communicating phages
thres_all = [0 0.66];       % thres values of bet-hedging (0) and communicating phages
ns = 2;                     % Number of strains
nphi = 2;                   % Number of different phi-values

InitS = 1.0;                % Initial density of S
InitL = 0.0*ones(1,ns);     % Initial density of lysogens
InitP = [Init_phage_dens eps*Init_phage_dens];    % Initial density of free phages

% Check length of initial vectors
if ( length(InitL) ~= ns || length(InitP) ~= ns)
    disp('Error: InitL and InitP should have length ns')
    return
end


%%%%% Time vectors %%%%%
if mod(Tmax, Ttransfer) ~= 0
    nt = floor(Tmax / Ttransfer);       % Number of transfers
else
    nt = Tmax/ Ttransfer - 1;
end
dt_tot = Tmax/ntsteps;
steps_trans = floor(Ttransfer / dt_tot) + 1;    % Number of steps per transfer
dt = Ttransfer / steps_trans;
tau = 0:dt:Ttransfer;               % Time vector for single transfer cycles
ltau1 = length(tau) - 1;                         
t = 0:dt:Tmax;                      % Total time vector
lt = length(t);


%%%%% MODEL %%%%%

% Matrices to store full dynamics over time
S = zeros(1,lt);
A = zeros(1,lt);
L = zeros(ns,lt);
P = zeros(ns,lt);

% Define intial conditions: Susceptibles, Arbitrium, Lysogens, Phages
Init = [InitS,0,InitL,InitP];

% Solve ODE system with ode-solver to simulate single transfer, store final time point, and restart after transfer
for i = 0:nt
    [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phimax_all,thres_all,nphi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
    S(1,(1+i*ltau1):(ltau1 + i*ltau1)) = S_temp(1,1:ltau1);
    A(1,(1+i*ltau1):(ltau1 + i*ltau1)) = A_temp(1,1:ltau1);
    L(:,(1+i*ltau1):(ltau1 + i*ltau1)) = L_temp(:,1:ltau1);
    P(:,(1+i*ltau1):(ltau1 + i*ltau1)) = P_temp(:,1:ltau1);
    % New initial condition
    Init = [InitS,0,zeros(1,ns),dilute*(P_temp(:,end))'];
end

%%%%% PLOT %%%%%

timeticks = 0:6:Tmax;       % location of time tickmarks

% Bacteria
fig = figure;
set(fig,'DefaultAxesColorOrder',parula((ns+1)))
plot(t,S(1:length(t)),'LineWidth',3, 'Color',[0,0,0]);
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,L(:,1:length(t)),'LineWidth',3);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Bacteria','Fontsize',14);
leg = legend('Susceptible','Bet-hedging','Communicating','Location','eastoutside');

%{
% Phages
fig = figure;
set(fig,'DefaultAxesColorOrder',parula((ns+1)))
plot(t,P(:,1:length(t)),'LineWidth',3);
set(gca, 'Xtick',timeticks,'Fontsize',12);
ylim([0, inf]);
xlabel('Time (hours)','Fontsize',14);
ylabel('Free Phages','Fontsize',14);

% Arbitrium
fig = figure;
plot(t,A(1:length(t)),'LineWidth',3, 'Color',[0,0,0]);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Arbitrium','Fontsize',14);
%}


