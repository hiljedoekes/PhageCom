%%% Plot results on short time scale (couple of transfers) %%%

clear all;

%%%%% PARAMETERS %%%%%

r = 1.0;                % Growth rate of bacteria (per hour)
K = 1.0;                % Scaled carrying capacity of bacteria
a = 10;                 % Scaled adsorption rate free phages to cells (per hour)
B = 2;                  % Effective burst size (per adsorption to susceptible cell)
alpha = 0.001;          % Reactivation rate of lysogens (per hour)
deltaP = 0.01;          % Decay rate of free phage particles (per hour)

mu = 0.001;             % Mutation rate (prob of mutation in strain j)
mut_type = "step";      % Mutation type ("step" for stepwise or "rand" for random)

cL = 0;                 % Scaled arbitrium production rate by lysogens (per hour)
deltaA = 0.0;           % Decay rate of free arbitrium (per hour)
u = 0.1;                % Scaled arbitrium uptake + break down rate by cells (per hour)

Ttransfer = 12;          % Period between serial transfers (hours)
Tmax = 18;              % Maximum number of hours "simulated"
ntsteps = 1e4;          % (approximate) Number of time steps plotted

Init_phage_dens = 1e-5;      % Initial total phage density
dilute = 0.01;               % Dilution factor of supernatant (phages) at each transfer

phimax = 1;             % Vector of phimax values of different strains
thres = 0.3:0.3:0.9;    % Vector of threshold values of different strains

nphi = length(phimax);  % Number of different phi values.
nthres = length(thres); % Number of different thres values.
ns = nphi*nthres;

InitS = K;                % Initial density of S
InitL = 0.0*ones(1,ns);     % Initial density of lysogens
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
if mod(Tmax, Ttransfer) ~= 0
    nt = floor(Tmax / Ttransfer);       % Number of transfers
else
    nt = Tmax/ Ttransfer - 1;
end
dt_tot = Tmax/ntsteps;
steps_trans = floor(Ttransfer / dt_tot) + 1;    % Number of steps per transfer
dt = Ttransfer / steps_trans;
tau = 0:dt:Ttransfer;               % Time vector for single transfers
ltau1 = length(tau) - 1;                         
t = 0:dt:Tmax;                      % Time vector for total time
lt = length(t);


%%%%% MODEL %%%%%

% Matrices to store full dynamics over time:
S = zeros(1,lt);
A = zeros(1,lt);
L = zeros(ns,lt);
P = zeros(ns,lt);

% Define intial conditions: Susceptibles, Arbitrium, Lysogens, Phages
Init = [InitS,0,InitL,InitP];

% Solve ODE system with ode-solver to simulate single transfer, store results, and restart after transfer
for i = 0:nt
    [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phimax_all,thres_all,nphi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
    S(1,(1+i*ltau1):(ltau1 + i*ltau1)) = S_temp(1,1:ltau1);
    A(1,(1+i*ltau1):(ltau1 + i*ltau1)) = A_temp(1,1:ltau1);
    L(:,(1+i*ltau1):(ltau1 + i*ltau1)) = L_temp(:,1:ltau1);
    P(:,(1+i*ltau1):(ltau1 + i*ltau1)) = P_temp(:,1:ltau1);
    % New initial condition
    Init = [InitS,0,zeros(1,ns),dilute*(P_temp(:,end))'];
end

% Shorten vectors if Tmax is not a multiple of Ttransfer
if mod(Tmax, Ttransfer) ~= 0
    S = S(1:lt);
    A = A(1:lt);
    L = L(:,1:lt);
    P = P(:,1:lt);
end

%%%%% PLOTS FOR EXPLORATION %%%%%

%%% Absolute numbers %%%
figure;
set(gcf,'DefaultAxesColorOrder',parula(ns+1))   % colour scheme
timeticks = 0:floor(Tmax/10):Tmax;  % Location of time tickmarks

% Bacteria
subplot(3,1,1)
plot(t,S,'LineWidth',2,'Color',[0 0 0]);
hold on;
ax = gca;
ax.ColorOrderIndex = 1;     % Reset colour scheme
plot(t,L(:,1:length(t)),'LineWidth',2);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Bacteria','Fontsize',14);
%Legend
legendentries = cell(ns+1,1);
legendentries(1,1) = {'Susceptibles'};
for j = 2:(ns+1)
    legendentries(j,1) = { ['phi =',num2str(phimax_all(j-1)),', thres = ',num2str(thres_all(j-1))] };
end
leg = legend(legendentries,'Location','eastoutside');

% Arbitrium
subplot(3,1,2)
plot(t,A,'LineWidth',2,'Color',[0 0 0]);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Arbitrium','Fontsize',14);
leg = legend('phi = X, thres = .X','Location','eastoutside'); % To make plot as wide as other panels

% Free phages
subplot(3,1,3)
plot(t,P,'LineWidth',2);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Free phages','Fontsize',14);
% Legend
leg = legend(legendentries(2:end,:),'Location','eastoutside');


%%% Relative frequencies %%%

figure;
set(gcf,'DefaultAxesColorOrder',parula(ns+1))
timeticks = 0:floor(Tmax/10):Tmax;

% Lysogens
subplot(2,1,1)
plot(t,L./sum(L,1),'LineWidth',2);
ylim([0 1]);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Bacteria','Fontsize',14);
%Legend
legendentries = cell(ns+1,1);
legendentries(1,1) = {'Susceptibles'};
for j = 2:(ns+1)
    legendentries(j,1) = { ['phi =',num2str(phimax_all(j-1)),', thres = ',num2str(thres_all(j-1))] };
end
leg = legend(legendentries(2:end,:),'Location','eastoutside');

% Free phages
subplot(2,1,2)
plot(t,P./sum(P,1),'LineWidth',2);
ylim([0 1]);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Free phages','Fontsize',14);
% Legend
leg = legend(legendentries(2:end,:),'Location','eastoutside');


%%%%% SEPARATE PLOTS FOR PUBLICATION %%%%%
%{
timeticks = 0:6:Tmax;       % Location of tickmarks on x-axis.

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