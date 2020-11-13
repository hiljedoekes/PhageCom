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
mut_type = "step";      % Mutation type ("step" for stepwise, "rand" for random)

cL = 0;                 % Scaled arbitrium production rate by lysogens (per hour)
deltaA = 0.0;           % Decay rate of free arbitrium (per hour)
u = 0.1;                % Scaled arbitrium uptake + break down rate by cells (per hour)

Ttransfer = 12;         % Period between serial transfers (hours)
Tmax = 36;              % Maximum number of hours simulated
ntsteps = 1e4;          % (approximate) Number of time steps plotted

phi = 0:0.02:0.1;       % Vector of phi values of different strains
ns = length(phi);       % Number of strains

Init_phage_dens = 1e-5; % Initial total phage density   
dilute = 0.01;          % Dilution factor of phages at transfer

InitS = 1.0;            % Initial density of susceptible cells
InitL = 0.0*ones(1,ns); % Initial densities of lysogens (vector)
InitP = (Init_phage_dens/ns)*ones(1,ns);    % Initial densities of free phages (vector)

% Check if initial vectors are correct length
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
dt_tot = Tmax/ntsteps;                          % Size of time steps
steps_trans = floor(Ttransfer / dt_tot) + 1;    % Number of steps per transfer
dt = Ttransfer / steps_trans;
tau = 0:dt:Ttransfer;                   % Time vector for single transfer cycle
ltau1 = length(tau) - 1;                         
t = 0:dt:Tmax;                          % Total time vector
lt = length(t);


%%%%% MODEL %%%%%

% Matrices to store full dynamics over time:
S = zeros(1,lt);        % Susceptible cells
A = zeros(1,lt);        % Arbitrium
L = zeros(ns,lt);       % Lysogens (different strains)
P = zeros(ns,lt);       % Phages (different strains)

% Define intial conditions: Susceptibles, Arbitrium, Lysogens, Phages
Init = [InitS,0,InitL,InitP];

% Solve ODE system with ode-solver within transfers, save results, and restart after transfer
for i = 0:nt
    [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
    S(1,(1+i*ltau1):(ltau1 + i*ltau1)) = S_temp(1,1:ltau1);
    A(1,(1+i*ltau1):(ltau1 + i*ltau1)) = A_temp(1,1:ltau1);
    L(:,(1+i*ltau1):(ltau1 + i*ltau1)) = L_temp(:,1:ltau1);
    P(:,(1+i*ltau1):(ltau1 + i*ltau1)) = P_temp(:,1:ltau1);
    % New initial condition
    Init = [InitS,0,zeros(1,ns),dilute*(P_temp(:,end))'];
end

%%%%%% PLOTS FOR EXPLORATION %%%%%%

timeticks = 0:6:Tmax;       % Location of time tickmarks

%%% Absolute numbers %%%
figure;
set(gcf,'DefaultAxesColorOrder',parula((ns)))   % colour scheme
% Bacteria
subplot(2,1,1)
plot(t,S,'LineWidth',3, 'Color',[0,0,0]);       % Susceptibles
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,L,'LineWidth',3);                        % Lysogens
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Bacteria','Fontsize',14);
%Legend
legendentries = cell(ns+1,1);
legendentries(1,1) = {'Susceptibles'};
for j = 2:(ns+1)
    legendentries(j,1) = { ['phi = ',num2str(phi(j-1))] };
end
leg = legend(legendentries,'Location','eastoutside');

% Phages
subplot(2,1,2)
plot(t,P,'LineWidth',3);
set(gca, 'Xtick',timeticks,'Fontsize',12);
ylim([0, inf]);
xlabel('Time (hours)','Fontsize',14);
ylabel('Free Phages','Fontsize',14);
% Legend
leg = legend(legendentries(2:end,:),'Location','eastoutside');


%%% Relative frequency %%%
figure;
set(gcf,'DefaultAxesColorOrder',parula((ns)))

% Lysogens
subplot(2,1,1)
plot(t,L./sum(L,1),'LineWidth',2);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Relative strain freq lysogens','Fontsize',14);
% Legend
leg = legend(legendentries(2:end,:),'Location','eastoutside');

% Free phages
subplot(2,1,2)
plot(t,(P./sum(P,1)),'LineWidth',2);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Relative strain freq Free phages','Fontsize',14);
% Legend
leg = legend(legendentries(2:end,:),'Location','eastoutside');


%%%%% SEPARATE PLOTS FOR PUBLICATION %%%%%
%{
timeticks = 0:6:Tmax;       % Location of time tickmarks

% Bacteria
fig = figure;
fig.Renderer = 'Painters';
set(fig,'DefaultAxesColorOrder',parula((ns)))   % colour scheme
plot(t,S,'LineWidth',3, 'Color',[0,0,0]);       % Susceptibles
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,L,'LineWidth',3);                        % Lysogens
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Bacteria','Fontsize',14);

% Phages
fig = figure;
fig.Renderer = 'Painters';
set(fig,'DefaultAxesColorOrder',parula((ns)))
plot(t,P,'LineWidth',3);
set(gca, 'Xtick',timeticks,'Fontsize',12);
ylim([0, inf]);
xlabel('Time (hours)','Fontsize',14);
ylabel('Free Phages','Fontsize',14);

% Relative frequency of lysogens
fig = figure;
fig.Renderer = 'Painters';
set(fig,'DefaultAxesColorOrder',parula((ns)))
plot(t,L./sum(L,1),'LineWidth',3);
set(gca, 'Xtick',timeticks,'Fontsize',12);
ylim([0, inf]);
xlabel('Time (hours)','Fontsize',14);
ylabel('Relative frequency of lysogens','Fontsize',14);

% Relative frequency of phages
fig = figure;
fig.Renderer = 'Painters';
set(fig,'DefaultAxesColorOrder',parula((ns)))
plot(t,P./sum(P,1),'LineWidth',3);
set(gca, 'Xtick',timeticks,'Fontsize',12);
ylim([0, inf]);
xlabel('Time (hours)','Fontsize',14);
ylabel('Relative frequency of phages','Fontsize',14);
%}

