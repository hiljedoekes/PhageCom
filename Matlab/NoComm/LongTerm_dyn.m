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
mut_type = "step";      % Mutation type ("step" for stepwise, "rand" for random)

cL = 0;                 % Scaled arbitrium production rate by lysogens (per hour)
deltaA = 0.0;           % Decay rate of free arbitrium (per hour)
u = 0.1;                % Scaled arbitrium uptake + break down rate by cells (per hour)

Ttransfer = 12;         % Period between serial transfers (hours)
ntausteps = 1e2;        % Number of time steps returned in calculation of a single transfer
nt = 100;                % Number of transfers simulated
Tmax = nt*Ttransfer;    % Total time simulated

phi = 0:0.02:0.1;          % Vector of phi values of different strains
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
dtau = Ttransfer / ntausteps;
tau = 0:dtau:Ttransfer;                 % Within-transfer time vector
ltau = length(tau);

t = 0:Ttransfer:(nt*Ttransfer);         % Long-term time vector, one time point per transfer
lt = length(t);


%%%%% MODEL %%%%%
tic;

% Matrices to store full dynamics over time:
S = zeros(1,lt);        % Susceptibles
A = zeros(1,lt);        % Arbitrium
L = zeros(ns,lt);       % Lysogens (different strains)
P = zeros(ns,lt);       % Phages (different strains)

% Define intial conditions: Susceptible, Arbitrium, Lysogens, Phage
Init = [InitS,0,InitL,InitP];
S(1,1) = InitS;
A(1,1) = 0;
L(:,1) = InitL;
P(:,1) = InitP;

% Solve ODE system with ode-solver within transfers, save results, and restart after transfer
for i = 2:lt
    [S_temp,A_temp,L_temp,P_temp] = solve_ode(tau,Init,phi,r,K,B,alpha,deltaP,a,cL,deltaA,u,mu,mut_type);
    S(1,i) = S_temp(1,ltau);        % Save last time point within transfer
    A(1,i) = A_temp(1,ltau);
    L(:,i) = L_temp(:,ltau);
    P(:,i) = P_temp(:,ltau);
    % New initial condition
    Init = [InitS,0,zeros(1,ns),dilute*(P_temp(:,end))'];
end

simultime = toc;
disp(['Simulation run in: ',num2str(simultime,'%.2f'),' seconds'])

%%%%% PLOTS FOR EXPLORATION %%%%%

figure;
set(gcf,'DefaultAxesColorOrder',parula((ns)))       % colour scheme

timeticks = 0:floor(Tmax/10):Tmax;      % Location of time tickmarks

% Lysogens
subplot(2,1,1)
plot(t,L,'LineWidth',2);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Lysogens','Fontsize',14);
%Legend
legendentries = cell(ns,1);
for j = 1:ns
    legendentries(j,1) = { ['phi = ',num2str(phi(j))] };
end
leg = legend(legendentries,'Location','eastoutside');

% Free phages
subplot(2,1,2)
plot(t(2:end),P(:,2:end),'LineWidth',2);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Free phages','Fontsize',14);
% Legend
leg = legend(legendentries,'Location','eastoutside');


%%%%% SEPARATE PLOTS FOR PUBLICATION %%%%%
%{
timeticks = 0:(10*Ttransfer):Tmax;  % Location of time tickmarks

% Lysogens
fig = figure;
set(fig,'DefaultAxesColorOrder',parula((ns)))
plot(t,L,'LineWidth',3);
set(gca, 'Xtick',timeticks,'Fontsize',12);
xlabel('Time (hours)','Fontsize',14);
ylabel('Bacteria','Fontsize',14);

% Phages
fig = figure;
set(fig,'DefaultAxesColorOrder',parula((ns)))
plot(t,P,'LineWidth',3);
set(gca, 'Xtick',timeticks,'Fontsize',12);
ylim([0, inf]);
xlabel('Time (hours)','Fontsize',14);
ylabel('Free Phages','Fontsize',14);

%}

