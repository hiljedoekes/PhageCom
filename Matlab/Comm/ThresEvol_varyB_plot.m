%%% Analyse 20191206_ThresEvol runs

clear all;

% Characteristics of phages
thres = 0:0.02:1;      % Vector of threshold values of different strains
ns = length(thres);    % Number of strains

%%%%% RUN 1 - beta varied

% Parameters
logB = 0.03:0.03:3;
B = 10.^(logB);
nsamples = length(B);

% Simulation results
p = dlmread('~/PhageCom/Results/20200214_ThresEvol/p.txt');
S = dlmread('~/PhageCom/Results/20200214_ThresEvol/S.txt');

% Check if all runs solved correctly : S should not be -1
min(S)


% Determine most prevalent strain
thressel = zeros(1,nsamples);
thressel_u0 = zeros(1,nsamples);
[max_p, max_p_ind] = max(p,[],1);
for j=1:nsamples
    thressel(j) = thres(max_p_ind(j));
end

% Predicted thres-value
threspred = ones(1,nsamples) ./ (2 - (ones(1,nsamples)./B));

% `Full' prediction - uses the fact that some susceptible cells are left
% over
%{
threspred_full = zeros(1,nsamples);
for i=1:nsamples
    fun = @(y) ( ((1-1/x(i))*y)/(1-y) - 1 + exp(-beta(i)*(((1-1/x(i))/(1-y)))) );
    threspred_full(i) = fzero(fun, 0.6);
end
%}


% Plot results - plot both thressel and threspred as functions of B
figure;
semilogx(B,thressel,'-o','LineWidth',2,'MarkerSize',3);
%hold on
%semilogx(x,thressel_u0,'-o','LineWidth',2);
hold on
semilogx(B,threspred,'LineWidth',2);
%legend({'Simulation, u=0.1','Simulation, u=0','Predicted'},'Location','eastoutside');
xlabel('B','Fontsize',14);
ylabel('Threshold','Fontsize',14);
ylim([0 1]);
xlim([1 1000]);


