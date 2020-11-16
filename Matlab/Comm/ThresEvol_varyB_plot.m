%%% Plot results of ThresEvol_varyB runs and predicted threshold values %%%

clear all;

% Simulation results
p = dlmread('../Data/20200214_ThresEvol_varyB/p.txt');

% Characteristics of phages included in the runs
thres = 0:0.02:1;      % Vector of threshold values of different strains
ns = length(thres);    % Number of strains

% Parameter value for the different runs
logB = 0.03:0.03:3;
B = 10.^(logB);
nsamples = length(B);

% Determine thres-value of most prevalent strain for each run
thressel = zeros(1,nsamples);
[max_p, max_p_ind] = max(p,[],1);
for j=1:nsamples
    thressel(j) = thres(max_p_ind(j));
end

% Calculate predicted thres-value for all values of B included
threspred = ones(1,nsamples) ./ (2 - (ones(1,nsamples)./B));


%%% Plot thressel and threspred as a function of B %%%
figure;
semilogx(B,thressel,'-o','LineWidth',2,'MarkerSize',3);
hold on
semilogx(B,threspred,'LineWidth',3,'Color',[0 0.7 0 0.5]);
legend({'Simulation results','Analytical prediction'},'Location','eastoutside');
xlabel('B','Fontsize',14);
ylabel('Threshold','Fontsize',14);
ylim([0 1]);
xlim([1 1000]);


