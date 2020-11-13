%%%% Analysis of Parameter sweep %%%%
clear all;

% Read data
S = dlmread('~/PhageCom/Results/20200210_ThresPhiEvol_paramsweep1/S.txt');
temp = dlmread('~/PhageCom/Results/20200210_ThresPhiEvol_paramsweep2/S.txt');
S = [S temp];
p = dlmread('~/PhageCom/Results/20200210_ThresPhiEvol_paramsweep1/p.txt');
temp = dlmread('~/PhageCom/Results/20200210_ThresPhiEvol_paramsweep2/p.txt');
p = [p temp];

nsample = length(S);

%%% Check if all ODE integrations were successful (S = -1 set if not)
min(S)

%%% Find the most prevalent strain, and corresponding phi and thres 
%%% + Find the mean value of phi and thres in simulations
phimax = 0:0.1:1;
thres = 0:0.1:1;
nphi = length(phimax);
nthres = length(thres);
ns = nphi*nthres;
% Construct phi and thres vector for all strains
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

phisel = zeros(1,nsample);
thressel = zeros(1,nsample);

% Store index of most prevalent phage
[max_p, max_p_ind] = max(p,[],1);
% Translate back to phi-value
indphi = mod(max_p_ind,nthres);
for j=1:nsample
    if indphi(j) == 0
        phisel(j) = phimax(nphi);
        indphi(j) = nphi;hpi
    else
        phisel(j) = phimax(indphi(j));
    end
end
% Translate back to thres-value
indthres = floor((max_p_ind-1)/nphi)+1;
for j=1:nsample
        thressel(j) = thres(indthres(j));
end

% Find mean value of phi and thres in simulations
phimean = phimax_all*p;
thresmean = thres_all*p;

%%% PLOTS %%%

% Data structure for heatmap and bubble plot
heatmapdata = zeros(nphi,nthres);
for j=1:nsample
    heatmapdata(indphi(j),indthres(j)) = heatmapdata(indphi(j),indthres(j)) + 1;
end

% Bubble plot instead of heatmap
[Xbins, Ybins] = meshgrid(0:0.1:1,0:0.1:1);
ind = heatmapdata(:) > 0;
figure;
scatter(Xbins(ind), Ybins(ind), 4*(heatmapdata(ind).^1.0),'filled')
grid(gca,'on')
xlim([0 1]);
ylim([0 1]);
xlabel('Response threshold (\theta)','Fontsize',14);
ylabel('Lysogeny propensity (\phi_{max})','Fontsize',14);

% Plot mean values against each other
figure;
scatter(thresmean,phimean,16,'filled','MarkerFaceAlpha',0.5);
xlim([0 1])
ylim([0 1])
xlabel('Mean Response threshold (\theta)','Fontsize',14);
ylabel('Mean Lysogeny propensity (\phi_{max})','Fontsize',14);
title('Results of 500 parameter sweep runs','Fontsize',14);

% How many runs yielded phisel = 1?
length(phisel(phisel==1))
% What is phisel for the exceptions?
phisel(phisel<1)

% Histogram of thres values for runs with phisel = 1
figure;
histogram(thressel(phisel==1));
xlabel('Selected thres', 'Fontsize',14);
ylabel(['Number of runs (out of ',num2str(length(thressel(phisel==1))), ')'], 'Fontsize',14);
title('phisel = 1','Fontsize',16);


%%%% PART 2: ThresEvol for those runs that yielded phisel = 1 %%%%
%{
% Read data
S_thres = dlmread('~/PhageCom/Results/20191122_ThresEvol/S.txt');
p_thres = dlmread('~/PhageCom/Results/20191122_ThresEvol/p.txt');

nsample_t = length(S_thres(phisel==1));

% Check if all ODE integrations were successful (i.e. S = -1 doesn't occur)
min(S_thres)

% Only use the data of runs that yielded phi = 1 in ThresPhiEvol
p_thres = p_thres(:,phisel==1);

% Characteristics of phages
thres_t = 0:0.02:1;      % Vector of threshold values of different strains
ns_t = length(thres_t);    % Number of strains
% Find most abundant thres value
thressel_t = zeros(1,nsample_t);
[max_p_t, max_p_t_ind] = max(p_thres,[],1);
for j=1:nsample_t
    thressel_t(j) = thres_t(max_p_t_ind(j));
end
% Find mean thres value
thressel_t_av = thres_t*p_thres;


% Histogram of most abundant thres-value (ESS?) in simulations
figure;
histogram(thressel_t);
xlabel('Selected thres', 'Fontsize',14);
ylabel('Number of runs (out of 482)', 'Fontsize',14);

% Predicted thres-values
p0 = df.*((alpha.*(1-alpha))./(deltaP + a.*(1-alpha)));

% Letterlijk van mathematica:
full_threspred = -( ( beta.*( beta - a + beta.*p0 + p0.*(a-2*beta).*( (-( (p0.*(a-2*beta))./(beta-a+beta.*p0) )).^(u./(beta-a)) ) ) )./( (a-2*beta).*(beta-a+u) ) );

eta = 1 - a ./ beta;

simple_threspred_1 = eta ./ ((1+eta).*(eta + u./beta));
simple_threspred_2 = (1+eta).^(-1);

% Plot different predictions vs thressel
figure;
subplot(1,3,1);
scatter(full_threspred(phisel==1),thressel_t,15,'filled','MarkerFaceAlpha',0.5);
xlim([0.45 0.7]);
ylim([0.45 0.7]);
l = refline([1 0]);
l.Color = "black";
l.LineWidth = 2;
title('Full prediction, including P0','Fontsize',18);
xlabel('Predicted threshold','Fontsize',14);
ylabel('Most abundant threshold in simulation','Fontsize',14);
subplot(1,3,2);
scatter(simple_threspred_1(phisel==1),thressel_t,15,'filled','MarkerFaceAlpha',0.5);
xlim([0.45 0.7]);
ylim([0.45 0.7]);
l = refline([1 0]);
l.Color = "black";
l.LineWidth = 2;
title('Simplfied prediction, including u','Fontsize',18);
xlabel('Predicted threshold','Fontsize',14);
ylabel('Most abundant threshold in simulation','Fontsize',14);
subplot(1,3,3);
scatter(simple_threspred_2(phisel==1),thressel_t,15,'filled','MarkerFaceAlpha',0.5);
xlim([0.45 0.7]);
ylim([0.45 0.7]);
l = refline([1 0]);
l.Color = "black";
l.LineWidth = 2;
title('Simplified prediction, 1 / (1+\eta)','Fontsize',18);
xlabel('Predicted threshold','Fontsize',14);
ylabel('Most abundant threshold in simulation','Fontsize',14);

% Comparision between predictions
figure;
subplot(1,3,1);
scatter(simple_threspred_1(phisel==1),full_threspred(phisel==1),15,'filled','MarkerFaceAlpha',0.5);
xlim([0.45 1]);
ylim([0.45 1]);
l = refline([1 0]);
l.Color = "black";
l.LineWidth = 1;
xlabel('Simplified pred with u','Fontsize',14);
ylabel('Full prediction','Fontsize',14);

subplot(1,3,2);
scatter(simple_threspred_2(phisel==1),full_threspred(phisel==1),15,'filled','MarkerFaceAlpha',0.5);
xlim([0.45 1]);
ylim([0.45 1]);
l = refline([1 0]);
l.Color = "black";
l.LineWidth = 1;
xlabel('Simplified pred 1/(1+\eta)','Fontsize',14);
ylabel('Full prediction','Fontsize',14);

subplot(1,3,3);
scatter(simple_threspred_2(phisel==1),simple_threspred_1(phisel==1),15,'filled','MarkerFaceAlpha',0.5);
xlim([0.45 1]);
ylim([0.45 1]);
l = refline([1 0]);
l.Color = "black";
l.LineWidth = 1;
xlabel('Simplified pred (1/1+\eta)','Fontsize',14);
ylabel('Simplified pred with u','Fontsize',14);

% Does using measured P0 (instead of estimated/calculated P0) improve the
% prediction?
P_t = dlmread('~/PhageCom/Results/20191122_ThresEvol/P.txt');
P0_t = df.*sum(P_t,1);

% Letterlijk van mathematica:
P0_threspred = -( ( beta.*( beta - a + beta.*P0_t + P0_t.*(a-2*beta).*( (-( (P0_t.*(a-2*beta))./(beta-a+beta.*P0_t) )).^(u./(beta-a)) ) ) )./( (a-2*beta).*(beta-a+u) ) );

% Plot prediction vs thressel
figure;
subplot(1,2,1);
scatter(P0_threspred(phisel==1),thressel_t,15,'filled','MarkerFaceAlpha',0.5);
xlim([0.45 0.7]);
ylim([0.45 0.7]);
l = refline([1 0]);
l.Color = "black";
l.LineWidth = 2;
title('Full prediction, with measured P0','Fontsize',18);
xlabel('Predicted threshold','Fontsize',14);
ylabel('Most abundant threshold in simulation','Fontsize',14);
subplot(1,2,2);
scatter(simple_threspred_1(phisel==1),thressel_t,15,'filled','MarkerFaceAlpha',0.5);
xlim([0.45 0.7]);
ylim([0.45 0.7]);
l = refline([1 0]);
l.Color = "black";
l.LineWidth = 2;
title('Simplified prediction, including u','Fontsize',18);
xlabel('Predicted threshold','Fontsize',14);
ylabel('Most abundant threshold in simulation','Fontsize',14);

% Plot simplest prediction vs observed thressel
figure;
scatter(simple_threspred_2(phisel==1),thressel_t,15,'filled','MarkerFaceAlpha',0.5);
xlim([0.49 0.7]);
ylim([0.49 0.7]);
%xlim([0 1]);
%ylim([0 1]);
l = refline([1 0]);
l.Color = "black";
l.LineWidth = 2;
title('Simplified prediction, 1 / (1+\eta)','Fontsize',18);
xlabel('Predicted threshold','Fontsize',14);
ylabel('Most abundant threshold in simulation','Fontsize',14);
%}
