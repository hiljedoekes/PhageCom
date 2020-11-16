%%%% Plot results of parameter sweep runs %%%%

clear all;

datafolder1 = '../Data/20200210_ThresPhiEvol_paramsweep1';
datafolder2 = '../Data/20200210_ThresPhiEvol_paramsweep2';

% Read data - combine results of both data folders in one matrix
S = dlmread([datafolder1,'/S.txt']);
temp = dlmread([datafolder2,'/S.txt']);
S = [S temp];
p = dlmread([datafolder1,'/p.txt']);
temp = dlmread([datafolder2,'/p.txt']);
p = [p temp];

nsample = length(S);

%%% Check if ODE integrations were successful (S = -1 if not)
check = 0;
for i=1:nsample
    if S(i) == -1
        disp(['ODE integration failed in run ',num2str(i)]);
        check = check+1;
    end
end
if check == 0
    disp('ODE integration successful in all runs');
end

%%% Characteristics of phage strains included in each run
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

%%% Find the most prevalent strain and corresponding phimax and thres values for each run. 
%%% + Find the mean value of phimax and thres in each run
phisel = zeros(1,nsample);
thressel = zeros(1,nsample);

% Find index of most prevalent phage for each run
[max_p, max_p_ind] = max(p,[],1);
% Translate back to phimax-value
indphi = mod(max_p_ind,nthres);
for j=1:nsample
    if indphi(j) == 0
        phisel(j) = phimax(nphi);
        indphi(j) = nphi;
    else
        phisel(j) = phimax(indphi(j));
    end
end
% Translate back to thres-value
indthres = floor((max_p_ind-1)/nphi)+1;
for j=1:nsample
        thressel(j) = thres(indthres(j));
end

% Find mean value of phi and thres for each run
phimean = phimax_all*p;
thresmean = thres_all*p;


%%% PLOT %%%

% Data structure for heatmap / bubble plot
heatmapdata = zeros(nphi,nthres);
for j=1:nsample
    heatmapdata(indphi(j),indthres(j)) = heatmapdata(indphi(j),indthres(j)) + 1;
end

% Bubble plot showing distribution of the dominant phage strain
[Xbins, Ybins] = meshgrid(0:0.1:1,0:0.1:1);
ind = heatmapdata(:) > 0;
figure;
scatter(Xbins(ind), Ybins(ind), 4*(heatmapdata(ind).^1.0),'filled')
grid(gca,'on')
xlim([0 1]);
ylim([0 1]);
xlabel('Response threshold (\theta)','Fontsize',14);
ylabel('Lysogeny propensity (\phi_{max})','Fontsize',14);


%%% EXTRA ANALYSIS %%%
%{
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


% Plot mean values of phi and thres - instead of most prevalent strain
figure;
scatter(thresmean,phimean,16,'filled','MarkerFaceAlpha',0.5);
xlim([0 1])
ylim([0 1])
xlabel('Mean Response threshold (\theta)','Fontsize',14);
ylabel('Mean Lysogeny propensity (\phi_{max})','Fontsize',14);
title('Results of 500 parameter sweep runs','Fontsize',14);
%}
