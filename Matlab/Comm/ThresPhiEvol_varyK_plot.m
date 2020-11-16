%%% Analyse and plot the results of ThresPhiEvol-runs with long time
%%% between transfers (T = 24h) and a varying degree of noise on the
%%% bacterial carrying capacity K.

clear all;
datafolder = '../Data/20200819_ThresPhiEvol_varyK';

%%% Read data
K = dlmread([datafolder,'/K.txt']);
ef = dlmread([datafolder,'/solve_exitflag.txt']);
Plong = dlmread([datafolder,'/P.txt']);

% Check if all runs were successful
check = 0;
for i=1:length(ef)
    if ef(i) ~= 0
        disp(['Problems in run ',num2str(i)']);
        check = check+1;
    end
end
if check == 0
    disp('All runs successful.');  
end

%%% Characteristics of runs and phage strains included

nt = 100;       % number of transfers that data were saved for

% Phage strains
phimax = 0:0.05:1;      % phimax-values included
nphi  = length(phimax);
thres = 0:0.05:1;       % thres-values included
nthres = length(thres);
ns = nphi*nthres;       % number of phage strains

varvec = [0.025:0.025:0.2,0.25:0.05:1]; % Variance of K for different runs
CVvec = sqrt(varvec);                   % Vector of coefficients of variation
n = length(varvec);                     % Number of runs

%%% Store phage distributions in a useful dataframe (distribution per run)
P = zeros(ns,nt,n);
for k=1:n
    P(:,:,k) = Plong(:,(nt*(k-1)+1):(nt*k));
end

%%% Calculate relative strain frequencies in the last nt timepoints
%%% Store these in a heatmap data format
p = zeros(ns,n);
for i=1:n
    p(:,i) = (sum(P(:,:,i),2))/(sum(sum(P(:,:,i),2)));
end
heatmapdata = zeros(nphi,nthres,n);
for i=1:nthres
    heatmapdata(:,i,:) = p((1+(i-1)*nphi):(i*nphi),:);
end

%%% PLOT results (supporting figure)

% Control the number of subplots:
nrows = 3;
ncols = 4;
nplots = nrows*ncols;

cfig = figure;
set(cfig,'Units','centimeters','Position',[2 2 44 30],'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[42, 30])
for i=[1:(nplots-1),n]      % Include increasing variance up to some point, and the distribution for the highest variance included
    if i==n
        subplot(nrows,ncols,nplots)
    else
        subplot(nrows,ncols,i)
    end
    colormap(flipud(bone));
    imagesc(heatmapdata(:,:,i));
    cbar=colorbar;
    set(cbar,'FontSize',12);
    set(gca,'YDir','normal');
    set(gca,'Xtick',1:10:nthres,'FontSize',12);
    set(gca,'XtickLabel',thres(1:10:nthres));
    set(gca,'Ytick',1:10:nphi,'FontSize',12);
    set(gca,'YtickLabel',phimax(1:10:nphi));
    title(['Var = ',num2str(varvec(i)),'; CV = ',num2str(round(CVvec(i),2))],'FontSize',14);
    xlabel('Response threshold (\theta)','Fontsize',14);
    ylabel('Lysogeny propensity (\phi_{max})','Fontsize',14);
end


%%% Separate panels for main text figure

% A) Var = 0.05 (subplot 2)
cfig = figure;
set(cfig,'Units','centimeters','Position',[2 2 12.5 10],'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[42, 30])
colormap(flipud(bone));
imagesc(heatmapdata(:,:,2));
cbar=colorbar;
set(cbar,'FontSize',12);
set(gca,'YDir','normal');
set(gca,'Xtick',1:10:nthres,'FontSize',12);
set(gca,'XtickLabel',thres(1:10:nthres));
set(gca,'Ytick',1:10:nphi,'FontSize',12);
set(gca,'YtickLabel',phimax(1:10:nphi));
title(['Var = ',num2str(varvec(2)),'; CV = ',num2str(CVvec(2))],'FontSize',14);
xlabel('Response threshold (\theta)','Fontsize',14);
ylabel('Lysogeny propensity (\phi_{max})','Fontsize',14);

% B) Var = 0.15 (Subplot 6)
cfig = figure;
set(cfig,'Units','centimeters','Position',[2 2 12.5 10],'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[42, 30])
colormap(flipud(bone));
imagesc(heatmapdata(:,:,6));
cbar=colorbar;
set(cbar,'FontSize',12);
set(gca,'YDir','normal');
set(gca,'Xtick',1:10:nthres,'FontSize',12);
set(gca,'XtickLabel',thres(1:10:nthres));
set(gca,'Ytick',1:10:nphi,'FontSize',12);
set(gca,'YtickLabel',phimax(1:10:nphi));
title(['Var = ',num2str(varvec(6)),'; CV = ',num2str(CVvec(6))],'FontSize',14);
xlabel('Response threshold (\theta)','Fontsize',14);
ylabel('Lysogeny propensity (\phi_{max})','Fontsize',14);

% C) Var = 0.35 (Subplot 11)
cfig = figure;
set(cfig,'Units','centimeters','Position',[2 2 12.5 10],'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[42, 30])
colormap(flipud(bone));
imagesc(heatmapdata(:,:,11));
cbar=colorbar;
set(cbar,'FontSize',12);
set(gca,'YDir','normal');
set(gca,'Xtick',1:10:nthres,'FontSize',12);
set(gca,'XtickLabel',thres(1:10:nthres));
set(gca,'Ytick',1:10:nphi,'FontSize',12);
set(gca,'YtickLabel',phimax(1:10:nphi));
title(['Var = ',num2str(varvec(11)),'; CV = ',num2str(CVvec(11))],'FontSize',14);
xlabel('Response threshold (\theta)','Fontsize',14);
ylabel('Lysogeny propensity (\phi_{max})','Fontsize',14);
