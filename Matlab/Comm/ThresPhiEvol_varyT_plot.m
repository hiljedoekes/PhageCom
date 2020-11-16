%%%%% Plot results of ThresPhiEvol runs with varying time between transfers %%%%

%data = dlmread('../Data/20191121_ThresPhiEvol_varyT/p.txt');
data = dlmread('../Data/20200819_ThresPhiEvol_fulltransfer/p.txt');

Ttransfer_range = [0.5:0.5:12,24];      % Range of Ttransfer values that were included
nTt = length(Ttransfer_range);
phimax = 0:0.05:1;                      % phimax-values included
nphi  = length(phimax);
thres = 0:0.05:1;                       % thres-values included
nthres = length(thres);

%%% Heatmap data format %%%
heatmapdata = zeros(nphi,nthres,nTt);
for i=1:nthres
    heatmapdata(:,i,:) = data((1+(i-1)*nphi):(i*nphi),:);
end

% 1. Heat map of strain distribution for each Ttransfer value 
% (supporting figures)

nrows = floor(sqrt(nTt));
ncols = ceil(nTt / nrows);

cfig = figure;
set(cfig,'Units','centimeters','Position',[2 2 42 30],'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[42, 30])
for i=1:nTt
    subplot(nrows,ncols,i)
    colormap(flipud(bone));         % Colour map
    imagesc(heatmapdata(:,:,i));    % Plot results
    cbar=colorbar;
    set(cbar,'FontSize',10);
    set(gca,'YDir','normal');
    set(gca,'Xtick',1:10:nthres,'FontSize',10);
    set(gca,'XtickLabel',thres(1:10:nthres));
    set(gca,'Ytick',1:10:nphi,'FontSize',10);
    set(gca,'YtickLabel',phimax(1:10:nphi));
    title(['T = ',num2str(Ttransfer_range(i)),' h'],'FontSize',12);
    xlabel('Response threshold (\theta)','Fontsize',10);
    ylabel('Lysogeny propensity (\phi_{max})','Fontsize',10);
end


%%% 2. Two single heatmaps to include in main text figure

% T = 2
fig = figure;
colormap(flipud(bone));
clims = [0 0.2];
imagesc(heatmapdata(:,:,4),clims);
cbar=colorbar;
set(cbar,'FontSize',10);
set(gca,'YDir','normal');
set(gca,'Xtick',1:10:nthres,'FontSize',10);
set(gca,'XtickLabel',thres(1:10:nthres));
set(gca,'Ytick',1:10:nphi,'FontSize',10);
set(gca,'YtickLabel',phimax(1:10:nphi));
title('T = 2','FontSize',12);
xlabel('Response threshold (\theta)','Fontsize',16);
ylabel('Lysogeny propensity (\phi_{max})','Fontsize',16);

% T = 12
fig = figure;
colormap(flipud(bone));
clims = [0 0.2];
imagesc(heatmapdata(:,:,24),clims);
cbar=colorbar;
set(cbar,'FontSize',10);
set(gca,'YDir','normal');
set(gca,'Xtick',1:10:nthres,'FontSize',10);
set(gca,'XtickLabel',thres(1:10:nthres));
set(gca,'Ytick',1:10:nphi,'FontSize',10);
set(gca,'YtickLabel',phimax(1:10:nphi));
title('T = 12','FontSize',12);
xlabel('Response threshold (\theta)','Fontsize',16);
ylabel('Lysogeny propensity (\phi_{max})','Fontsize',16);


%%% 3. Plot means with s.d. %%%
%{
meanphis = zeros(1,nTt);
sdphis = zeros(1,nTt);
for i=1:nTt
    totphifreq = sum(heatmapdata(:,:,i),2);
    meanphis(i) = phimax*totphifreq;
    sdphis(i) = sqrt(((phimax-meanphis(i)).^2)*totphifreq);
end

meanthres = zeros(1,nTt);
sdthres = zeros(1,nTt);
for i=1:nTt
    totthresfreq = sum(heatmapdata(:,:,i),1);
    meanthres(i) = thres*totthresfreq';
    sdthres(i) = sqrt(((thres-meanthres(i)).^2)*totthresfreq');
end

cfig=figure;
set(cfig,'Units','centimeters','Position',[2 2 30 20])
hold on;
errorbar(1:nTt,meanphis,sdphis,'-o','Color','blue','MarkerSize',8,'MarkerFaceColor','blue')
errorbar(1:nTt,meanthres,sdthres,'-o','Color','red','MarkerSize',8,'MarkerFaceColor','red')
set(gca,'Xtick',[1,2:2:24,25],'FontSize',12);
set(gca,'XtickLabel',[0.5,1:12,24]);
set(gca,'Ytick',0:0.2:1.0,'FontSize',12);
xlabel('Time between transfers (hours)','Fontsize',14);
ylabel('Mean Phi max (blue) and Response threshold (red), +/- s.d.','Fontsize',14);
%}




