%%%%% Figures to show 2D selection of Phi_max and Thres values %%%%
%{
data = dlmread('~/PhageCom/Results/20191021_ThresPhiEvol_runs/Run2/p.txt');
% For Tt=7 and Tt=10, take newly solved (= correct) equilibrium
data(:,14) = dlmread('~/PhageCom/Results/20191021_ThresPhiEvol_runs/Run2/Extra_T7/p.txt');
data(:,20) = dlmread('~/PhageCom/Results/20191021_ThresPhiEvol_runs/Run2/Extra_T10/p.txt');
%}
exitflags = dlmread('~/PhageCom/Results/20191121_ThresPhiEvol_vary_Ttransfer/solve_exitflag.txt');
%data = dlmread('~/PhageCom/Results/20191121_ThresPhiEvol_vary_Ttransfer/p.txt');
data = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200819_ThresPhiEvol_fulltransfer/p.txt');

Ttransfer_range = [0.5:0.5:12,24];
nTt = length(Ttransfer_range);
phimax = 0:0.05:1;
nphi  = length(phimax);
thres = 0:0.05:1;
nthres = length(thres);

%%% Multiple heat maps %%%

nrows = floor(sqrt(nTt));
ncols = ceil(nTt / nrows);

heatmapdata = zeros(nphi,nthres,nTt);
for i=1:nthres
    heatmapdata(:,i,:) = data((1+(i-1)*nphi):(i*nphi),:);
end

cfig = figure;
set(cfig,'Units','centimeters','Position',[2 2 42 30],'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[42, 30])
for i=1:nTt
    subplot(nrows,ncols,i)
    colormap(flipud(bone));
    imagesc(heatmapdata(:,:,i));
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

%print(cfig,'~/PhageCom/Figures/Result_plots/thresphi_heatmaps.pdf','-dpdf','-r0')


%%% Two single heatmaps to include in main text figure

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
xlabel('Response threshold (\theta)','Fontsize',10);
ylabel('Lysogeny propensity (\phi_{max})','Fontsize',10);

saveas(fig,'~/PhageCom/Figures/Result_plots/thresphi_heatmap_T2.svg')

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
xlabel('Response threshold (\theta)','Fontsize',10);
ylabel('Lysogeny propensity (\phi_{max})','Fontsize',10);

saveas(fig,'~/PhageCom/Figures/Result_plots/thresphi_heatmap_T12.svg')

%%% Means with s.d. %%%
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




