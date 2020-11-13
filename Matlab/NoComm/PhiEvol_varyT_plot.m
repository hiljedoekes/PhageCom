%%%%%%%% Construct heatmap figure of PhiEvol data %%%%%%%%

plotdata = dlmread('../Data/20191101_PhiEvol_varyT/p.txt');
%plotdata = dlmread('../Data/20200819_PhiEvol_fulltransfer/p.txt');

cfig=figure;
set(cfig,'Units','centimeters','Position',[2 2 25 10])
clf;
colormap(flipud(bone));
imagesc(plotdata);
cbar=colorbar;
cbar.Label.String = 'Relative strain frequency';
set(cbar,'FontSize',14);
set(gca,'YDir','normal');
set(gca,'Xtick',[2:2:24,25],'FontSize',14);
set(gca,'XtickLabel',[1:12,24]);
set(gca,'Ytick',1:10:41,'FontSize',14);
set(gca,'YtickLabel',0:0.05:0.2);
xlabel('Time between transfers (hours)','Fontsize',16);
ylabel('Lysogeny propensity','Fontsize',16);
