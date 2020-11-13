%%% Analysis of 20200619 data

% Read data
K = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200619_ThresPhiEvol_varyingK/K.txt');
exitflag = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200619_ThresPhiEvol_varyingK/solve_exitflag.txt');
P = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200619_ThresPhiEvol_varyingK/P.txt');
L = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200619_ThresPhiEvol_varyingK/L.txt');

% Remove column(s) of data that failed 
% (because gamma distribution with var=0 doesn't work)
ns = 441;
nt = 1000;
nruns = 20;
K = K(:,2:end);
Pfull = zeros(ns,nt,nruns);
Lfull = zeros(ns,nt,nruns);
for k=1:nruns
    Pfull(:,:,k) = P(:,(1000*k+1):(1000*(k+1)));
    Lfull(:,:,k) = L(:,(1000*k+1):(1000*(k+1)));
end

% For a given time point, plot the phage distribution for each variance
plottime = 900;
vars = 0.05:0.05:1;

p = zeros(ns,nruns);
l = zeros(ns,nruns);
for j=1:nruns
    p(:,j) = Pfull(:,plottime,j) / sum(Pfull(:,plottime,j));
    l(:,j) = Lfull(:,plottime,j) / sum(Lfull(:,plottime,j));
end

phimax = 0:0.05:1;
nphi  = length(phimax);
thres = 0:0.05:1;
nthres = length(thres);

nrows = floor(sqrt(nruns));
ncols = ceil(nruns / nrows);

heatmapdata = zeros(nphi,nthres,nruns);
for i=1:nthres
    heatmapdata(:,i,:) = p((1+(i-1)*nphi):(i*nphi),:);
end

cfig = figure;
set(cfig,'Units','centimeters','Position',[2 2 42 30],'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[42, 30])
for i=1:nruns
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
    title(['Var = ',num2str(vars(i))],'FontSize',12);
    xlabel('Response threshold (\theta)','Fontsize',10);
    ylabel('Lysogeny propensity (\phi_{max})','Fontsize',10);
end
%print(cfig,'~/PhageCom/Figures/Result_plots/thresphi_heatmaps.pdf','-dpdf','-r0')


% NEWNEWNEWNEWNEWNEWNEWNEWNEW %

ns = 441;
nt = 1000;
phimax = 0:0.05:1;
nphi  = length(phimax);
thres = 0:0.05:1;
nthres = length(thres);

varvec = [0.025:0.025:0.2,0.25:0.05:1];
CVvec = sqrt(varvec);
n = length(varvec);
Kfinal = zeros(nt,n);
Pfinal = zeros(ns,nt,n);
Lfinal = zeros(ns,nt,n);


% Construct a useful dataframe

% Read data 1
K = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200619_ThresPhiEvol_varyingK/K.txt');
exitflag = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200619_ThresPhiEvol_varyingK/solve_exitflag.txt');
P = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200619_ThresPhiEvol_varyingK/P.txt');
L = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200619_ThresPhiEvol_varyingK/L.txt');
% Remove column(s) of data that failed 
% (because gamma distribution with var=0 doesn't work)
nruns = 20;
K = K(:,2:end);
Pfull = zeros(ns,nt,nruns);
Lfull = zeros(ns,nt,nruns);
for k=1:nruns
    Pfull(:,:,k) = P(:,(1000*k+1):(1000*(k+1)));
    Lfull(:,:,k) = L(:,(1000*k+1):(1000*(k+1)));
end

% Store these values in the correct columns of final K, P and L
for i=1:4
    Kfinal(:,2*i) = K(:,i);
    Pfinal(:,:,2*i) = Pfull(:,:,i);
    Lfinal(:,:,2*i) = Lfull(:,:,i);
end
for i=5:nruns
    Kfinal(:,i+4) = K(:,i);
    Pfinal(:,:,i+4) = Pfull(:,:,i);
    Lfinal(:,:,2+4) = Lfull(:,:,i);
end


% Read data 2
K = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200819_ThresPhiEvol_varyingK/K.txt');
exitflag = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200819_ThresPhiEvol_varyingK/solve_exitflag.txt');
P = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200819_ThresPhiEvol_varyingK/P.txt');
L = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20200819_ThresPhiEvol_varyingK/L.txt');
for i=0:3
    Kfinal(:,(1+2*i)) = K(:,(i+1));
    Pfinal(:,:,(1+2*i)) = P(:,(1000*i+1):(1000*(i+1)));
    Lfinal(:,:,(1+2*i)) = L(:,(1000*i+1):(1000*(i+1)));
end

% Calculate relative frequencies
timevec = 900:1000;
lt = length(timevec);
pfinal = zeros(ns,n);
lfinal = zeros(ns,n);
for i=1:n
    pfinal(:,i) = (sum(Pfinal(:,timevec,i),2))/(sum(sum(Pfinal(:,timevec,i),2)));
    lfinal(:,i) = (sum(Lfinal(:,timevec,i),2))/(sum(sum(Lfinal(:,timevec,i),2)));
end

% Heatmap data
heatmapdata = zeros(nphi,nthres,n);
for i=1:nthres
    heatmapdata(:,i,:) = pfinal((1+(i-1)*nphi):(i*nphi),:);
end

%%% Figuur - supplement
nrows = 3;
ncols = 4;
nplots = nrows*ncols;

cfig = figure;
set(cfig,'Units','centimeters','Position',[2 2 44 30],'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[42, 30])
for i=[1:(nplots-1),n]
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
%print(cfig,'~/PhageCom/Figures/Result_plots/thresphi_heatmaps.pdf','-dpdf','-r0')


% Seperate panels for main text figure

% A) Var = 0.05
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

% B) Var = 0.15
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

% C) Var = 0.35
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



%%%% Analysis of 20201105 data

% Read data
K = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20201105_ThresPhiEvol_varyingK_u1/K.txt');
exitflag = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20201105_ThresPhiEvol_varyingK_u1/solve_exitflag.txt');
P = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20201105_ThresPhiEvol_varyingK_u1/P.txt');
L = dlmread('~/Documents/PhD/PhageCom/Revision_NewSimuls/20201105_ThresPhiEvol_varyingK_u1/L.txt');

ns = 441;
nt = 1000;
phimax = 0:0.05:1;
nphi  = length(phimax);
thres = 0:0.05:1;
nthres = length(thres);

varvec = [0.025:0.025:0.2,0.25:0.05:1];
CVvec = sqrt(varvec);
n = length(varvec);

for k=0:(n-1)
    Pfinal(:,:,(k+1)) = P(:,(1000*k+1):(1000*(k+1)));
    Lfinal(:,:,(k+1)) = L(:,(1000*k+1):(1000*(k+1)));
end

% Calculate relative frequencies
timevec = 900:1000;
lt = length(timevec);
p = zeros(ns,n);
l = zeros(ns,n);
for i=1:n
    p(:,i) = (sum(Pfinal(:,timevec,i),2))/(sum(sum(Pfinal(:,timevec,i),2)));
    l(:,i) = (sum(Lfinal(:,timevec,i),2))/(sum(sum(Lfinal(:,timevec,i),2)));
end

% Heatmap data
heatmapdata = zeros(nphi,nthres,n);
for i=1:nthres
    heatmapdata(:,i,:) = p((1+(i-1)*nphi):(i*nphi),:);
end

%%% Figuur - supplement
nrows = 3;
ncols = 4;
nplots = nrows*ncols;

cfig = figure;
set(cfig,'Units','centimeters','Position',[2 2 44 30],'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[42, 30])
for i=[1:(nplots-1),n]
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