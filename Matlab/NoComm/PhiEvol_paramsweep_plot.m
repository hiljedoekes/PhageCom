%%%% Plot results of paramsweep against predicted phi* %%%%

datafolder = '../Data/20200210_PhiEvol_paramsweep';

% Read data
a = dlmread([datafolder,'/a.txt']);
B = dlmread([datafolder,'/B.txt']);
dilution_factor = dlmread([datafolder,'/dilution_factor.txt']);
alpha = dlmread([datafolder,'/alpha.txt']);
deltaP = dlmread([datafolder,'/deltaP.txt']);

P = dlmread([datafolder,'/P.txt']);
phisel = dlmread([datafolder,'/phisel.txt']);

% Calculate estimate of P0
P0 = dilution_factor.*((alpha.*(1-alpha))./(deltaP + a.*(1-alpha)));
%p0 = dilution_factor.*sum(P,1);  % Alternatively, use data to calculate actual P0

% Calculate predicted phi*
eta = 1 - ones(1,500)./B;
phipred_simple = - eta ./ (log(P0));
%phipred = ( eta.*(eta+P0).*log(eta./P0 + 1) ) ./ ( -eta + (eta+P0).* (log(eta./P0 + 1)).^2 ) ; % Alternatively, use this full prediction (see Mathematica nb for derivation)

% Plot
figure;
scatter(phipred_simple, phisel,36,'filled','MarkerFaceAlpha',0.5);
xlim([0.0 0.15]);
ylim([0.0 0.15]);
l = refline([1 0]);
l.Color = "black";
l.LineWidth = 2;
xlabel('Predicted phi','Fontsize',14);
ylabel('Selected phi in simulation','Fontsize',14);

%%%%% ADDITIONAL ANALYSIS %%%%%
%{
% Explanatory effect of eta and P0
figure;
subplot(1,2,1);
scatter(eta, phisel,36,'filled','MarkerFaceAlpha',0.4);
xlim([0 1]);
ylim([0 0.15]);
xlabel('eta = 1 - 1/B', 'Fontsize',14);
ylabel('Selected phi in simulation','Fontsize',14);
subplot(1,2,2);
scatter(-(log(P0)).^(-1), phisel,36,'filled','MarkerFaceAlpha',0.4);
%xlim([0 1]);
ylim([0 0.15]);
xlabel('- 1/log(P0)', 'Fontsize',14);
ylabel('Selected phi in simulation','Fontsize',14);


% Perform linear regression
x = (phipred(phisel<1))';           % Exclude the two values that have phisel == 1
y = (phisel(phisel<1))';
X = [ones(length(x),1) x];
b = X\y
yCalc = X*b;
Rsq = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2)

% Regression without intercept
b1 = x\y
yCalc1 = b1*x;
Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2) 

%}

