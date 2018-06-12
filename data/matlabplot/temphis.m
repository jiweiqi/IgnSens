% clear;clc;close all;
% D1 = dlmread('test0.out')';
% 
% D2 = dlmread('test.out' )';
% 
% diff1 = diff(D1(:,2))./diff(D1(:,1));
% index = find( diff1 == max(diff1) );
% ign1 = D1(index, 1)
% 
% diff1 = diff(D2(:,2))./diff(D2(:,1));
% index = find( diff1 == max(diff1) );
% ign2 = D2(index, 1)
% 
% figure();
% 
% i1 = find( D1(:,2) > 2000 & D1(:,2) );
% i2 = find( D2(:,2) > 2000 & D2(:,2) );
% plot(D1(i1,1), D1(i1,2), '-r');
% hold all;
% plot(D2(i2,1), D2(i2,2), '--b');
% 
% return;

clear;clc;close all;
%
% datadir = 'C:\Users\weiqi\Dropbox\UQ_Manuscript\Manuscript\PCI_Sensitivity\data\CH420atm1000K';
% condstr = ' CH4, \phi=1, 20atm, 1000K ';
% gas = GRI30;

% datadir = ...
%     'C:\Users\weiqi\Dropbox\UQ_Manuscript\Manuscript\PCI_Sensitivity\data\DME20atm800K';
% condstr = ' DME, \phi=1, 20atm, 800K ';
% gas = IdealGasMix( '../mech/kin20285-chem.cti');

% datadir = ...
%     'C:\Users\weiqi\Dropbox\UQ_Manuscript\Manuscript\PCI_Sensitivity\data\NC7H1620atm800K';
% condstr = ' NC7H16, \phi=1, 20atm, 800K ';
% gas = IdealGasMix( '../mech/nc7sk88.cti');

datadir = ...
    'C:\Users\weiqi\Dropbox\UQ_Manuscript\Manuscript\PCI_Sensitivity\data\IC8H1820atm750K';
condstr = ' IC8H18, \phi=1, 20atm, 750K ';
gas = IdealGasMix( '../mech/ic8sk143.cti');

addpath(datadir);
% read time and temeprature
D1 = dlmread('senstime.out')';
time = D1(:,1);
T1 = D1(:,2);
diff1 = diff(D1(:,2))./diff(D1(:,1));
ignpos = find( diff1 == max(diff1) );
ign1 = D1(ignpos, 1);
nt = time./ign1;

iend = find( T1 > max(T1)-1.0 );
tmax = time(iend(1))/ign1;

% local temperature sensitivity
sens = dlmread('sens.out');
sensnorm = normr(sens);
sensign = sensnorm( ignpos, : );

% local species sensitivity
sensY = dlmread('sensH2O');
sensYnorm = normr(sensY);

% for i=1:size(sens,1)
%     if sensnorm(i, 242) < 0
%         sensnorm(i, :) = -sensnorm(i, :);
%     end
% end

cos = [];
cosY = [];
amp = [];
for i=1:length( time )
    u = sensnorm(i, :)* sign( max(sensnorm(i, :))*max(sens(i, :)) ) ;
    uY = sensYnorm(i, :)* sign( max(sensYnorm(i, :))*max(sensY(i, :)) );
    
    cos = [cos; sensign*u' ];
    amp = [amp; norm( sens(i , :) )];
    
    cosY = [cosY; sensign*uY'];
end

%% validation vs finite difference
%% Sensitvity bar plot
sensFD = dlmread('idtFD.txt');
sensFDnorm = sensFD./norm(sensFD);

H = figure_pci(14.4, 7); 
figure( H );
A=axes('Parent',H, ...
'Position',[0.18 0.18 0.79 0.32], ...
'FontName','Times', ...
'YColor','k', ...
'YAxisLocation','left', ...
'Box','off');
n=25;
[~, sortedindex] = sort(abs(sensign), 'descend');
barwidth = 0.4;

if sensign(sortedindex(1)) > 0
    h = bar(A, [1:n]-0.22, sensign(sortedindex(1:n)), barwidth);
else
    h = bar(A,[1:n]-0.22,  -sensign(sortedindex(1:n)), barwidth);
end
hold all;
h1 = bar(A, [1:n]+0.22, sensFDnorm(sortedindex(1:n)), barwidth);

A.Position = [0.08,0.57,0.89,0.42];
ylabel(A, 'Component', 'Position', [-1.0, 0]);
A.YLim = [-0.70, 0.70];
A.YTick = -0.5:0.5:0.5;
A.XLim = [0.5, n+0.5];
A.XTick = 1:n;
legend(A, 'u_\tau', 'u_{T,ign}');
A.Legend.Box='off';
% xlim([-0.6, 0.6]);
% xlabel('Reaction index');
% ylim( [0.5, n+0.5] );
A.FontSize = 8;A.FontName='Times';

B=axes('Parent',H, ...
'Position',[0.08,0.16,0.89,0.42], ...
'FontName','Times', ...
'YColor','k', ...
'YAxisLocation','left', ...
'Box','off');

h2 = bar(B, -(sensFDnorm( sortedindex(1:n) )' + sensign(sortedindex(1:n)))./sensign(sortedindex(1:n)) );
B.XTick = 1:n;
B.XTickLabel = sortedindex(1:n);
B.XTickLabelRotation=90;
xlabel('Reaction index');
B.Position = [0.08,0.14,0.89,0.42];
B.YLim = [ -0.1, 0.1 ];
B.YTick = -0.1:0.1:0.1;
B.XLim = [0.5, n+0.5];
ylabel(B, 'Relatvie error', 'Position', [-1.0, 0]);
B.FontSize = 8;B.FontName='Times';
annotation(H,'textbox',[0.45 0.85 0.1 0.1],...
    'String',condstr,  'linestyle','none','FontSize',9,'FontName','Times');

%% Sensitvity bar plot
H = figure_pci(14.4, 5); 
figure( H );
A=axes('Parent',H, ...
'Position',[0.18 0.18 0.79 0.79], ...
'FontName','Times', ...
'YColor','k', ...
'YAxisLocation','left', ...
'Box','on');
n=9;
[~, sortedindex] = sort(abs(sensign), 'descend');
if sensign(sortedindex(1))< 0
    h = barh(A, sensign(sortedindex(1:n)));
else
    h = barh(A, -sensign(sortedindex(1:n)));
end
yticklabel =  strrep(reactionEqn(gas,sortedindex(1:n)), ' ', '');
for i=1:n
    yticklabel{i} = strcat( yticklabel{i}, ' (R' ,num2str(sortedindex(i)), ')' );
end
A.YTickLabel = yticklabel;
xlim([-0.6, 0.6]);
xlabel('component');
A.Position = [0.5,0.18,0.48,0.79];
ylim( [0.5, n+0.5] );
A.FontSize = 8;A.FontName='Times';
% annotation(H,'textbox',[0.12 0.1 0.1 0.1],...
%     'String',{'CH4, \phi=1,20 atm,1000K'},  'linestyle','none','FontSize',9,'FontName','Times');
annotation(H,'textbox',[0.15 0.1 0.1 0.1],...
    'String',condstr,  'linestyle','none','FontSize',9,'FontName','Times');

annotation(H,'textbox',[0.15 0.2 0.1 0.1],...
    'String','(a)',  'linestyle','none','FontSize',9,'FontName','Times');

%% Evolution of components
H = figure_pci(6.7, 6); 
figure( H );
A=axes('Parent',H, ...
'Position',[0.18 0.18 0.79 0.79], ...
'FontName','Times', ...
'YColor','k', ...
'YAxisLocation','left', ...
'Box','on');
n=5;
% stylelist = { '-k', '--r', '-b', ':c', '--+m' };
stylelist = { '-', '--', '-.', ':', '--', '-','-', '--', '-', ':', '--.', '-' };
[sortedvariable, sortedindex] = sort(abs(sensign), 'descend');

for i=1:n
    j = sortedindex(i);
%     plot( nt, sensnorm(:, j).*(-sign(max(sensnorm(:,271)))), stylelist{i} ,'LineWidth', 2 ,'DisplayName', ['R', num2str( j )] );
%     plot( nt, sens(:, j).*(-sign(max(sensnorm(:,2)))), stylelist{i} ,'LineWidth', 2 ,'DisplayName', ['R', num2str( j )] );
    plot( nt, -sensnorm(:, j), stylelist{i} ,'LineWidth', 2 ,'DisplayName', ['R', num2str( j )],'MarkerSize',3 );
    hold all;
end
legend show;
legend boxoff;
A.XLim = [0, 1.05];
A.FontSize = 8;
xlabel( 'NT');
ylabel('component');
% annotation(H,'textbox',[0.12 0.1 0.1 0.1],...
%     'String',{'CH4, \phi=1,20 atm,1000K'},  'linestyle','none','FontSize',9,'FontName','Times');
annotation(H,'textbox',[0.35 0.17 0.1 0.1],...
    'String',condstr,  'linestyle','none','FontSize',9,'FontName','Times');

%% Amp of temeprature sensitivity
H = figure_pci(7.2, 5.0); 
figure( H );
A=axes('Parent',H, ...
'Position',[0.18 0.18 0.79 0.79], ...
'FontName','Times', ...
'YColor','k', ...
'YAxisLocation','left', ...
'Box','on');
h = plot(nt(2:end), diff1, '-k' ,'LineWidth', 2,'MarkerSize', 4);
hold all;
h1 = plot(nt(1:end), amp, '--r' ,'LineWidth', 2,'MarkerSize', 4);
A.XLim = [0,1.05];
A.YScale = 'log';
A.YLim = [1E-10, 1E10]; %DME 800 K
% A.YLim = [1E-5, 1E10]; %CH4 1000K
xlabel(A, 'NT');
ylabel(A, 'dT/dt, \lambda', 'Position', [-0.15, 1] );

yyaxis right
h2 = plot(nt(1:end), T1, '-.b' ,'LineWidth', 1,'MarkerSize', 4);
% ylabel('T (K)', 'Position', [1.30, 1800]);
ylabel('T (K)', 'Position', [1.23, 2000]);
legend(A, 'dT/dt', '\lambda', 'T');
legend boxoff;
A.Position = [0.16, 0.17, 0.66, 0.80];
annotation(H,'textbox',[0.22 0.87 0.1 0.1],...
    'String',condstr,  'linestyle','none','FontSize',9,'FontName','Times');
annotation(H,'textbox',[0.22 0.87 0.1 0.1],...
    'String','(a)',  'linestyle','none','FontSize',9,'FontName','Times');
A.FontSize = 8; A.FontName = 'Times';

%% plot the inner product
H = figure_pci(7.2, 5.0); 
figure( H );
A=axes('Parent',H, ...
'Position',[0.18 0.18 0.66 0.66], ...
'FontName','Times', ...
'YColor','k', ...
'YAxisLocation','left', ...
'Box','on');
% yyaxis left
set(gca,'linestyleorder',{'-','--','-.',':','*','+'},...
'colororder',[0 0 1;0 .5 0;1 0 0],'nextplot','add')

h1 = plot(nt, cos, '-k' ,'LineWidth', 2,'MarkerSize', 4, 'DisplayName', 'T');
A.YLim = [0,1];
A.XLim = [0,tmax];
hold all;

stylelist = { '-', '--', '-.', ':', '-', '--','-.', ':', '-', '--', '-.', ':' };
j = 1;
for sp = {'CO2', 'CO', 'H2O', 'OH','H2O2'}
% for sp = {'CH4', 'CO2', 'CO', 'H2O', 'OH', 'H2O2'}
% for sp = {'CH3OCH3', 'CH2O', 'CO2', 'H2O', 'OH', 'H2O2'}
% for sp = {'H2O'}
    j = j+1;
    fname = [ 'sens', sp{1}];
    sensY = dlmread(fname);
    sensYnorm = normr(sensY);
    cosY = [];
    ampY = [];
    for i=1:length( time )
        uY = sensYnorm(i, :);
        cosY = [cosY; sensign*uY'];
        ampY = [ampY;  norm(sensY(i, :)) ];
    end
    
    h = plot(nt, abs(cosY), stylelist{j},'LineWidth', 1,'MarkerSize', 4, 'DisplayName', sp{1});
    hold all;
end
legend show;
legend boxoff;
xlabel(A, 'NT' );
ylabel(A, '<u, u_{T,ign}>');
A.FontSize = 8;A.FontName = 'Times';
annotation(H,'textbox',[0.35 0.18 0.1 0.1],'String',condstr,  'linestyle','none','FontSize',9,'FontName','Times');
A.Position = [0.16, 0.17, 0.82, 0.81];
annotation(H,'textbox',[0.18 0.87 0.1 0.1],...
    'String','(a)',  'linestyle','none','FontSize',9,'FontName','Times');

yyaxis right
plot(nt(2:end), log(diff1) , '-r' ,'LineWidth', 2,'MarkerSize', 4, 'DisplayName', 'd[OH]/dt');
hold all;
plot( [0,1], [0, 0], '-m',  'DisplayName', 'y=0' );

return

% a=sensnorm(:,[270, 271, 272, 247, 50, 239, 19, 243, 273, 43, 17, 56, 51, 45]+1);
% a=sensnorm(:,[241	247	270	51	19	46	243	271	56	43	244	239	272	45	273]+1);

figure();
yyaxis left
plot( time(1:end)./ign1, abs(cos), '-ob' );
hold all;
plot( time(1:end)./ign1, abs(cosch4), '--sr' );
yyaxis right
plot( time(1:end)./ign1, T1 );
xlim([0,1.05]);

return;

%% Components for HCCI Natural gas
datadir = 'C:\Users\weiqi\Dropbox\UQ_Manuscript\Manuscript\PCI_Sensitivity\data\ch4hcci';
condstr = 'Natural gas, HCCI';
gas = GRI30;
DATA = dlmread(fullfile(datadir, 'idthcci.txt') );
sensFD = DATA(:,1);
sensFDnorm = sensFD./norm(sensFD);
sensign = DATA(:,2)'./norm(DATA(:,2));
% sensignnorm = sensign./norm(sensign);
H = figure_pci(14.4, 7); 
figure( H );
A=axes('Parent',H, ...
'Position',[0.18 0.18 0.79 0.32], ...
'FontName','Times', ...
'YColor','k', ...
'YAxisLocation','left', ...
'Box','off');
n=25;
[~, sortedindex] = sort(abs(sensign), 'descend');
barwidth = 0.4;

if sensign(sortedindex(1)) < 0
    h = bar(A, [1:n]-0.2, sensign(sortedindex(1:n)), barwidth);
else
    h = bar(A,[1:n]-0.2,  -sensign(sortedindex(1:n)), barwidth);
end
hold all;
h1 = bar(A, [1:n]+0.2, sensFDnorm(sortedindex(1:n)), barwidth);

A.Position = [0.08,0.57,0.89,0.42];
ylabel(A, 'Component', 'Position', [-1.0, 0]);
A.YLim = [-0.65, 0.65];
A.YTick = -0.5:0.5:0.5;
A.XLim = [0.5, n+0.5];
A.XTick = 1:n;
legend(A, 'u_\tau', 'u_{T,ign}');
A.Legend.Box='off';
% xlim([-0.6, 0.6]);
% xlabel('Reaction index');
% ylim( [0.5, n+0.5] );
A.FontSize = 8;A.FontName='Times';

B=axes('Parent',H, ...
'Position',[0.08,0.16,0.89,0.42], ...
'FontName','Times', ...
'YColor','k', ...
'YAxisLocation','left', ...
'Box','off');

Error = -(sensFDnorm( sortedindex(1:n) )' - sensign(sortedindex(1:n)))./sensign(sortedindex(1:n));
h2 = bar(B, Error , barwidth);
B.XTick = 1:n;
B.XTickLabel = sortedindex(1:n);
B.XTickLabelRotation=90;
xlabel('Reaction index');
B.Position = [0.08,0.14,0.89,0.42];
B.YLim = [ -0.2, 0.2 ];
B.YTick = -0.1:0.1:0.1;
B.XLim = [0.5, n+0.5];
ylabel(B, 'Relatvie error', 'Position', [-1.0, 0]);
B.FontSize = 8;B.FontName='Times';
annotation(H,'textbox',[0.45 0.85 0.1 0.1],...
    'String',condstr,  'linestyle','none','FontSize',9,'FontName','Times');

%%
% figure();
% % yyaxis left
% plot( time(1:end)./ign1, amp./max(amp) );
% hold all;
% % yyaxis right
% plot( time(2:end)./ign1, diff1./max(diff1) );
% xlim([0,1.05]);
% ylim([0,1]);


% % plot the amp of temeprature sensitivity
% H = figure_pci(6.7, 6.0); 
% figure( H );
% A=axes('Parent',H, ...
% 'Position',[0.18 0.18 0.79 0.79], ...
% 'FontName','Times', ...
% 'YColor','k', ...
% 'YAxisLocation','left', ...
% 'Box','on');
% heatrelease = dlmread('heatrelease.out');
% h = plot(nt, -amp./heatrelease, '-k' ,'LineWidth', 2,'MarkerSize', 4);YOH
% hold all;
% % h1 = plot(nt(1:end), amp, '--r' ,'LineWidth', 2,'MarkerSize', 4);
% A.XLim = [0,1.05];
% A.YScale = 'log';
% % A.YLim = [1E-10, 1E10]; %DME 800 K
% % A.YLim = [1E-5, 1E10]; %CH4 1000K
% xlabel(A, 'NT');
% ylabel(A, 'dt/dT(\lambda)', 'Position', [-0.15, 1] );
% 
% yyaxis right
% h2 = plot(nt(1:end), T1, '-.b' ,'LineWidth', 1,'MarkerSize', 4);
% ylabel('T (K)', 'Position', [1.30, 1800]);
% legend(A, 'dT/dt', '\lambda', 'T');
% legend boxoff;
% A.FontSize = 8;
% A.Position = [0.16, 0.16, 0.66, 0.80];
% annotation(H,'textbox',[0.35 0.17 0.1 0.1],...
%     'String',condstr,  'linestyle','none','FontSize',9,'FontName','Times');
