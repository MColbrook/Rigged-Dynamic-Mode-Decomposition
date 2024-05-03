clear
close all

%% Extract the vorticity field, KE and time

load('Cavity19k.mat') % data from https://ucsb.app.box.com/s/tbtdij1r4z5qt7o4cbgmolx28do98xci
[Grid,Operators] = CavityGridOperators(VelocityField.N);
VORT = - Operators.del2*real(VelocityField.Psi);
VORT = VORT - mean(VORT(:));
KE = VelocityField.KE; KE = KE - mean(KE);
Time = VelocityField.Time - VelocityField.Time(1);
T = range(Time);

%% Extract Koopman modes
M = 10000;
ind1 = 1:M;
DATA = VORT(:,[ind1,M+1]);
DATA = DATA - mean(DATA,2);
NN =randn(size(DATA));%/mean(std(DATA'));
a = std(DATA(:));
NN = a'.*NN;
DATA2 = DATA+0.2*NN;
[~,S,V]=svd(transpose(DATA(:,ind1)),'econ');
[~,S2,V2]=svd(transpose(DATA2(:,ind1)),'econ');

%%
N =200;
PX=transpose(DATA(:,ind1))*V(:,1:N);
PY=transpose(DATA(:,ind1+1))*V(:,1:N);
g_coeffs = PX\transpose(DATA(:,ind1));

PX2=transpose(DATA2(:,ind1))*V2(:,1:N);
PY2=transpose(DATA2(:,ind1+1))*V2(:,1:N);
g_coeffs2 = PX2\transpose(DATA2(:,ind1));

%% Run riggedDMD
epsilon = 0.01;
order = 6;

TH = (1:4)*2*pi*(Time(2)-Time(1))*0.1598;
gMODES = riggedDMD_MODES(PX,PY,1/M,epsilon,TH,[],'order',order,'g_coeffs',g_coeffs);
gMODES2 = riggedDMD_MODES(PX2,PY2,1/M,epsilon,TH,[],'order',order,'g_coeffs',g_coeffs2);

%
[xx2,yy2] = meshgrid(Grid.x,Grid.x);
xxx = xx2; yyy = yy2;
[xxx,yyy] = meshgrid(-1:.005:1,-1:.005:1);
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(2,4,"TileSpacing","compact", 'TileIndexing', 'columnmajor')


for jj=1:length(TH)
    C = squeeze(gMODES(:,:,jj));
    cc = transpose(mean(transpose(C)));
    u = cc'*C; u = u - mean(u(:));

    C = squeeze(gMODES2(:,:,jj));
    cc = transpose(mean(transpose(C)));
    u2 = cc'*C; u2 = u2 - mean(u2(:));

    aa = mean(angle(u(abs(u)>std(u))./u2(abs(u)>std(u))));
    u = real(u); u = u - mean(u(:));
    u2 = real(exp(1i*aa)*u2); u2 = u2 - mean(u2(:));


    u(u(:)>3*std(u(:))+mean(u(:))) = 3*std(u(:))+mean(u(:));
    u(u(:)<-3*std(u(:))+mean(u(:))) = -3*std(u(:))+mean(u(:));
    u = reshape(real(u),VelocityField.N+1,VelocityField.N+1);
    u = interp2(xx2,yy2,u,xxx,yyy,'cubic');

    u2(u2(:)>3*std(u2(:))+mean(u2(:))) = 3*std(u2(:))+mean(u2(:));
    u2(u2(:)<-3*std(u2(:))+mean(u2(:))) = -3*std(u2(:))+mean(u2(:));
    u2 = reshape(real(u2),VelocityField.N+1,VelocityField.N+1);
    u2 = interp2(xx2,yy2,u2,xxx,yyy,'cubic');

    
    nexttile
    contourf(xxx,yyy,u,120,'LineStyle','None')     % vorticty
    axis square equal; 
    box on;
    colormap(brewermap([],'RdYlBu'))
    clim([-3*std(u(:))+mean(u(:)),3*std(u(:))+mean(u(:))])
    ax=gca; ax.FontSize=14;

    nexttile
    contourf(xxx,yyy,u2,120,'LineStyle','None')     % vorticty
    axis square equal;  
    box on;
    colormap(brewermap([],'RdYlBu'))
    clim([-3*std(u2(:))+mean(u2(:)),3*std(u2(:))+mean(u2(:))])
    ax=gca; ax.FontSize=14;

end
exportgraphics(gcf,sprintf('cavity_noise_3.png',jj),'ContentType','vector','BackgroundColor','none')
