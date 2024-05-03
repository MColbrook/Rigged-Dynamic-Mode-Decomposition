clear
close all
rng(2); % set random seed to get identical snapshots each time

%% Set parameters
options = odeset('RelTol',1e-15,'AbsTol',1e-15); % for the numerical solver
SIGMA=10;   BETA=8/3;   RHO=28;
ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];

N=1000;                                 % number of delay embeddings
g = @(x,y,z) tanh((x.*y-z*5)/10);       % observable
M = 10^5;                               % number of data points for testing
M2 = 10^4;                              % number of snapshots
dt = 0.05;                              % time step for trajectory sampling

%% Produce the data - slowest part of this script!
Y0=(rand(3,1)-0.5)*4;
[~,Y0]=ode45(ODEFUN,[0.000001 1, 100],Y0,options); Y0 = Y0(end,:)'; % sample after when on the attractor
[~,DATA]=ode45(ODEFUN,[0.000001 (1:((M+(N+1))))*dt],Y0,options);

%% Use delay embedding
PX1=zeros(M,N); PX1(:,1)=DATA(1:M,1); PX2=zeros(M,N); PX2(:,1)=DATA(1:M,2); PX3=zeros(M,N); PX3(:,1)=DATA(1:M,3);
PY1=zeros(M,N); PY1(:,1)=DATA((1:M)+1,1); PY2=zeros(M,N); PY2(:,1)=DATA((1:M)+1,2); PY3=zeros(M,N); PY3(:,1)=DATA((1:M)+1,3);
for j=2:N
    PX1(:,j)=DATA((1:M)+(j-1),1); PX2(:,j)=DATA((1:M)+(j-1),2); PX3(:,j)=DATA((1:M)+(j-1),3);
    PY1(:,j)=DATA((1:M)+1+(j-1),1); PY2(:,j)=DATA((1:M)+1+(j-1),2); PY3(:,j)=DATA((1:M)+1+(j-1),3);
end

%% Run riggedDMD
PX = g(PX1,PX2,PX3);
gmean = mean(PX(:));
PX = PX - gmean;
PY = g(PY1,PY2,PY3)-gmean;

order=6;
epsilon=0.05;
TH = 0.4137;
TH2 = -pi:0.0002:pi;
g_coeffs = zeros(N,1); g_coeffs(1)=1;

[gTH,xi] = riggedDMD(PX(1:M2,:),PY(1:M2,:),1/M2,epsilon,TH,[],'order',order,'g_coeffs',g_coeffs,'TH2',TH2);
xi = xi/sum(xi*(TH2(2)-TH2(1))); % normalise g to ||g||=1

%% Plot the results

figure
plot(TH2,xi,'linewidth',2)
ax=gca; ax.FontSize=18;
grid on
box on

for jj = 1:length(TH)
    C = PX(:,1:end)*gTH(:,jj);
    u = real(C(1:min(10^6,M),:));
    [~,I]=sort(u,'ascend');
    figure
    scatter3(PX1(I,jj),PX2(I,jj),PX3(I,jj),6,u(I),'filled');
    colormap(magma)
    clim([-2*std(u)+mean(u),2*std(u)+mean(u)])
    view(gca,[13.1786087602293 -1.28469255513244]);
    axis tight; grid off; axis off
    axis equal  tight
    set(gca,'DataAspectRatio',[1 1 1]);
    xlabel('$X$','interpreter','latex','fontsize',14)
    ylabel('$Y$','interpreter','latex','fontsize',14)
    zlabel('$Z$','interpreter','latex','fontsize',14,'rotation',0)
    hold off
end

