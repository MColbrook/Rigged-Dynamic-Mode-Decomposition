clear
close all

%% Set up the experiment

N = 300; % number of delay embeddings + 1
M1 = 500; % number of x data points
M2 = 500; % number of y data points

g = @(x,y) exp(1i*x)./cosh(y);

%% Collect the data
[PX,PY,W,x1,x2] = pendulum_matrix(M1,M2,N,4,g);

%% Perform riggedDMD
TH = pi/5;                              % angle for generalised eigenfunction
TH2 = -pi:0.01:pi;                      % angles for spectral measure
epsilon = 0.05;                         % smoothing parameter
order = 6;                              % order of kernel
g_coeffs = zeros(N,1); g_coeffs(1)=1;   % cofficients of g in dictionary expansion

[gTH,xi] = riggedDMD(PX(:,1:N),PY(:,1:N),W,epsilon,TH,[],'order',order,'g_coeffs',g_coeffs,'TH2',TH2);
xi = (max(xi,0.000001)/sum(xi*(TH2(2)-TH2(1))));

%% Plot the results
figure
plot(TH2/pi,xi,'linewidth',2)
xlim([0,1])
ylim([0,2])
ax=gca; ax.FontSize=18;
grid on
box on

C = PX(:,1:200)*gTH(1:200,1);

figure
toPlot = real((reshape(C,length(x2),length(x1))));
contourf(x1, x2 , toPlot ,60,'edgecolor','none');
colormap(brighten(brewermap([],'RdYlBu'),0))
axis equal on;   view(0,90);    ylim([-4,4]);   xlim([-pi,pi])
clim([-max(abs(toPlot(:)))*0.8,max(abs(toPlot(:)))*0.8])
ax=gca; ax.FontSize=16;
set(gca,'YDir','normal')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PSI_X,PSI_Y,W,x1,x2] = pendulum_matrix(M1,M2,N,L,g)
% computes data matrices for the pendulum
delta_t=1;
x1=linspace(-pi,pi,M1+1);   x1=x1+(x1(2)-x1(1))/2; x1(end)=[];
x2=linspace(-L,L,M2+1); x2=x2+(x2(2)-x2(1))/2; x2(end)=[];
[X1,X2] = meshgrid(x1,x2);  X1=X1(:); X2=X2(:);
M=length(X1); % number of data points

Y0=[X1(:)';X2(:)'];
[~,Y]=ode45(@(t,y) pensystem(y,length(X1)),[0.000001 (1:N)*delta_t],Y0);

XX = mod([X1(:)';Y(2:end,1:2:end)]+pi,2*pi)-pi;
YY = [X2(:)';Y(2:end,2:2:end)];

Y = transpose(g(XX,YY));
PSI_X = Y(:,1:end-1);
PSI_Y = Y(:,2:end);

W=zeros(M,1)+(x1(2)-x1(1))*(x2(2)-x2(1));
end



function dydt = pensystem(y,n)
y = reshape(y,[],n);
dydt = [y(2,:);-sin(y(1,:))];

% Linearize output.
dydt = dydt(:);
end



