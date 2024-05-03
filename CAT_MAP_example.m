clear
close all

%% Set up the experiment

F = @(z) mod(2*real(z(:))+imag(z(:)),2*pi) + 1i*mod(real(z(:))+imag(z(:)),2*pi);
N = 500; % number of delay embeddings + 1
M1 = 256; % number of x data points
M2 = 256; % number of x data points

g = chebfun2(@(x,y)   sin((x)) + sin((2*x+y))/2 + 1i*sin((5*x+3*y))/4,[-pi,pi],[-pi,pi],'trig');

%% Collect data and compute the correlations

PX1 = linspace(0,2*pi,M1+1); PX1(1) = [];
PX2 = linspace(0,2*pi,M2+1); PX2(1) = [];
PX = kron(PX1,ones(length(PX2),1))+1i*kron(ones(1,length(PX1)),PX2(:));
PX0 = PX(:);

PX = zeros(length(PX0),N+1);
PX(:,1) = g(real(PX0(:)),imag(PX0(:)));

for j = 1:N
    PX0 = F(PX0);
    PX(:,j+1) = g(real(PX0(:)),imag(PX0(:)));
end

PY = PX(:,2:end);
PX = PX(:,1:end-1);
G = (PX'*PX)/size(PX,1);

COEFFS = G(1,1:3)/G(1,1);
COEFFS = conj(COEFFS(:));
COEFFS = [COEFFS;zeros(N-2,1)];

G=toeplitz(COEFFS,COEFFS');
A = G(1:N,2:N+1);
G = G(1:N,1:N);

%% Perform riggedDMD
TH = 1;                                 % angle for generalised eigenfunction
TH2 = -pi:0.01:pi;                      % angles for spectral measure
epsilon = 0.5;                          % smoothing parameter
order = 2;                              % order of kernel
g_coeffs = zeros(N,1); g_coeffs(1)=1;   % cofficients of g in dictionary expansion

[gTH,xi] = riggedDMD([],[],[],epsilon,TH,[],'QR',0,'G',G,'A',A,'order',order,'g_coeffs',g_coeffs,'TH2',TH2);

%% Plot the results

figure
plot(TH2,xi,'k','linewidth',2)
xlim([-pi,pi])
ax=gca; ax.FontSize=12;
xlabel('$\theta$','interpreter','latex','fontsize',18)
title('Spectral Measure','interpreter','latex','fontsize',18)

C1 = PX(:,1:50)*gTH(1:50,:);
    
figure % increase number of data points for better resolution of plot if wanted
imagesc([-pi,pi],[-pi,pi],reshape(real(C1),M2,M1))
set(gca,'YDir','normal')
colormap(brighten(brewermap([],'RdYlBu'),0))
colorbar
clim([-3*std(real(C1))+mean(real(C1)),3*std(real(C1))+mean(real(C1))])
xlabel('$x$','interpreter','latex','fontsize',14)
ylabel('$y$','interpreter','latex','fontsize',14)
title(T2,'interpreter','latex','fontsize',22)
axis equal tight