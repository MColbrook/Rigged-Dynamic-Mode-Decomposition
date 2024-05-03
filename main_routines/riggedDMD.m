function [gTH,xi] = riggedDMD(PX,PY,W,epsilon,TH,g,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:                   PX: dictionary evaluated at snapshots,
%                           PY: dictionary evaluated at snapshots one time step later
%                           W:  vector of weights for the quadrature
%                           epsilon: smoothing parameter
%                           TH: angles at which we evaluate the generalized eigenfunction
%                           (wave packet approximation)
%                           g: observable evaluated at data points
%
% OPTIONAL LABELLED INPUTS: order of kernel used, default is 2
%                           PoleType: where to place poles, default is equispaced
%                           TH2: angles at which we evluate the spectral
%                           measures
%                           QR: tells us whether to use QR version of
%                           mpEDMD (1 or 0, default is 1). If
%                           0 the matrices G and A must be provided
%                           g_coeffs: coefficients of g in the dictionary
%
% OUTPUTS:   Koopman matrix mpK and its eigendecomposition [mpV,mpD]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLEASE CITE:  1.  Colbrook, Matthew J., Catherine Drysdale, and Andrew Horning.
%                   "Rigged Dynamic Mode Decomposition: Data-Driven Generalized Eigenfunction Decompositions for Koopman Operators."
%               2.  Colbrook, Matthew J. 
%                   "The mpEDMD algorithm for data-driven computations of measure-preserving dynamical systems." 
%                   SIAM Journal on Numerical Analysis 61.3 (2023): 1585-1608.
%               3.  Colbrook, Matthew J., and Alex Townsend.
%                   "Rigorous data‐driven computation of spectral properties of Koopman operators for dynamical systems."
%                   Communications on Pure and Applied Mathematics 77.1 (2024): 221-283.
% Author and copyright:  Matthew Colbrook, https://www.damtp.cam.ac.uk/user/mjc249/home.html
% Comments and questions: m.colbrook@damtp.cam.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
p.CaseSensitive = true;
addRequired(p,'PX',@isnumeric);
addRequired(p,'PY',@isnumeric);
addRequired(p,'W',@isnumeric);
addRequired(p,'epsilon',@isnumeric);
addRequired(p,'TH',@isnumeric);
addRequired(p,'g',@isnumeric);

validPole = {'cheb','roots','extrap','equi'};
checkPole = @(x) any(validatestring(x,validPole));

addParameter(p,'order',2,@(x) x==floor(x))
addParameter(p,'PoleType','equi',checkPole)
addParameter(p,'TH2',[],@isnumeric)
addParameter(p,'QR',1,@(x) (x==1)+(x==0) )
addParameter(p,'G',[],@isnumeric)
addParameter(p,'A',[],@isnumeric)
addParameter(p,'g_coeffs',[],@isnumeric)

parse(p,PX,PY,W,epsilon,TH,g,varargin{:});

%% Stage A: apply mpEDMD

if p.Results.QR == 1
    [~,mpV,mpD] = mpEDMDqr(PX,PY,W);
else
    [~,mpV,mpD] = mpEDMD(p.Results.G,p.Results.A);
end
mpE = diag(mpD); clear mpD

[~,I]=sort(abs(angle(mpE)),'ascend'); % sort eigenvalues by closest to lambda = 1
mpV = mpV(:,I);
mpE = mpE(I);

if isempty(p.Results.g_coeffs)
    c = (sqrt(W)*PX*mpV)\(sqrt(W)*g);
else
    c = mpV\p.Results.g_coeffs;
end

clear PX PY W

%% Stage B: apply smoothing kernels

[poles,res]=rational_kernel(p.Results.order,p.Results.PoleType);

% compute generalized eigenfunctions

gTH = zeros(size(mpV,1),length(TH));


for mm=1:p.Results.order
    EXP = exp(1i*TH(:)-1i*epsilon*poles(mm));
    gp = c.*(mpE(:) + transpose(EXP))./(mpE(:) - transpose(EXP));
    EXP = exp(1i*TH(:)-1i*epsilon*conj(poles(mm)));
    gm = c.*(mpE(:) + transpose(EXP))./(mpE(:) - transpose(EXP));

    gTH = gTH + res(mm)*gp - conj(res(mm))*gm;
end
gTH = -mpV*gTH/(4*pi);

% compute spectral measures

if isempty(p.Results.TH2)
    xi = [];
else
    xi = zeros(1,length(p.Results.TH2));

    for mm=1:p.Results.order
        EXP = exp(1i*p.Results.TH2(:)-1i*epsilon*poles(mm));
        gp = c.*(mpE(:) + transpose(EXP))./(mpE(:) - transpose(EXP));
    
        xi = xi - real(res(mm)*(c'*gp))/pi/2;
    end
end



end