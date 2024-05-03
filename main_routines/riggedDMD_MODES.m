function [gMODES] = riggedDMD_MODES(PX,PY,W,epsilon,TH,g,varargin)
% Same as riggedDMD but now we compute the modes (end algorithm of paper)
% for a matrix of g_coeffs
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

gMODES = zeros(size(mpV,1),size(c,2),length(TH));

for jj = 1:length(TH)
    th = TH(jj);
    gTH = zeros(size(mpV,1),size(c,2));
    for mm=1:p.Results.order
        EXP = exp(1i*th-1i*epsilon*poles(mm));
        gp = c.*(mpE(:) + transpose(EXP))./(mpE(:) - transpose(EXP));
        EXP = exp(1i*th-1i*epsilon*conj(poles(mm)));
        gm = c.*(mpE(:) + transpose(EXP))./(mpE(:) - transpose(EXP));
        gTH = gTH + res(mm)*gp - conj(res(mm))*gm;
    end
    gMODES(:,:,jj) = -gTH/(4*pi);
end



end