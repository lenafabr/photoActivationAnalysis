function [concvals, eiglist, eigfunc] = cylDiffInnerFixOuterRef(rin,tin,options)
% calculate diffusive spread in cylinder with fixed concentration at inner
% bound, reflecting outer bound
% returns concvals(tc,rc) has the value for time tin(tc), radius rin(rc)

opt.a = 1; % inner bound
opt.b = 10; % outer bound
opt.D = 1; % diffusion constant
opt.c0 = 1; % concentration at inner bound
opt.nmax = 500; % number of eigenvalues to use

if (exist('options','var'))
    opt = copyStruct(options,opt);
end

%% ------
% nondimensionalize r and t
rvals = rin/opt.a;
tvals = tin*opt.D/opt.a^2;

if (size(tvals,2)>1) 
    error('tin must be a column vector')
end
if (size(rvals,1)>1)
    error('rin must be a row vector')
end

%% calculate cumulative distribution of hitting times
% between an outer reflecting and inner absorbing cylinder


%% get eigenvalues
[eiglist,eigfunc] = getEigsCyl(opt.a,opt.b,opt.nmax,[0,1],0.1);
eiglist = eiglist(eiglist>1e-8);


%% 
b = opt.b/opt.a;

J0a = besselj(0,eiglist)';
J1a = besselj(1,eiglist)';
Y1a = bessely(1,eiglist)';
J1b = besselj(1,eiglist*b)';
Y1b = bessely(1,eiglist*b)';

betarlist = eiglist'*rvals; % nmax x nr matrix

% nmax by nr matrix of coefficients 
coeff = (J0a.^2.*eiglist'./(J0a.^2-J1b.^2)).*...
    (besselj(0,betarlist).*Y1b - bessely(0,betarlist).*J1b).*...
    (Y1a.*J1b - J1a.*Y1b)*pi^2/2;
%
% evaluate at different times
expvals = exp(-tvals*eiglist.^2);

% Gvals is a ntime by nr matrix of values for conc difference from steady
% state
Gvals = expvals*coeff;

% Fvals is the actual concentration
concvals = 1 - Gvals;

concvals = concvals*opt.c0;
