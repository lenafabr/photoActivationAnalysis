function [concvals, dtconc, dt2conc] = cylDiffInnerFix(rin,tin,options)
% calculate diffusive spread in cylinder with fixed concentration at inner
% bound, no outer bound
% NOTE: umin, minmaxu should be provided as *dimensionless* quantities
% rin, tin are provided with real dimensions

opt.a = 1; % inner bound
opt.D = 1; % diffusion constant
opt.c0 = 1; % concentration at inner bound
opt.nu1 = 100; % number of u for low limit numerical integration (u<1)
opt.nu2 = 5000; % number of u for high limit numerical integration
opt.minmaxu = 100; % go up at least this high in u
opt.maxuscl = 2; % how high up to go in u, relative to exponential cutoff
opt.umin = 1e-6; % minimal u for numerical integration, use asymptotic form below this
opt.mint = 1e-2; % below this cutoff in time, use the low-time asymptic approximation


if (exist('options','var'))
    opt = copyStruct(options,opt);
end

%% ------
% nondimensionalize r and t
rvals = rin/opt.a;
tlist = tin*opt.D/opt.a^2;

concvals = zeros(length(tlist),length(rvals));
dtconc = concvals;
dt2conc = concvals;

%% ------
% nondimensionalize r and t
indhigh = find(tlist>opt.mint*min(rvals).^2);
tvals = tlist(indhigh);

if (size(tvals,2)>1)
    error('tin must be a column vector')
end

% euler gamma
gam = double(eulergamma);

% u values over which to do numerical integral
umax = max(opt.maxuscl*sqrt(-log(eps))/(min(tvals)),opt.minmaxu);
ulist = [logspace(log10(opt.umin),0,opt.nu1) linspace(1,umax,opt.nu2)]';
du = diff(ulist);

uR = ulist*rvals;
denom = (besselj(0,ulist).^2 + bessely(0,ulist).^2).*ulist;
integrand = (besselj(0,uR).*bessely(0,ulist) - bessely(0,uR).*besselj(0,ulist))./denom;

% exponential prefactor involving time and u
expvals = exp(-tvals*ulist'.^2);

% do the integral
intnum1 = expvals(:,1:end-1)*(integrand(1:end-1,:).*du);
intnum2 = expvals(:,2:end)*(integrand(2:end,:).*du);
intnum = 0.5*(intnum1+intnum2);

% get the small u integral
intsmall =  - (atan(2/pi*(gam + log(opt.umin/2))) + pi/2)*log(rvals);
intsmall = ones(length(tvals),length(rvals)).*intsmall;

concvals(indhigh,:) = opt.c0*(1 + 2/pi*(intnum+intsmall));

%% Now get first derivative wrt time
if (nargout>1)
    integrandDt = -integrand.*ulist.^2;
    intnum1 = expvals(:,1:end-1)*(integrandDt(1:end-1,:).*du);
    intnum2 = expvals(:,2:end)*(integrandDt(2:end,:).*du);
    intnum = 0.5*(intnum1+intnum2);
    
    dtconc(indhigh,:) = opt.c0*2/pi*intnum*opt.D/opt.a^2;
end

%% Now get second derivative wrt time
if (nargout>2)
    integrandDt = integrand.*ulist.^4;
    intnum1 = expvals(:,1:end-1)*(integrandDt(1:end-1,:).*du);
    intnum2 = expvals(:,2:end)*(integrandDt(2:end,:).*du);
    intnum = 0.5*(intnum1+intnum2);
    
    dt2conc(indhigh,:) = opt.c0*2/pi*intnum*(opt.D/opt.a^2)^2;
end

% ------------------------
%% do the low t asymptotics separately for each r
% ------------------------
for rc = 1:length(rvals)
    r = rvals(rc);
    indlow = find(tlist<opt.mint*r.^2);
    tvals = tlist(indlow);
    
    ierfc = @(x) -x.*erfc(x) + exp(-x.^2)/sqrt(pi);
    i2erfc = @(x) 0.25*(erfc(x) - 2*x.*ierfc(x));
    
    % first derivatives
    derfc = @(x) -2*exp(-x.^2)/sqrt(pi);
    d2erfc = @(x) 4*x.*exp(-x.^2)/sqrt(pi);
    
    rt = (1./sqrt(tvals))*(r-1)/2;
    
    % r-dependent prefactors
    rterm1 =  1./sqrt(r);
    rterm2 = ((r-1)./(4*r.^1.5));
    rterm3 = (9 - 2*r- 7*r.^2)./(32*r.^2.5);
    erfcvals = erfc(rt);
    ierfcvals = ierfc(rt);
    i2erfcvals = i2erfc(rt);
    
    concvals(indlow,rc) = rterm1.*erfcvals + sqrt(tvals).*ierfcvals.*rterm2...
        + tvals.*i2erfcvals.*rterm3;
    
    concvals(indlow,rc) = concvals(indlow,rc)*opt.c0;
    
    % first time derivative
    if (nargout>1)
        derfcvals = derfc(rt);
        dtconc(indlow,rc) = (rterm1.*derfcvals + sqrt(tvals).*erfcvals.*rterm2...
            + tvals.*ierfcvals.*rterm3).*rt./(-2*tvals)...
            +0.5./sqrt(tvals).*ierfcvals.*rterm2 + i2erfcvals.*rterm3;
    end
    
    % second time derivative
    if (nargout>2)
        dt2conc(indlow,rc) = (rterm1.*d2erfc(rt) + sqrt(tvals).*derfcvals.*rterm2...
            + tvals.*erfcvals.*rterm3).*rt.^2./(2*tvals).^2 ...
            + (rterm1.*derfcvals + sqrt(tvals).*erfcvals.*rterm2...
            + tvals.*ierfcvals.*rterm3).*(rt*0.75./tvals.^2)...
            +(0.5./sqrt(tvals).*erfcvals.*rterm2 + ierfcvals.*rterm3).*rt./(-2*tvals)...
            -0.25./tvals.^1.5.*ierfcvals.*rterm2;
    end
    
end
end