function results = fmgaussfit(xx,yy,zz)
%, zfit, fiterr, zerr, resnorm, rr]
% FMGAUSSFIT Create/alter optimization OPTIONS structure.
%   [fitresult,..., rr] = fmgaussfit(xx,yy,zz) uses ZZ for the surface 
%   height. XX and YY are vectors or matrices defining the x and y 
%   components of a surface. If XX and YY are vectors, length(XX) = n and 
%   length(YY) = m, where [m,n] = size(Z). In this case, the vertices of the
%   surface faces are (XX(j), YY(i), ZZ(i,j)) triples. To create XX and YY 
%   matrices for arbitrary domains, use the meshgrid function. FMGAUSSFIT
%   uses the lsqcurvefit tool, and the OPTIMZATION TOOLBOX. The initial
%   guess for the gaussian is places at the maxima in the ZZ plane. The fit
%   is restricted to be in the span of XX and YY.
%   See:
%       http://en.wikipedia.org/wiki/Gaussian_function
%          
%   Examples:
%     To fit a 2D gaussian:
%       [fitresult, zfit, fiterr, zerr, resnorm, rr] =
%       fmgaussfit(xx,yy,zz);
%   See also SURF, OMPTMSET, LSQCURVEFIT, NLPARCI, NLPREDCI.

%   Copyright 2013, Nathan Orloff.

%% Condition the data
[xData, yData, zData] = prepareSurfaceData( xx, yy, zz );
xyData = {xData,yData};
%% Set up the startpoint
[amp, ind] = max(zData); % amp is the amplitude.
xo = xData(ind); % guess that it is at the maximum
yo = yData(ind); % guess that it is at the maximum
ang = 0; % angle in radians.
sx = 0.5;
sy = 0.5;
% zo = median(zData(:))-std(zData(:));
xmax = max(xData)+0.5;
ymax = max(yData)+0.5;
xmin = min(xData)-0.5;
ymin = min(yData)-0.5;

%% Set up fittype and options.
% Lower = [0, 0.0001, 0, 0, xmin, ymin, 0];
% Upper = [255,pi/2+0.0001, 3, 3, xmax, ymax, Inf]; % angles greater than 90 are redundant
% StartPoint = [amp, ang, sx, sy, xo, yo, zo];%[amp, sx, sxy, sy, xo, yo, zo];

Lower = [0, 0.0001, 0, 0, xmin, ymin];
Upper = [255,pi/2+0.0001, 3, 3, xmax, ymax]; % angles greater than 90 are redundant
StartPoint = [amp, ang, sx, sy, xo, yo];%[amp, sx, sxy, sy, xo, yo, zo];

tols = 1e-16;
% options = optimset('Display','off',...
%     'MaxFunEvals',5e2,...
%     'MaxIter',5e2,...
%     'TolX',tols,...
%     'TolFun',tols,...
%     'TolCon',tols );
options = optimset('Display','off','TolFun',tols,'LargeScale','off');
%% perform the fitting
results = lsqcurvefit(@gaussian2D,StartPoint,xyData,zData,Lower,Upper,options);
% b = 1/(2*fitresult(3)^2); c = 1/(2*fitresult(4)^2); 
% results = [fitresult(1), fitresult(2), b, c, fitresult(5),fitresult(6)];


end


function z = gaussian2D(par,xy)
% compute 2D gaussian
xx=(xy{1}-par(5)).*cos(par(2))+(xy{2}-par(6)).*sin(par(2));
yy=-(xy{1}-par(5)).*sin(par(2))+(xy{2}-par(6)).*cos(par(2));
% z = par(7) + ...
%     par(1)*exp(-(xx.^2/(2*par(3)^2)+yy.^2./(2*par(4)^2)));
z = par(1)*exp(-(xx.^2 * par(3) + yy.^2 * par(4)));
end

