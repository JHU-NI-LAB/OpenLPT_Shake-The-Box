function OTFParam = OTF(data, image_dir, calib_dir, num_grid, x_range,y_range,z_range, particle_radius, num_calibpart)
camParaCalib = load(calib_dir);
camParaCalib = camParaCalib.camParaCalib;

dx = (x_range(2) - x_range(1)) / num_grid; 
dy = (y_range(2) - y_range(1)) / num_grid;
dz = (z_range(2) - z_range(1)) / num_grid;

ncams = 4;
mkdir([image_dir, 'OTFfiles']);
a = zeros(ncams,num_grid, num_grid, num_grid);
b = zeros(ncams,num_grid, num_grid, num_grid);
c = zeros(ncams,num_grid, num_grid, num_grid);
alpha = zeros(ncams,num_grid, num_grid, num_grid);

for i = 1 : num_grid
    for j = 1 : num_grid
        for k = 1 : num_grid
            % get the range of the sub grid.
            box_range_x = [x_range(1) + (i - 1) * dx, x_range(1) + i * dx];
            box_range_y = [y_range(1) + (j - 1) * dy, y_range(1) + j * dy];
            box_range_z = [z_range(1) + (k - 1) * dz, z_range(1) + k * dz];
            % get particles inside this sub grid
            particle_inside_box = data(data(:,1) > box_range_x(1) & data(:,1) < box_range_x(2) & ...
                data(:,2) > box_range_y(1) & data(:,2) < box_range_y(2) & ...
                data(:,3) > box_range_z(1) & data(:,3) < box_range_z(2), :);

            num_particle = size(particle_inside_box, 1);
            if num_particle > num_calibpart
                particle_inside_box = particle_inside_box(randperm(num_particle, num_calibpart), :); % randomly choose particle
                num_particle = num_calibpart;
            end
            OTFParam_array = zeros(ncams, num_particle, 6);
            for l = 1 : num_particle
                for m = 1 : ncams
                    % read the image
                    nthframe = particle_inside_box(l, 4);
                    img = imread([image_dir 'cam' num2str(m) '/cam' num2str(m) 'frame' num2str(nthframe,'%05d') '.tif']);
                    [ysize,xsize]=size(img);
                    % check the particle is overlapping or not
                    % get the original projection of particle on image
                    center_2D = round(calibProj(camParaCalib(m), particle_inside_box(l, 1:3)));
                    if sum(center_2D(1) > xsize - particle_radius | center_2D(1) < particle_radius + 1 | ...
                            center_2D(2) > ysize - particle_radius | center_2D(2) < particle_radius + 1) 
                        continue; %skip this particle because it's out of image
                    end
                    image_range_x = center_2D(1) - particle_radius : center_2D(1) + particle_radius;
                    image_range_y = center_2D(2) - particle_radius : center_2D(2) + particle_radius;
                    
                    particle_projection = double(img(image_range_y, image_range_x));
                    OTFParam_array(m, l, :)= fmgaussfit(image_range_x,image_range_y,particle_projection);
                end
            end 
            save([image_dir, 'OTFfiles/Grid_' num2str(i) '_' num2str(j) '_' num2str(k) '.mat'], 'OTFParam_array');
            for m = 1 : ncams
                if ~isempty(nonzeros(OTFParam_array(m, :, :))) 
                    a(m, i, j , k) = mean(nonzeros(OTFParam_array(m, :, 1)));
                    alpha(m, i, j , k) = mean(nonzeros(OTFParam_array(m, :, 2)));
                    b(m, i, j , k) = mean(nonzeros(OTFParam_array(m, :, 3)));
                    c(m, i, j , k) = mean(nonzeros(OTFParam_array(m, :, 4)));
                end
            end
        end
    end
end
OTFParam = {a; b; c; alpha};
end

function Xtest_proj = calibProj(camParaCalib, Xtest3D)
% Use the calibrated camera parameters to predict the particle position
% projected onto the image plane.
%
% inputs:
%   camParaCalib    --  calibrated camera parameters
%   Xtest3D         --  test particle coordinates in 3D world system
%
% output:
%   Xtest_proj      --  projected particle position on image plane (in
%   pixels)

    Xc = Xtest3D * (camParaCalib.R)';
    Xc(:,1) = Xc(:,1) + camParaCalib.T(1);
    Xc(:,2) = Xc(:,2) + camParaCalib.T(2);
    Xc(:,3) = Xc(:,3) + camParaCalib.T(3);
    dummy = camParaCalib.f_eff ./ Xc(:,3);
    Xu = Xc(:,1) .* dummy;  % undistorted image coordinates
    Yu = Xc(:,2) .* dummy;
    ru2 = Xu .* Xu + Yu .* Yu;
    dummy = 1 + camParaCalib.k1 * ru2;
    Xd = Xu ./ dummy;
    Yd = Yu ./ dummy;
    
    Np = size(Xtest3D,1);
    Xtest_proj = zeros(Np,2);
    Xtest_proj(:,1) = Xd/camParaCalib.wpix + camParaCalib.Noffw + camParaCalib.Npixw/2;
    Xtest_proj(:,2) = camParaCalib.Npixh/2 - camParaCalib.Noffh - Yd/camParaCalib.hpix;
end

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
% results = lsqcurvefit(@GaussianIntegration,StartPoint,xyData,zData,Lower,Upper,options);
% b = 1/(2*fitresult(3)^2); c = 1/(2*fitresult(4)^2); 
% results = [fitresult(1), fitresult(2), b, c, fitresult(5),fitresult(6)];


end

% function z = GaussianIntegration(par, xy)
% z = zeros(size(xy{1}, 1), 1);
% for i = 1 : size(xy{1}, 1)
%     z(i) = integral2(@(x, y)gaussian2D(par, x, y), xy{1}(i) - .5, xy{1}(i) + .5, xy{2}(i) -.5, xy{2}(i) + .5);
% end
% end

function z = gaussian2D(par,xy)
% compute 2D gaussian
xx=(xy{1}-par(5)).*cos(par(2))+(xy{2}-par(6)).*sin(par(2));
yy=-(xy{1}-par(5)).*sin(par(2))+(xy{2}-par(6)).*cos(par(2));
% z = par(7) + ...
%     par(1)*exp(-(xx.^2/(2*par(3)^2)+yy.^2./(2*par(4)^2)));
z = par(1)*exp(-(xx.^2 * par(3) + yy.^2 * par(4)));
end

% function z = gaussian2D(par,x, y)
% % compute 2D gaussian
% xx=(x-par(5)).*cos(par(2))+(y-par(6)).*sin(par(2));
% yy=-(x-par(5)).*sin(par(2))+(y-par(6)).*cos(par(2));
% % z = par(7) + ...
% %     par(1)*exp(-(xx.^2/(2*par(3)^2)+yy.^2./(2*par(4)^2)));
% z = par(1)*exp(-(xx.^2 * par(3) + yy.^2 * par(4)));
% end

