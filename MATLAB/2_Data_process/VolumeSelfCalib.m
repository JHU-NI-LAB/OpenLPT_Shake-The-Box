function [camParaCalib, particles_info, good_particles, good_tracks] = VolumeSelfCalib(tracks, view_size, image_path, start_frame, calib_path, save_name_label, load_file)
if ~exist('load_file', 'var')
    % Pick good tracks
     good_tracks = PickGoodTracks(tracks, view_size);

    % Choose 4000 particles evenly 
    good_particles = PickParticleEvenly(good_tracks, 8000);

    %start_frame = 0;
    good_particles(:, 4) = good_particles(:, 4) + start_frame - 1;

    % Get the 2D positions and 3D poisitions for each particle
    warning('off','all');
    particles_info = GetParticleInfo(good_particles, image_path, calib_path);
    warning('on','all');
    save([image_path 'VSCdata_' save_name_label '.mat'], 'good_tracks', 'good_particles', 'particles_info');
else
   load([image_path 'VSCdata_' save_name_label '.mat']);
end

%     fig = figure;
%     PlotTracks(good_tracks, fig, 'g.');
%     hold on
%     plot3(good_particles(:, 1), good_particles(:, 2), good_particles(:, 3), 'b.');
%     plot3(particles_info(:, 1), particles_info(:, 2), particles_info(:, 3), 'r*');
%     pause(1);
% Do the VSC
camParaCalib = VSC(particles_info(:, 1:11), calib_path);

% output the calibration file
CalibMatToTXT(camParaCalib, [image_path 'VSC_Calib_' save_name_label '.txt']);
save([image_path 'VSC_Calib_' save_name_label '.mat'], 'camParaCalib');
end

function good_particles = PickParticleEvenly(tracks, num_particles)
% the number of segments on each axis
num_seg = 10;
num_particle_seg = ceil(num_particles / num_seg ^ 3); % number of cubes
x_min = min(min(tracks(:, 3))); x_max = max(max(tracks(:, 3)));
y_min = min(min(tracks(:, 4))); y_max = max(max(tracks(:, 4)));
z_min = min(min(tracks(:, 5))); z_max = max(max(tracks(:, 5)));
select_label = zeros(size(tracks, 1), 1);

for i = 1 : num_seg
    for j = 1 :  num_seg
        for k = 1 : num_seg
            x_lo = x_min + (i - 1) * (x_max - x_min) / num_seg;
            y_lo = y_min + (j - 1) * (y_max - y_min) / num_seg;
            z_lo = z_min + (k - 1) * (z_max - z_min) / num_seg;
            x_up = x_min + i * (x_max - x_min) / num_seg;
            y_up = y_min + j * (y_max - y_min) / num_seg;
            z_up = z_min + k * (z_max - z_min) / num_seg;
            index = find(tracks(:, 3) > x_lo & tracks(:, 3) < x_up & ...
                tracks(:, 4) > y_lo & tracks(:, 4) < y_up & ...
                tracks(:, 5) > z_lo & tracks(:, 5) < z_up);
            total_num = length(index);
            if total_num > num_particle_seg
                % if the number of particles in this cube is more than the
                % number to be selected, then choose them randomly
                ind = 1 : total_num;
                index_selected = datasample(ind, num_particle_seg, 'Replace',false);
                select_label(index(index_selected)) = 1;
                    % indicate this particle is selected
            elseif total_num > 0
                % if the number of particles in this cube is less than the
                % number to be selected, then select all of them
                select_label(index) = 1;
                    % indicate this particle is selected
            end
        end
    end
end

% for the number of selected particle is less than the required number,
% then select them again from the orginal data randomly
num_left = num_particles - sum(select_label);
while num_left > 0
    select_index = datasample(1 : size(select_label, 1), num_left);
    for i = 1 : num_left
        if select_label(select_index(i)) == 0
           select_label(select_index(i)) = 1;
           num_left = num_left - 1;
        end
    end
end

% fill in the good_particle matrix
index_tracks = select_label == 1;
good_particles = tracks(index_tracks, :);
good_particles = [good_particles(:, 3 : 5), good_particles(:, 2)];
end

function particles_info = GetParticleInfo(particles, image_path, calib_path)
num_particles = size(particles, 1);
load(calib_path);
particles_info = zeros(num_particles, 15);

ii = 1; % index for particles_info
for i = 1 : num_particles
% Project the particle on each image
   xc = zeros(1, 4); yc = zeros(1, 4);
   for j = 1 : 4
       X2D = calibProj(camParaCalib(j), particles(i, 1:3));
       xc(j) = X2D(1); yc(j) = X2D(2);
   end
   
   delete_particle  = 0;
% Get the 2D positions around the projected point
    for j = 1 : 4
        img = imread([image_path 'cam' num2str(j) '/cam' num2str(j) 'frame' num2str(particles(i, 4), '%05.0f') '.tif']);
        [ymax, xmax] = size(img);
        particle_size = 4;
        search_range = particle_size * 3 / 4; % searching range on the image
        img_searcharea = img(max(1, floor(yc(j) - search_range)) : min(ymax, ceil(yc(j) + search_range)), ...
            max(1, floor(xc(j) - search_range)) : min(xmax, ceil(xc(j) + search_range))); % get the search area
        position2D_candidate = Get2DPosOnImage(img_searcharea);
        if ~isempty(position2D_candidate)
            position2D_candidate(:, 1) = position2D_candidate(:,1) + max(1, floor(xc(j) - search_range)) - 1;
            position2D_candidate(:, 2) = position2D_candidate(:,2) + max(1, floor(yc(j) - search_range))- 1;
            switch j
                case 1
                    position2D_candidate_1 = position2D_candidate;
                case 2
                    position2D_candidate_2 = position2D_candidate;
                case 3
                    position2D_candidate_3 = position2D_candidate;
                case 4
                    position2D_candidate_4 = position2D_candidate;
            end
        else
            delete_particle  = 1;
            break;
        end
    end
    
    if delete_particle % delete the particle 
        particles_info(ii, :) = [];
        continue;
    end

% Do triangulation for each combination of 2D positions
    len1 = size(position2D_candidate_1, 1); len2 = size(position2D_candidate_2, 1);
    len3 = size(position2D_candidate_3, 1); len4 = size(position2D_candidate_4, 1);
    particle = zeros(len1 * len2 * len3 * len4, 15);
    error = zeros(1, len1 * len2 * len3 * len4);
    for j = 1 : len1
        for k = 1 : len2
            for m = 1 : len3
                for n = 1 : len4
                    position2D = [position2D_candidate_1(j, :), position2D_candidate_2(k, :), ...
                        position2D_candidate_3(m, :), position2D_candidate_4(n, :)];
                    [position3D, e] = Triangulation(camParaCalib, position2D);
                    error((j - 1) * len2 * len3 * len4 + (k - 1) * len3 * len4 + (m -1) * len4 + n) = e;
                    particle((j - 1) * len2 * len3 * len4 + (k - 1) * len3 * len4 + (m -1) * len4 + n, :) = ...
                        [position3D', position2D, norm(position3D' - particles(i, 1:3)), particles(i, 1:3)];
                    % the last one is the distance from the projected particle
                end
            end
        end
    end
    
% Choose the closest one to the projected particle as the real particle
    thred_3D = 2.5; 
    ind = find(error < thred_3D);
    if isempty(ind)
        particles_info(ii, :) = []; % delete the particle since there is no suitable candidate
        continue;
    end
    
    particle = particle(ind, :);
    error = error(ind);
    [~, index] = min(particle(:, 12) / mean(particle(:, 12)) ...
        + error' / mean(error)); % combination of goal achieving the minimum distance as well as minimum error
%     [~, index] = min(error);
%     [~, index] = min(particle(:, 12));
    particles_info(ii, :) = particle(index, :);
    ii = ii + 1;
    if ~(mod(i, 100)) 
        ['Processing the ' num2str(i) 'th particle...']
    end
end

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

function position2D = Get2DPosOnImage(img) 
[rows, cols] = size(img);
threshold = 35;
position2D = [];
for i = 2 : rows - 1
    for j = 2 : cols - 1
        if (img(i, j) >= threshold && IsLocalMax(img, i, j))
            x1 = j - 1; x2 = j; x3 = j + 1;
            y1 = i - 1; y2 = i; y3 = i + 1;
            ln1 = NoInfLog(img(i, j - 1));
            ln2 = NoInfLog(img(i, j));
            ln3 = NoInfLog(img(i, j + 1));
        
            xc = -.5 * (ln1 * (x2 ^ 2 - x3 ^ 2) - ln2 * (x1 ^ 2 - x3 ^ 2) + ln3 * (x1 ^ 2 - x2 ^ 2)) / ...
            (ln1 * (x3 - x2) - ln3 * (x1 - x2) + ln2 * (x1 - x3));
            
            ln1 = NoInfLog(img(i - 1, j));
            ln2 = NoInfLog(img(i, j));
            ln3 = NoInfLog(img(i + 1, j));
            yc = -.5 * (ln1 * (y2 ^ 2 - y3 ^ 2) - ln2 * (y1 ^ 2 - y3 ^ 2) + ln3 * (y1 ^ 2 - y2 ^ 2)) / ...
            (ln1 * (y3 - y2) - ln3 * (y1 - y2) + ln2 * (y1 - y3));
            
            if (~isinf(xc) && ~isinf(yc))
                position2D = [position2D; xc, yc];
            end
        end
    end
end

end

function result = IsLocalMax(img, i, j)
   result = 1;
   if (img(i - 1, j) > img(i, j) || img(i + 1, j) > img(i, j) || ...
           img(i, j - 1) > img(i, j) || img(i, j + 1) > img(i, j))
       result = 0;
   end
end

function y = NoInfLog(x)
if x == 0
    x = .0001;
end
y = log(double(x));
end

function [position3D, error] = Triangulation(camParaCalib, position2D) 
A = zeros(3, 3);
B = zeros(3, 1);
D = 0;
for i = 1 : 4
    % get the world position of the point on the image
    SIpos = Img2World(camParaCalib(i), UnDistort(position2D((i - 1) * 2 + 1 : (i - 1) * 2 + 2), camParaCalib(i))); 
    % get the vector 
    sight = SIpos - camParaCalib(i).Tinv;
    sight = sight / norm(sight); 
    C = eye(3) - sight *  sight';
    A = A + C;
    B = B + C * camParaCalib(i).Tinv;
    D = D + camParaCalib(i).Tinv' * C * camParaCalib(i).Tinv;
end
position3D = A \ B;
error = (position3D' * A * position3D - 2 * position3D' * B + D) ^ .5 / 4;
end

function camParaCalib = VSC(particles_info, calib_path)
%% define parameters
ncams = 4;
fixcamera=[]; % which camera is fixed. fixcamera=[1,4]; camera 1 and 4 are fixed. fixcamera=[]; all camera can be changed. 
fmin_options.Display = 'iter'; %'final' to display only the final mean distance, 'iter' to show the steps along the way.
fmin_options.MaxFunEvals = 3000; %350;

load (calib_path);
camParaCalib = camParaCalib([1,2,3,4]);

%% pick data to shake
data = particles_info;
data1 = data(1:end,:);

%% all data to test
cam2d = zeros(size(data1,1),2,ncams);

%% organize calibration file
for icam = 1 : ncams
    cam2d(:, :, icam) = data1(:,(icam - 1) * 2 + 4 : (icam - 1) * 2 + 5);
    cam2d(:, 1, icam) = (cam2d(:, 1, icam) - camParaCalib(icam).Npixw / 2 ...
        - camParaCalib(icam).Noffw) * camParaCalib(icam).wpix;
    cam2d(:, 2, icam) = (-cam2d(:,2,icam) + camParaCalib(icam).Npixh / 2 ...
        - camParaCalib(icam).Noffh) * camParaCalib(icam).hpix;  %vertical coordinate needed to be switched in sign to make it work--I thought the relection in the rotation matrix took care of this--but it doesn't work without this sign negative
end

data2 = data(1 : end, :);

cam2d_check = zeros(size(data2,1), 2, ncams);

for icam = 1 : ncams
    cam2d_check(:,:,icam) = data2(:, (icam-1) * 2 + 4 : (icam-1) * 2 + 5);
    cam2d_check(:,1,icam) = (cam2d_check(:, 1, icam) - camParaCalib(icam).Npixw / 2 ...
        - camParaCalib(icam).Noffw) * camParaCalib(icam).wpix;
    cam2d_check(:,2,icam) = (-cam2d_check(:, 2, icam) + camParaCalib(icam).Npixh / 2 ...
        - camParaCalib(icam).Noffh) * camParaCalib(icam).hpix;  %vertical coordinate needed to be switched in sign to make it work--I thought the relection in the rotation matrix took care of this--but it doesn't work without this sign negative
end

calinitial=zeros(8,ncams);
for icam = 1 : ncams
    eul = rotm2eul(camParaCalib(icam).R);
    calinitial(1:3,icam) = eul';    
    calinitial(4:6,icam) = camParaCalib(icam).T;
    calinitial(7,icam) = camParaCalib(icam).f_eff;
    calinitial(8,icam) = camParaCalib(icam).k1;
end

%%

%calopt=calarray(1:6,2:ncams);
% calopt=calinitial([1:5,8],:);
% fixpara=calinitial(6:7,:);
calopt = calinitial;

%% fix cameras?

calall = calopt;

if isempty(fixcamera)==0
    calfix = calall(:, fixcamera);
    left = setdiff(1: ncams, fixcamera);
    calnofix = calall(:, left);
else
    calfix = [];
    calnofix = calall;
end


%% remove incorrect ones
ray3mismatch = ray_mismatch(calnofix, cam2d,  ncams);
cam2d(isnan(ray3mismatch),:,:) = [];

ray3mismatch = ray_mismatch(calnofix, cam2d_check,  ncams);
cam2d_check(isnan(ray3mismatch),:,:) = [];

%% key step: nonlinear search for the best parameters
calout = fminsearch(@(x) fitfunc(x,calfix,cam2d,ncams,fixcamera), calnofix, fmin_options);
calarray = calall;

if isempty(fixcamera) == 0
    left = setdiff(1:ncams,fixcamera);
    calarray(:,left) = calout;
else
    calarray = calout;
end

%% Check the quality of the final calibration
[dist3, dist1, pall] = ray_mismatch(calarray,cam2d,ncams);

% figure(5);
% hist(dist3,50);
% title('mismatch distribution after optimization on good matches');
% xlabel('mismatch (mm)')
allcams = mean(dist3)
cammismatch2 = zeros(ncams);
for n1 = 1 : ncams
cammismatch2(n1) = mean(dist1(n1,:));
end

[dist3, dist1, pall] = ray_mismatch(calarray, cam2d_check, ncams);

% figure(6);
% hist(dist3, 50);
% title('final mismatch distribution of data not used for the optimization');
% xlabel('mismatch (mm)')
allcams_check_with_mismatches = mean(dist3)


% Refill the camParaCalib structure and write it to a .cfg file.
%  (There turns out to be a very small change in the first camera even though it was not
%    dynamically calibrated since the rotation matrix is projected onto a matrix with
%     determinant exactly -1, so it changes just a bit.  The new R needs to be kept since the
%     dynamic calibration was done for this rotation matrix)
for icam = 1:ncams
    camParaCalib(icam).R = eul2rotm(calarray(1:3,icam));
    camParaCalib(icam).Rinv = inv(camParaCalib(icam).R);
    camParaCalib(icam).T = calarray(4:6,icam);
    camParaCalib(icam).Tinv = camParaCalib(icam).Rinv * (-1 * camParaCalib(icam).T);
    camParaCalib(icam).f_eff = calarray(7,icam);
    camParaCalib(icam).k1 = calarray(8,icam);
    camParaCalib(icam).err_x = cammismatch2(icam);
    camParaCalib(icam).err_y = cammismatch2(icam);
end

%fname_cfg=strcat(dist3D_path,'PTVSetup_optimized_',dist3D_stem,'.cfg');
%gv_write_calib_cfg(camParaCalib, ncams, fname_cfg);
% %% save the result
% save ([direc 'VSC_calib.mat'] , 'camParaCalib');
end

function rotm = eul2rotm(eul, sequence)
    if (size(eul,1) ~= 3)
        error('eul2rotm: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end

    if ~exist('sequence', 'var')
        % use the default axis sequence ...
        sequence = 'ZYX';
    end
    rotm = zeros(3,3);

    s_1 = sin(eul(1,1)); % theta_z or theta_z1
    c_1 = cos(eul(1,1));
    s_2 = sin(eul(2,1)); % theta_y
    c_2 = cos(eul(2,1));
    s_3 = sin(eul(3,1)); % theta_x or theta_z2
    c_3 = cos(eul(3,1));

    %% Convert the given Euler angles theta for the x, y and z-axis into the corresponding
    %  direction cosine rotation matrix R, in dependency of the axis rotation sequence for
    %  the multiplication order of the rotation factors:
    % For further details see:
    %   [1] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/EulerAngles.pdf>, p. 9 & 16.
    %   [2] MATLAB TOOLBOX FOR RIGID BODY KINEMATICS, Hanspeter Schaub & John L. Junkins,
    %       9th AAS/AIAA Astrodynamics Specialist Conference, AAS 99-139, 1999, <http://hanspeterschaub.info/Papers/mtb1.pdf>, p. 4.
    %   [3] GitHub: ShoolCode/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox,
    %       <https://github.com/ZachDischner/SchoolCode/tree/master/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox/>
    %   [4] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
    %       pp. 31-32, formulas (2.18) and (2.20).
    switch sequence
        case 'ZYX'
            %            |c_1*c_2    c_1*s_2*s_3 - s_1*c_3    c_1*s_2*c_3 + s_1*s_3|
            % R(Theta) = |s_1*c_2    s_1*s_2*s_3 + c_1*c_3    s_1*s_2*c_3 - c_1*s_3|
            %            |   -s_2                  c_2*s_3                  c_2*c_3|
            rotm(1,1) =  c_1*c_2;
            rotm(1,2) =  c_1*s_2*s_3 - s_1*c_3;
            rotm(1,3) =  c_1*s_2*c_3 + s_1*s_3;

            rotm(2,1) =  s_1*c_2;
            rotm(2,2) =  s_1*s_2*s_3 + c_1*c_3;
            rotm(2,3) =  s_1*s_2*c_3 - c_1*s_3;

            rotm(3,1) = -s_2;
            rotm(3,2) =  c_2*s_3;
            rotm(3,3) =  c_2*c_3;
        case 'ZYZ'
            %            |c_1*c_2*c_3 - s_1*s_3   -c_1*c_2*s_3 - s_1*c_3    c_1*s_2|
            % R(Theta) = |s_1*c_2*c_3 + c_1*s_3   -s_1*c_2*s_3 + c_1*c_3    s_1*s_2|
            %            |             -s_2*c_3                  s_2*s_3        c_2|
            rotm(1,1) =  c_1*c_2*c_3 - s_1*s_3;
            rotm(1,2) = -c_1*c_2*s_3 - s_1*c_3;
            rotm(1,3) =  c_1*s_2;

            rotm(2,1) =  s_1*c_2*c_3 + c_1*s_3;
            rotm(2,2) = -s_1*c_2*s_3 + c_1*c_3;
            rotm(2,3) =  s_1*s_2;

            rotm(3,1) = -s_2*c_3;
            rotm(3,2) =  s_2*s_3;
            rotm(3,3) =  c_2;
        otherwise
            error('eul2rotm: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end

function deviation = fitfunc(calpara, calfix, cam2d, ncams, fixcamera)

% This function is for use with dynamic calibration.  It combines the calibration
% parameters into a single array and then calls the routine to calculate the distances
%  between the rays corresponding to given 2d image plane coordinates.  
%
% inputs:

%   calarray     --  full 8 by 3 calibration matrix.  Columns are for each camera.  1:3 are angles, 4:6 are T, 7 is effective focal length (f_eff) and 8 is distortion (k1) 
%   cam2d        --  N by 2 by 3 array containing 2d image plane coordinates of
%                       N matched particles.   Last index is camera id.  Middle index is coordinate (h or w).
%
% outputs:
%   deviation    -- average of the 3D ray mismatch. 
%
%  called functions:
%    gv_calc_ray_mismatch
%change Sep6,2011, calarray(1:6,2:3)=calopt(1:6,1:2) for 3cameras
%to calarray(1:6,2:ncams)=calopt(1:6,:); for 4 cameras
%calarray(1:6,2:ncams)=calopt(1:6,:);   %here we assume camera 1 is fixed.  If this changes, it needs to be changed here also.
%calarray(1:6,3:4)=calopt(1:6,1:2);   %here we assume camera 1,2 are fixed.  If this changes, it needs to be changed here also.


%calarray = [calpara(1:5,:);fixpara;calpara(6,:)];

if isempty(fixcamera)==0
    calarray(:,fixcamera) = calfix;
    left = setdiff(1:ncams,fixcamera);
    calarray(:,left) = calpara;

else
    calarray = calpara;
end



ray3mismatch=ray_mismatch(calarray, cam2d,  ncams);

deviation = mean(ray3mismatch);
end

function [ray3mismatch,h,pall] = ray_mismatch(calarray, cam2d, ncams)

for icam = 1:ncams
    R(:,:,icam) = eul2rotm(calarray(1:3,icam));
    Rinv(:,:,icam) = inv(R(:,:,icam));
    Tinv(:,icam) = Rinv(:,:,icam) * (-1* calarray(4:6,icam));
end
npoints=size(cam2d,1);
ray3mismatch=zeros(npoints,1); %average deviation over all cameras
h = zeros(npoints, ncams); %deviation from each camera
point3D= zeros(npoints,3);

for np=1:npoints 
    M = zeros(3,3);
    pM = zeros(3,ncams);
    u = zeros(3, ncams);
for icam = 1:ncams
        % then find the unit vector designated by the camera position and
        % the 2D pixel coordinates. 
    %u(:,icam) = gv_imgplane2unitvector(calarray(:,icam), Rinv(:,:,icam), cam2d(np,:,icam));
    u(:,icam) = img2unitv(calarray(:,icam), Rinv(:,:,icam), cam2d(np,:,icam));
    uM = eye(3) - u(:,icam) * (u(:,icam))';
    pM(:,icam) = uM * Tinv(:,icam);
    M = M + uM;
end

if (det(M) < 10)
  det(M);  
end
    % find the point minimizing the distance from all rays
    p = M \ sum(pM,2);  % sums pm x together for all three cameras.  Makes a column vector, then does inv(M)*sum(pM,2)
    %find the distances from each ray.

    pall(:,np)=p;
for icam = 1:ncams
    temp = p - ((p') * u(:,icam)) * u(:,icam) - pM(:,icam);
    h(np,icam) = sqrt(temp' *temp);
end
    ray3mismatch(np) = sqrt(mean(h(np,:).*h(np,:)));
end
end

function u = img2unitv(calarray, Rinv, dat2din)

% function to find unit vectors pointing in the direction corresponding to
% 2d image plane coordinates. 
%Similar to gv_pixel2unitvector but accepts a calibration matrix rather
%than a structure, and accepts undistorted image plane coordinates rather than 
% pixels.  These differences are needed to use dyanmic calibration.
%
% inputs:
%    calarray  -- 8 by ncams array with calibration parameters
%   Rinv  -- Inverse rotation matrix--could be calculated from calarray,
%       but this is faster
%   dat2din   -- array(N by 2) containing undistorted image plane
%       coordinates

% turns a 2d position into a unit vector.  The position of a particle in 3D space
% is along the ray given by the point cal.Tinv + lambda*(unit vector) where
% lambda is any real number.

%create arrays 
L=size(dat2din,1);
pos = zeros(L,3);
u=zeros(L,3);

pos(:,1:2)=dat2din;
%include distortion.  k1 seems to be defined using the radial distance
%after distortion
%rather than using the undistorted radius (because calibProj_Tsai iterates) so we
%don't iterate going this way.  I still have some questions about this
%since http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/DIAS1/
%puts the (1+k1r_d^2) term as multiplying the distorted position rather than dividing as calibProjTsai
%does.  It is all conventions in the end though--just need to be
%consistent.

x = pos(:,1);
y = pos(:,2);
a = x./y;

y = (1-sqrt(1-4.*y.^2.*calarray(8).*(a.^2+1)))./(2.*y.*calarray(8).*(a.^2+1));
x = a.*y;


% radius2 = pos(:,1)^2 + pos(:,2)^2;
% pos(:,1:2) = pos(:,1:2)/(1+calarray(8)*radius2);

%Now scale by the effective focal length and the z coordinates

pos(:,1:2)=[x,y] .* (calarray(6) / calarray(7));
pos(:,3)=calarray(6);

%we also need to map the pinhole which is a common point on the ray for any
%pixel.  This is the same calculation as for pixin=[0 0], but the initial z
%coordinate is zero rather than cal.T(3)
 
zpoint = [0, 0, 0]; 
zpoint = (zpoint - (calarray(4:6))') * (Rinv)';

for i=1:L
    pos(i,:)=(pos(i,:) - (calarray(4:6))') * (Rinv)'; % I found this by reversing the steps in calibProj_Tsai
    u(i,:)=(pos(i,:)-zpoint)/norm(pos(i,:)-zpoint); %here we subtract zpoint and normalize to a unit vector
end
end

function eul = rotm2eul(rotm, sequence)
    if ( (size(rotm,1) ~= 3) || (size(rotm,2) ~= 3) )
        error('rotm2eul: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end

    if ~exist('sequence', 'var')
        % use the default axis sequence ...
        sequence = 'ZYX';
    end
    eul = zeros(3,1);

    %% Compute the Euler angles theta for the x, y and z-axis from a rotation matrix R, in
    %  dependency of the specified axis rotation sequence for the rotation factorization:
    % For further details see:
    %   [1] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/EulerAngles.pdf>, pp. 9-10, 16-17.
    %   [2] Computing Euler angles from a rotation matrix, Gregory G. Slabaugh, <http://www.staff.city.ac.uk/~sbbh653/publications/euler.pdf>
    %   [3] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
    %       pp. 30-33, formulas (2.19), (2.19'), (2.21) and (2.21').
    switch sequence
        case 'ZYX'
            % convention used by (*) and (**).
            % note: the final orientation is the same as in XYZ order about fixed axes ...
            if (rotm(3,1) < 1)
                if (rotm(3,1) > -1) % case 1: if r31 ~= �1
                    % Solution with positive sign. It limits the range of the values
                    % of theta_y to (-pi/2, pi/2):
                    eul(1,1) = atan2(rotm(2,1), rotm(1,1)); % theta_z
                    eul(2,1) = asin(-rotm(3,1));            % theta_y
                    eul(3,1) = atan2(rotm(3,2), rotm(3,3)); % theta_x
                else % case 2: if r31 = -1
                    % theta_x and theta_z are linked --> Gimbal lock:
                    % There are infinity number of solutions for theta_x - theta_z = atan2(-r23, r22).
                    % To find a solution, set theta_x = 0 by convention.
                    eul(1,1) = -atan2(-rotm(2,3), rotm(2,2));
                    eul(2,1) = pi/2;
                    eul(3,1) = 0;
                end
            else % case 3: if r31 = 1
                % Gimbal lock: There is not a unique solution for
                %   theta_x + theta_z = atan2(-r23, r22), by convention, set theta_x = 0.
                eul(1,1) = atan2(-rotm(2,3), rotm(2,2));
                eul(2,1) = -pi/2;
                eul(3,1) = 0;
            end
        case 'ZYZ'
            % convention used by (*)
            if (rotm(3,3) < 1)
                if (rotm(3,3) > -1)
                    % Solution with positive sign, i.e. theta_y is in the range (0, pi):
                    eul(1,1) = atan2(rotm(2,3),  rotm(1,3)); % theta_z1
                    eul(2,1) = acos(rotm(3,3));              % theta_y (is equivalent to atan2(sqrt(r13^2 + r23^2), r33) )
                    eul(3,1) = atan2(rotm(3,2), -rotm(3,1)); % theta_z2
                else % if r33 = -1:
                    % Gimbal lock: infinity number of solutions for
                    %   theta_z2 - theta_z1 = atan2(r21, r22), --> set theta_z2 = 0.
                    eul(1,1) = -atan2(rotm(2,1), rotm(2,2)); % theta_z1
                    eul(2,1) = pi;                           % theta_y
                    eul(3,1) = 0;                            % theta_z2
                end
            else % if r33 = 1:
                % Gimbal lock: infinity number of solutions for
                %    theta_z2 + theta_z1 = atan2(r21, r22), --> set theta_z2 = 0.
                eul(1,1) = atan2(rotm(2,1), rotm(2,2)); % theta_z1
                eul(2,1) = 0;                           % theta_y
                eul(3,1) = 0;                           % theta_z2
            end
        % case 'ZYZ-'
        %     % convention used by (**)
        %     if (rotm(3,3) < 1)
        %         if (rotm(3,3) > -1)
        %             % Variant with negative sign. This is a derived solution
        %             % which produces the same effects as the solution above.
        %             % It limits the values of theta_y in the range of (-pi,0):
        %             eul(1,1) = atan2(-rotm(2,3), -rotm(1,3)); % theta_z1
        %             eul(2,1) = -acos(rotm(3,3));              % theta_y (is equivalent to atan2(-sqrt(r13^2 + r23^2), r33) )
        %             eul(3,1) = atan2(-rotm(3,2),  rotm(3,1)); % theta_z2
        %         else % if r33 = -1:
        %             % Gimbal lock: infinity number of solutions for
        %             %   theta_z2 - theta_z1 = atan2(-r12, -r11), --> set theta_z2 = 0.
        %             eul(1,1) = -atan2(-rotm(1,2), -rotm(1,1)); % theta_z1  (correct ???)
        %             eul(2,1) = -pi;                            % theta_y
        %             eul(3,1) = 0;                              % theta_z2
        %         end
        %     else % if r33 = 1:
        %         % Gimbal lock: infinity number of solutions for
        %         %    theta_z2 + theta_z1 = atan2(-r12, -r11), --> set theta_z2 = 0.
        %         eul(1,1) = atan2(-rotm(1,2), -rotm(1,1)); % theta_z1  (correct ???)
        %         eul(2,1) = 0;                             % theta_y
        %         eul(3,1) = 0;                             % theta_z2
        %     end
        otherwise
            error('rotm2eul: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end

function x2D = UnDistort(X2D,camParaCalib) 
    
    kr = camParaCalib.k1;
    X = (X2D(:,1) - 512)*0.02;
    Y = (-X2D(:,2) + 512)*0.02;

    if (kr ~= 0)
        a = X ./ Y;

        Y = (1 - sqrt(1 - 4 .* Y.^2 .* kr.*(a.^2 + 1))) / (2.* Y.* kr.* (a.^2 + 1));
        X = a.*Y;
    end
    
    x2D = [X Y];
end

function tracks = PickGoodTracks(tracks, view_size)
    num_tracks = max(tracks(:, 5));
    std_error = view_size / 1000;
    num_good = ceil(num_tracks * .1);
    if num_tracks * .1 > 1000,  nun_good = 1000; end
    error_pool = zeros(1, num_good);
    no_good = 1;
    for i = 1 :  num_tracks
        track = tracks(tracks(:, 5) == i, 1 : 3);      
        if size(track, 1) < 20 
            tracks(tracks(:, 1) == i, :) = []; % delete short tracks
            continue; 
        end
        error0 = zeros(1, 3);
        num_fit = 10;
        for j = 1 : 3
            p = polyfit(1:num_fit, track(end - num_fit + 1:end, j)', 3);
            track_est = polyval(p, 1:num_fit);
            error0(j) = mean(abs(track_est - track(end - num_fit + 1:end, j)'));
        end
        error = norm(error0);
        
        if error > std_error
            tracks(tracks(:, 5) == i, :) = []; % delete bad tracks
        end
        
        if no_good <= num_good
            error_pool(1, no_good) = error;
            no_good = no_good + 1;
        elseif no_good == num_good + 1
            std_error = mean(error_pool) + 3 *  std(error_pool);
            no_good = num_good + 1;
        end
            
        if ~(mod(i, 1000))
            i
        end
    end
end

function X3D = Img2World(camParaCalib,X2D)

    s = [size(X2D,1);size(X2D,2)];
%     X2D = [X2D zeros(s(1),1)];
    tmp = X2D.*(camParaCalib.T(3)/camParaCalib.f_eff);
    proj = [tmp(:,1) tmp(:,2) repmat(camParaCalib.T(3),[s(1) 1])]';
    X3D = camParaCalib.Rinv*(proj - camParaCalib.T);

end

function CalibMatToTXT(camParaCalib, save_path)

fileID = fopen(save_path,'w');
fprintf(fileID, '# Camera configuration file\n');
fprintf(fileID, '# generated %s\n \n', datetime);

fprintf(fileID, '4    # camera number\n');

for i = 1 : 4
    fprintf(fileID, '\n#camera %d\n', i );
    fprintf(fileID,'%d    #Noffh\n',camParaCalib(i).Noffh);
    fprintf(fileID,'%d    #Noffw\n',camParaCalib(i).Noffw);
    fprintf(fileID,'%d    #Npixw\n',camParaCalib(i).Npixw);
    fprintf(fileID,'%d    #Npixh\n',camParaCalib(i).Npixh);
    fprintf(fileID,'%f    #wpix\n',camParaCalib(i).wpix);
    fprintf(fileID,'%f    #hpix\n',camParaCalib(i).hpix);
    fprintf(fileID,'%f    #f_eff\n',camParaCalib(i).f_eff);
    fprintf(fileID,'%f    #kr\n',camParaCalib(i).k1);
    fprintf(fileID,'%d    #kx\n',1);
    fprintf(fileID,'%f    #R\n',camParaCalib(i).R');
    fprintf(fileID,'%f    #T\n',camParaCalib(i).T);
    fprintf(fileID,'%f    #Rinv\n',camParaCalib(i).Rinv');
    fprintf(fileID,'%f    #Tinv\n',camParaCalib(i).Tinv);
end

fclose(fileID);
end
