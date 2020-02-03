function VSCOTF(data)
% clc,clear, close all

%% OTF User Inputs (for 800 frames)
% Volume and Subvolumes configuration
% The range of actual x position will be from -10 to 10 if you set x=20
% x=60;        % [mm] physical x dimension of a volume
% y=60;        % [mm] physical y dimension of a volume
% z=50;        % [mm] physical z dimension of a volume
x_range = [-80, 80];
y_range = [-80, 80];
z_range = [-50, 50];
n=3;         % define the number of subvolume (1D) in each xyz direction

% Gaussian fitting parameters
x_size=2;    % x range for 2D gaussian fit. so it will be +1 and -1 from the center which gives you 3x3 matrix for fitting
y_size=2;    % y range for 2D gaussian fit. so it will be +1 and -1 from the center which gives you 3x3 matrix for fitting
ncam=4;      % Number of camera
nframe = 200; % number of frames 
% name for image, for example: the 4 images from 4 cameras for the first...
... frame should be named as frame1temp1, frame1temp2, frame1temp3,frame1temp4
% Naming of the image files
% Modify the image name in OTFParam
% framename='frame';
% imagename ='cam'; 
% pos3Dname = 'pos3D';
% Naming for the output file
outputfiletitle='OTFParameter_frame';
image_filepath = '/home/tanshiyong/Documents/Data/Mass-Transfer/2018-11-12/Run1/';
calib_filepath = '/home/tanshiyong/Documents/Data/Mass-Transfer/2018-11-12/Run1/VSC_Calib_111218Run1.mat';
OTF_filepath = [image_filepath 'OTF_files/'];
averageOTFfilename = [OTF_filepath 'OTFParameters_VSC.txt'];
camParaCalib = load(calib_filepath);
camParaCalib = camParaCalib.camParaCalib;
if ~exist(OTF_filepath,'dir')
    mkdir(OTF_filepath);
end


%% VSC User Inputs 
%{ 
Variables that define the disparity map
dx=-epsln_x:stepsize:epsln_x;
dy=-epsln_y:stepsize:epsln_y;
exsize=numel(dx);
eysize=numel(dy);
%}
% epsln_x=3;
% epsln_y=3;
% stepsize=0.25;
% % Debug
% dx=-epsln_x:stepsize:epsln_x;
% dy=-epsln_y:stepsize:epsln_y;
% exsize=numel(dx);
% eysize=numel(dy);
% %

% Define the range of peak fitting 
% VSC_xrange = 2;
% VSC_yrange = 2;

%% Functions that evaluate OTF Parameters 
% Create .mat file for OTF Parameters of each frame
fileID = fopen([OTF_filepath 'data_bin.mat'], 'w');
fwrite(fileID, data, 'double'); % save the data
fclose(fileID); 
[row, col] = size(data);
% map a variable is to save memory and enable to get more workers for
% parallelization
 data_map = memmapfile([OTF_filepath 'data_bin.mat'], 'Format',{'double',[row col],'tracks'}); 
 data = [];
%  addpath /home/tanshiyong/Documents/Code/MATLAB/Post_analysis/SoundZone_Tools-master;
%  showTimeToCompletion; startTime=tic;
%  percent = parfor_progress(nframe);
tic
parfor i = 1:nframe 
%         percent = parfor_progress;
%         showTimeToCompletion( percent/100, [], [], startTime );
        int = num2str(i);
%         Modify this name below to match 
%         load([pos3Dname framename int  '.mat']); %.mat file generated from c++ that has actual 3D positions and 2D positions from n numbers of cameras
%         pos3D=cell2mat(eval([pos3Dname framename int])); %This will convert cell array to 2D matrix
        pos3D = data_map.Data.tracks(data_map.Data.tracks(:,4) == i, 1:3); % get the particles in frame i
        for j = 1 : ncam
            pos3D(:, (j - 1) * 2 + 4 : (j - 1) * 2 + 5) = calibProj(camParaCalib(j), pos3D(:, 1:3)); % get the 2D posistion
        end
        pos3D = DeleteOverlapParticle(pos3D, image_filepath, i);
        outputfilename = [OTF_filepath outputfiletitle int '.mat'];
        OTFParam(x_size,y_size,ncam,x_range,y_range,z_range,n,pos3D,image_filepath,outputfilename,i); 
%         [subvol_partdispar,num,temp,particlesubvolumedata] = OTFParam(xrange,yrange,ncam,x,y,z,n,pos3D,image_filepath,outputfilename,i); 
%         [I.(['frame' num2str(i)]),exsize,eysize] = VSC_disparitymap(subvol_partdispar,ncam,num,n,particlesubvolumedata,temp,epsln_x,epsln_y,stepsize); % disparity map 
end
% close(h)
toc

% Examine and create .txt file for average OTF Parameters and XYZ subvolume vectors       
OTFaverage(nframe,ncam,[OTF_filepath outputfiletitle],averageOTFfilename)

% % Save intensity value for each frame of different cameras
% save 'I_eachframe.mat' I
% % Sum disparity maps for each subvolume and all cameras over many frames 
% I_sumframe = VSCsum(nframe,ncam,I,exsize,eysize);
% 
% % Save the disparity map just in case want to run 2D gaussian fitting to
% % find dx and dy for the peak for each subvolume 
% save 'I_sumframe.mat' I_sumframe
% %% Generate image to view the intensity map for each subvolume 
% % Create matrices that combine disparity plot for all the subvolumes on each xy plane
% for i = 1:ncam
%     I_matrices = cell2mat(I_sumframe.(['cam' num2str(i)]));
%     for j = 1:n
%         figure('Name',['Disparity map for xy plane ' num2str(j) ' of cam ' num2str(i)],'NumberTitle','off');
%         imshow(I_matrices(:,:,j))
%     end
% end
% 
% % Determine the subpixel disparity peak location which gives dx and dy
% I_subpix = VSCsubdisparpeak(exsize,eysize,stepsize,ncam,nframe,VSC_xrange,VSC_yrange,I_sumframe);
end