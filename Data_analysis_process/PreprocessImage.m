function PreprocessImage(dirr, save_dir, totalImgs)
if ~exist('save_dir', 'var')
    save_dir = dirr;
end

if ~exist('totalImgs', 'var')
    a = dir([dirr 'C001_H001S0001/*.tif']);
    totalImgs = numel(a); 
end
%% Inputs
camsave_dir = {[save_dir 'cam1/'],[save_dir 'cam2/'],[save_dir 'cam3/'],[save_dir 'cam4/'], [save_dir 'cam5/'],[save_dir 'cam6/'] }';
ncams = 4;
for i = 1 : ncams
    if ~exist(camsave_dir{i}, 'dir')
        mkdir(camsave_dir{i});
    end
end
print_dir = {['cam1/'],['cam2/'],['cam3/'],['cam4/'],['cam5/'],['cam6/'] }';

log_filepath = [dirr 'Log.txt'];
start = 1;
if exist(log_filepath, 'file')
if contains(fileread(log_filepath), 'cam1 has been processed!')
    if contains(fileread(log_filepath), 'cam2 has been processed!')
        if contains(fileread(log_filepath), 'cam3 has been processed!')
            start = 4;
        else
            start = 3;
        end
    else
        start = 2;
    end
end
end
%% Getting the background
totalImgs_for_background = 1000;
% ncams = 4;
Npixh = 1024;
Npixw = 1024;
bak = zeros(Npixh,Npixw,ncams);
%% Background images (using normal images)
% in the form of bak(:,:,cam)

for cam = start : ncams
    num = 0;
    camdir = [dirr 'C00' num2str(cam) '_H001S0001/'];
    for I = 1:totalImgs_for_background
        if (mod(I,20) == 0)
            imgdir = [camdir 'C00' num2str(cam) '_H001S00010' num2str(I,'%05.0f') '.tif'];
            img = imread(imgdir);
            bak(:,:,cam) = bak(:,:,cam) + double(img);
            num = num + 1;
        end
    end
    bak(:,:,cam) = bak(:,:,cam)/num;
    bak1 = uint8(bak(:,:,cam));
    bak(:,:,cam) = bak1;
end
bak = uint8(bak);

% figure;
% for i = 1 : 4
%     subplot(2,2,i);
%     imshow(bak(:, :, i));
% end

save([save_dir 'background.mat'],'bak')

%% processing
correction = [1 1 0 0 0 0];
threshold = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5];
for cam = start:ncams
    camdir = [dirr 'C00' num2str(cam) '_H001S0001/'];
    % adjust threshold
    low = 0;
    high = 1;
    h_end = 3;
    while h_end > 2 || h_end == 0 
        imgdir = [camdir 'C00' num2str(cam) '_H001S0001000001.tif'];
        img = imread(imgdir);
        img1 = double(bak(:,:,cam)) - double(img);
        img3 = uint8(img1);
        img4 = imadjust(img3,[0 threshold(cam)]);
        img4 = LaVision_ImgProcessing(img4);
        [h,~] = imhist(img4);
        h_end = h(end);
         if h_end > 1
            low = threshold(cam);
            threshold(cam) = (threshold(cam) + high) / 2;
        elseif h_end == 0
            high = threshold(cam);
            threshold(cam) = (threshold(cam) + low) / 2;
        end
    end
    parfor I = 1:1:totalImgs
        imgdir = [camdir 'C00' num2str(cam) '_H001S00010' num2str(I,'%05.0f') '.tif'];
%         imgdir = [camdir 'cam' num2str(cam) 'frame' num2str(I,'%04.0f') 'raw.tif'];
        img = imread(imgdir);
%         img1 = -double(bak(:,:,cam)) + double(img);
        img1 = double(bak(:,:,cam)) - double(img);
%         img2=img1./max(max(img1)).*255;
        img3 = uint8(img1);
        img4 = imadjust(img3,[0 threshold(cam)]);
        img4 = LaVision_ImgProcessing(img4);
        if ~exist(camsave_dir{cam}, 'dir')
            mkdir(camsave_dir{cam});
        end
        imwrite(img4,[camsave_dir{cam} 'cam' num2str(cam) 'frame' num2str(I - correction(cam) ,'%05.0f') '.tif']);           
    end
    disp(['cam' num2str(cam) ' has been preprossed.'])
    % make image sequence txt files
    fID = fopen([save_dir 'cam' num2str(cam) 'ImageNames.txt'],'w');
    a = dir(camdir);
    frame = 1;
    for i = 1:size(a,1)
        if (a(i).isdir == 0)
%           fprintf(fID, ['C00' num2str(cam) 'H001S0001/' '%s\n'], a(i).name);
            fprintf(fID, [print_dir{cam} 'cam' num2str(cam) 'frame' num2str(frame,'%05.0f') '.tif\n']);
            frame = frame + 1;
        end
    end
    fclose(fID);
        file_ID = fopen(log_filepath, 'a');
    fprintf(file_ID, ['cam' num2str(cam) ' has been processed!\n']);
end
end

function outImg = LaVision_ImgProcessing(a)
%     for cam = 1:4
%         for I = 1:200
%             a = [cam_dir{cam} 'cam' num2str(cam) 'frame' num2str(I,'%04.0f') '.tif'];
%             a = double(imread(a));
            % subtrack sliding minimum
            a = double(a);
            b = imerode(a, true(3));
            c = a - b;
            b = imerode(a, true(3));
            c = c - b;

            % Gaussian smoothing filter
            d = imgaussfilt(c);

            % image normalization
            filt_size = 100;
            filt = 1/filt_size^2 .*true(filt_size);
            e = imfilter(d,filt);

            f = a - e;

            % sharpen the image
            outImg = uint8(imsharpen(f));
end
