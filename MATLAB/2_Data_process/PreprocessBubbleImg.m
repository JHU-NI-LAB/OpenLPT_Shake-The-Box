function PreprocessBubbleImg(dirr, mask_dir, savedir, totalImgs)
%PREPROCESSBUBBLEIMG Summary of this function goes here
if ~exist('savedir', 'var')
    savedir = dirr;
end
% determine which section it is
section_no = extractBetween(dirr, '/Bubbles', '/');
if ~exist('totalImgs', 'var')
    a = dir([dirr 'C001_H001S000' section_no{1} '/*.tif']);
    totalImgs = numel(a) ; 
end

camsave_dir = {[savedir 'cam1/'],[savedir 'cam2/'],[savedir 'cam3/'],[savedir 'cam4/'], [savedir 'cam5/'],[savedir 'cam6/'] }';
for i = 1 : 4
    if ~exist(camsave_dir{i}, 'dir')
        mkdir(camsave_dir{i});
    end
end
print_dir = {['cam1/'],['cam2/'],['cam3/'],['cam4/'],['cam5/'],['cam6/'] }';
ncams = 4;

log_filepath = [dirr 'Log.txt'];
start = 1;
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
    
%% get the background
totalImgs_for_background = 500;
ncams = 4;
Npixh = 1024;
Npixw = 1024;
bak = zeros(Npixh,Npixw,ncams);
%% Background images (using normal images)
% in the form of bak(:,:,cam)
for cam = start : ncams
    missed_pixels = zeros(Npixh,Npixw);
    num = 0;
    camdir = [dirr 'C00' num2str(cam) '_H001S000' section_no{1} '/'];
    for I = 1:totalImgs_for_background
        if (mod(I,10) == 0)
            load([mask_dir 'Cam' num2str(cam) '.mat'],['Image' num2str(I,'%06.0f')]);
            
             % getting the mask
            eval(['mask = Image' num2str(I,'%06.0f') ';']);

            % modifying the mask to remove extra residual around the bubble
            mask1 = mask;
            [row, col] = find(mask == 1);
            rows = zeros(size(row,1),20); cols = zeros(size(col,1),20); % extending the radius of bubbles by 10 px
            for i = 1:10
                rows(:,i+10) = min(row + i,Npixh);
                rows(:,i) = max(1,row - i);
            end
            for i = 1:10
                cols(:,i+10) = min(Npixw,col + i);
                cols(:,i) = max(1,col - i);
            end
            linearInd = sub2ind(size(mask1),rows,cols);
            mask1(linearInd) = 1;
        
            imgdir = [camdir 'C00' num2str(cam) '_H001S000' section_no{1} '0' num2str(I,'%05.0f') '.tif'];
            img = imread(imgdir);
            img = uint8(double(img).*(~mask1));
            missed_pixels = missed_pixels + double(mask1);
            bak(:,:,cam) = bak(:,:,cam) + double(img);
            num = num + 1;
        end
    end
    tot_num = ones(Npixh,Npixw).*num;
    tot_num = tot_num - missed_pixels;
    bak(:,:,cam) = bak(:,:,cam)./tot_num;
%     bak(:,:,cam) = bak(:,:,cam)/num;
    bak1 = uint8(bak(:,:,cam));
    bak(:,:,cam) = bak1;
end
bak = uint8(bak);
clear Image*; % clear the image variable to free the memory

% figure;
% for i = 1 : 4
%     subplot(2,2,i);
%     imshow(bak(:, :, i));
% end

save([savedir 'background.mat'],'bak');

%% Background images
% in the form of bak(:,:,cam)
% correction = [1 1 0 0 0 0];
% threshold = [0.35; 0.45; 0.42; 0.4; 0.18; 0.21];
%% processing
for cam = start : ncams
    h = load([mask_dir 'Cam' num2str(cam) '.mat']);
    totalImgs = numel(fieldnames(h));
    clear h;
    camdir = [dirr 'C00' num2str(cam) '_H001S000' section_no{1} '/'];
    parfor I = 1: 1 : totalImgs
        %load([mask_dir 'Cam' num2str(cam) '.mat'],['Image' num2str(I,'%06.0f')]);
        imgdir = [camdir 'C00' num2str(cam) '_H001S000' section_no{1} '0' num2str(I,'%05.0f') '.tif'];
%         imgdir = [camdir 'cam' num2str(cam) 'frame' num2str(I,'%04.0f') 'raw.tif'];
        img = imread(imgdir);
%         img1 = -double(bak(:,:,cam)) + double(img);
        img1 = double(bak(:,:,cam)) - double(img);
        
        h = load([mask_dir 'Cam' num2str(cam) '.mat'],['Image' num2str(I,'%06.0f')]);
        h = struct2cell(h);
        mask = h{1};
%         % getting the mask
%         if exist(['Image' num2str(I,'%06.0f')], 'var')
%             eval(['mask = Image' num2str(I,'%06.0f') ';']);
%         else
%             break; % skip those frame without masks
%         end

        % modifying the mask to remove extra residual around the bubble
%         mask1 = mask;
%         [row, col] = find(mask == 1);
%         rows = zeros(size(row,1),20); cols = zeros(size(col,1),20); % extending the radius of bubbles by 10 px
%         for i = 1:10
%             rows(:,i+10) = min(row + i,1024);
%             rows(:,i) = max(1,row - i);
%         end
%         for i = 1:10
%             cols(:,i+10) = min(1024,col + i);
%             cols(:,i) = max(1,col - i);
%         end
%         linearInd = sub2ind(size(mask1),rows,cols);
%         mask1(linearInd) = 1;
%         img2 = (double(img1).*(~mask1));
%         
%         img3=img2./max(max(img2)).*255;
%         img3 = uint8(img3);
% %         img4 = imadjust(img3,[0 threshold(cam)]);
%         img4 = LaVision_ImgProcessing(img3);
        mask1 = mask;
%         [row, col] = find(mask == 1);
%         rows = zeros(size(row,1),20); cols = zeros(size(col,1),20); % extending the radius of bubbles by 10 px
%         for i = 1:10
%             rows(:,i+10) = min(row + i,1024);
%             rows(:,i) = max(1,row - i);
%         end
%         for i = 1:10
%             cols(:,i+10) = min(1024,col + i);
%             cols(:,i) = max(1,col - i);
%         end
%         linearInd = sub2ind(size(mask1),rows,cols);
%         mask1(linearInd) = 1;
        img2 = (double(img1).*(~mask1));
        
        img3=img2./max(max(img2)).*255;
        mask2 = bwareaopen(img3>30,70);
        img3 = (img3.*(~mask2));
        img3 = uint8(img3);
%         img4 = imadjust(img3,[0 threshold(cam)]);
        img4 = LaVision_ImgProcessing(img3);
        
        if (cam == 1 || cam == 2)
            imwrite(img4,[camsave_dir{cam} 'cam' num2str(cam) 'frame' num2str(I - 1 ,'%05.0f') '.tif']);
        else
            imwrite(img4,[camsave_dir{cam} 'cam' num2str(cam) 'frame' num2str(I ,'%05.0f') '.tif']);            
        end
%             [h,f] = imhist(img4);
%         if (mod(I,50) == 0)
%             I
%         end
    end
       disp(['cam' num2str(cam) ' has been preprossed.'])
    % make image sequence txt files
    fID = fopen([savedir 'cam' num2str(cam) 'ImageNames.txt'],'w');
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

