function RunSTBForVSC(work_dir)
% dir is the address of data
dir_code = '/home-4/stan26@jhu.edu/work/Code/LocalCode/ShakeTheBox/Release/';
% system(['cd ' dir_code]);
a = dir([work_dir '/Tracks/ConvergedTracks/ActiveLongTracks*.txt']);
num_frame = size(a, 1); % to get how many frames have been processed
fileID = fopen([work_dir 'STBforVSC.sh'], 'w');
if num_frame == 0
    txt = ['#!/bin/bash -l \n' 'module load gcc/6.4.0 \n' ...
        'cd ' dir_code '\n' './ShakeTheBox ' work_dir 'trackingconfig1.txt > ' work_dir 'result1.txt << EOF \n'...
        '0\n' '0\n' 'EOF\n'];
else
    txt = ['#!/bin/bash -l \n' 'module load gcc/6.4.0 \n' ...
        'cd ' dir_code '\n' './ShakeTheBox ' work_dir 'trackingconfig1.txt > ' work_dir 'result1.txt << EOF \n'...
        '7\n' num2str(num_frame + 4) '\n' '0\n' 'EOF\n'];
end
fprintf(fileID, txt);
fclose(fileID);
system(['cd ' work_dir ' && chmod +x ./STBforVSC.sh && ./STBforVSC.sh &']);
end

