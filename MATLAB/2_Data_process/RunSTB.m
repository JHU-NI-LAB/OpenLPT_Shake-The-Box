function RunSTB(work_dir)
% dir is the address of data
dir_code = '/home-4/stan26@jhu.edu/work/Code/LocalCode/ShakeTheBox/Release/';
% system(['cd ' dir_code]);

a = dir([work_dir '/Tracks/ConvergedTracks/ActiveLongTracks*.txt']);
if ~isempty(a)
    num_frame = extractBetween(a(end).name, 'Tracks', '.txt'); % to get how many frames have been processed
else 
    num_frame = [];
end
fileID = fopen([work_dir 'STBforVSC.sh'], 'w');
if isempty(num_frame)
    txt = ['#!/bin/bash -l \n' 'module load gcc/6.4.0 \n' ...
        'cd ' dir_code '\n' './ShakeTheBox ' work_dir 'trackingconfig1.txt > ' work_dir 'result1.txt << EOF \n'...
        '0\n' '0\n' 'EOF\n'];
else
    txt = ['#!/bin/bash -l \n' 'module load gcc/6.4.0 \n' ...
        'cd ' dir_code '\n' './ShakeTheBox ' work_dir 'trackingconfig1.txt > ' work_dir 'result1.txt << EOF \n'...
        '7\n' num_frame{1} '\n' '0\n' 'EOF\n'];
end
fprintf(fileID, txt);
fclose(fileID);
system(['cd ' work_dir ' && chmod +x ./STBforVSC.sh && ./STBforVSC.sh &']);


end

