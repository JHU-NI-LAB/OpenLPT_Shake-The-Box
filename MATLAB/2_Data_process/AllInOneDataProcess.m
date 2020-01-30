function AllInOneDataProcess(data_path, parentfolder, initial_Calib_filepath, job_num)
% data_path, remember to put '/' after
%initial_Calib_filepath format: ~/VSC_Calib_XXXXXX, don't include the file
%extension like: .txt or .mat
%parentfolder is the folder for the project
%% Upload data from file server to MARCC
% create the work path on MARCC
workpath = extractAfter(data_path, parentfolder);
workpathfolder = strsplit(workpath, '/');

bubble_label = 0;
if contains(data_path, 'Bubbles')
    bubble_label = 1;
end

% check whether the folder exist
[~, result] = system('cd ~/scratch && pwd');
marcc_orgpath = [result(1 : end - 1) '/'];
if ~exist([marcc_orgpath workpathfolder{1}], 'dir') 
    mkdir([marcc_orgpath workpathfolder{1}]);
end

if bubble_label
    if ~exist([marcc_orgpath workpathfolder{1} '/' workpathfolder{2}], 'dir') 
        mkdir([marcc_orgpath workpathfolder{1} '/' workpathfolder{2}]);
    end
end

marcc_workpath = [marcc_orgpath, workpath];
if ~exist(marcc_workpath, 'dir')
    mkdir(marcc_workpath);
end

log_filepath = [marcc_workpath 'Log.txt'];
if exist(log_filepath, 'file') 
    text = fileread([marcc_workpath 'Log.txt']);
else
    text = 'The data process has been started!';
    file_ID = fopen(log_filepath, 'w');
    fprintf(file_ID, 'The data process has been started!\n');
    %fclose(file_ID);
end

if ~contains(text, 'TestData has been uploaded to MARCC!')
    data_task_id = UploadDataToMARCC(data_path, marcc_workpath);
    if bubble_label 
        mask_path = replace(data_path, 'Data', 'Processed_Images');
        mask_task_id = UploadDataToMARCC(mask_path, marcc_workpath);
    end
end

%% Check whether VSC has been done. If not, do the VSC
VSC_file_path = [marcc_orgpath workpathfolder{1} '/ForVSC/VSC_Calib_' workpathfolder{1} '.txt'];
VSC_folder_path = [marcc_orgpath workpathfolder{1} '/ForVSC/'];
if ~exist(VSC_file_path, 'file') 
    %Upload data for VSC on MARCC
    if ~exist(VSC_folder_path, 'dir') 
        mkdir(VSC_folder_path);
    end
    VSC_rawimage_path = [extractBefore(data_path, workpathfolder{1}), workpathfolder{1}, '/ForVSC/'];
    if ~contains(fileread(log_filepath), 'VSCData has been uploaded to MARCC and processed!')
        VSC_data_task_id = UploadDataToMARCC(VSC_rawimage_path, VSC_folder_path);
        WaitDataTransfer(VSC_data_task_id);
        %Preprocess the image
        PreprocessImage(VSC_folder_path, 'Cam', '.tiff', VSC_folder_path, 1000);
        file_ID = fopen(log_filepath, 'a');
        fprintf(file_ID, 'VSCData has been uploaded to MARCC and processed!\n');
        %fclose(file_ID);
    end
    copyfile([initial_Calib_filepath '.txt'], VSC_folder_path);
    %Generate Configure files
    calibration_name = ['VSC_Calib' extractAfter(initial_Calib_filepath, 'VSC_Calib') '.txt'];
    GenerateConfigFile(VSC_folder_path, 1, 1000, calibration_name);
    %Run STB
    if ~exist([VSC_folder_path 'result1.txt'], 'file') || ~CheckSTB([VSC_folder_path 'result1.txt'])
        RunSTB(VSC_folder_path);
%         GenerateJobs(VSC_folder_path, ['PT' workpathfolder{1} 'VSC'], 1, 1);
%         [~, result] = system(['cd ' VSC_folder_path ' && ' 'sbatch ' 'PT' workpathfolder{1} 'VSC_1.sh']);
%         if ~contains(result, 'Submitted batch job')
%             disp('Job for VSC submision failed!');
%         else
%             disp('Job for VSC submision succeeded!');
%             job_id = extractAfter(result, 'job ');
%             pause(120);
%             [~, result] = system(['qstat -a ' job_id]);
%             if contains(result, '00:00') 
%                 % it means the job is still in queue then try to run the code in
%                 % the local computer
%                 RunSTB(VSC_folder_path);
%                 system(['scancel ' job_id]); %cancel the job
%             end
%         end
         pause(60);
         % Wait until it is finished
         while ~exist([VSC_folder_path 'result1.txt'], 'file') || ~CheckSTB([VSC_folder_path 'result1.txt'])
             pause(60); %Check every minute
         end
    end
    % Do VSC
    copyfile([initial_Calib_filepath '.mat'], VSC_folder_path);
    calibration_file = [VSC_folder_path 'VSC_Calib' extractAfter(initial_Calib_filepath, 'VSC_Calib'), '.mat'];
    RunVSC(VSC_folder_path, 1000, calibration_file, workpathfolder{1});
end


%% Copy the VSC calibration file to local folder
% copy the VSC file to the work path
copyfile([VSC_folder_path 'VSC_Calib_' workpathfolder{1} '.txt'], marcc_workpath);
copyfile([VSC_folder_path 'VSC_Calib_' workpathfolder{1} '.mat'], marcc_workpath);
copyfile([VSC_folder_path 'VSC_Calib_' workpathfolder{1} '.txt'], '~/scratch/VSC/');
copyfile([VSC_folder_path 'VSC_Calib_' workpathfolder{1} '.mat'], '~/scratch/VSC/');

%% Preprocess the raw image
% Upload data to work directory on MARCC
% UploadDataToMARCC(data_path, marcc_workpath);
if ~contains(text, 'TestData has been uploaded to MARCC!')
     WaitDataTransfer(data_task_id);
     if bubble_label
         WaitDataTransfer(mask_task_id);
     end
     file_ID = fopen(log_filepath, 'a');
     fprintf(file_ID, 'TestData has been uploaded to MARCC!\n');
     %fclose(file_ID);
end

%% Generate the configure files
%Preprocess the image
if ~contains(fileread(log_filepath), 'Image preprocessing is complete!')
    if ~bubble_label
        PreprocessImage(marcc_workpath, 'Cam', '.tiff');
    else
        PreprocessBubbleImg(marcc_workpath, marcc_workpath);
    end
     file_ID = fopen(log_filepath,'a');
     fprintf(file_ID, 'Image preprocessing is complete!\n');
     %fclose(file_ID);
end

%% Run STB for 100 frames and do VSC
% Generate Configure files
if ~contains(fileread(log_filepath), 'STB for data VSC is complete!')
    GenerateConfigFile(marcc_workpath, 1, 100, ['VSC_Calib_' workpathfolder{1} '.txt']);
    % RunSTB
    if ~exist([marcc_workpath 'result1.txt'], 'file') || ~CheckSTB([marcc_workpath 'result1.txt'])
    %     if ~contains(fileread(log_filepath), 'STB for data VSC is running!')
            RunSTBForVSC(marcc_workpath);
    %         fileID = fopen(log_filepath,'a');
    %         fprintf(fileID, 'STB for data VSC is running!\n');
    %         fclose(file_ID);
    %     end
        % GenerateJobs(marcc_workpath, ['PT' workpathfolder{1}], 1, 2);
        % [~, result] = system(['cd ' marcc_workpath ' && ' 'sbatch ' 'PT' workpathfolder{1} '_1.sh']);
        % 
        % if ~contains(result, 'Submitted batch job')
        %     disp('Job for VSC submision failed!');
        % else
        %     disp('Job for VSC submision succeeded!');
        %     job_id = extractAfter(result, 'job ');
        %     pause(120);
        %     [~, result] = system(['qstat -a ' job_id]);
        %     if contains(result, '00:00') 
        %         % it means the job is still in queue then try to run the code in
        %         % the local computer
        %         RunSTB(marcc_workpath);
        %         system(['scancel ' job_id]); %cancel the job
        %     end
        % end
        pause(60);
        % Wait until it is finished
        while ~exist([marcc_workpath 'result1.txt'], 'file') || ~CheckSTB([marcc_workpath 'result1.txt'])
             pause(60); %Check every minute
        end
        file_ID = fopen(log_filepath,'a');
        fprintf(file_ID, 'STB for data VSC is complete!\n');
        %fclose(file_ID);
    end
end
% Do VSC
if ~contains(fileread(log_filepath), 'VSC for data is complete!')
    RunVSC(marcc_workpath, 100, [marcc_workpath 'VSC_Calib_' workpathfolder{1} '.mat'], [workpathfolder{1} '_' workpathfolder{2}]);
    file_ID = fopen(log_filepath,'a');
    fprintf(file_ID, 'VSC for data is complete!\n');
    system(['rm -rf ' marcc_workpath 'Tracks/ConvergedTracks/*.txt']); %delete all the txt file for VSC
    %fclose(file_ID);
end

%% Creat the jobs
% regenerate config file
ncams = 4;
for i = 1 : ncams
    a = dir([marcc_workpath 'cam' num2str(i) '/*.tif']);
    total_imgs(i) = numel(a);
end
GenerateConfigFile(marcc_workpath, 1, min(total_imgs), ['VSC_Calib_' workpathfolder{1} '_' workpathfolder{2} '.txt']);
GenerateJobs(marcc_workpath, ['PT' workpathfolder{1} '_' workpathfolder{2}], job_num);
% Submit the jobs
for i = 1 : job_num
    [~, result] = system(['cd ' marcc_workpath ' && ' 'sbatch ' 'PT' workpathfolder{1} '_' workpathfolder{2} '_' num2str(i) '.sh']);
    if ~contains(result, 'Submitted batch job')
        disp(['Job' num2str(i) ' submision failed!']);
    else
        disp(['Job' num2str(i) ' submision succeeded!']);
    end
end

% %% When jobs are done, download the result to file server
% % Wait until it is finished
% job_done = 0;
% while ~job_done
%      label = 1;
%      for i = 1 : 7
%          if ~exist([marcc_workpath 'result' num2str(i) '.txt'], 'file')
%               label = 0;
%              continue;
%          end
%          label = CheckSTB([marcc_workpath 'result' num2str(i) '.txt']) & label;
%      end
%      if label == 1
%          job_done = 1;
%      end
%      pause(60); %Check every minute
% end
% 
% % Download to file sever
% save_path = ['/E/share2/Projects/1-Bubble/CamConfig_of_07.13.18/Processed_tracks/' extractAfter(data_path, 'Data/')];
% DownloadDataToFileServer([marcc_workpath '/Tracks'], save_path);
% disp('Finished!')

end


function finish = WaitDataTransfer(task_id)
finish = 0;
while ~finish
    [~, result] = system(['module load globus-cli && ' ...
    'globus login &&' ...
    'globus task show ' task_id]);
    if contains(result, 'SUCCEEDED')
        finish = 1;
        disp('Transfer finished!');
    else 
        pause(60);
    end
end
end
