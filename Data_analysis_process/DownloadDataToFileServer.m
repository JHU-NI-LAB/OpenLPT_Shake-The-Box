function DownloadDataToFileServer(path_marcc,path_file_server)
% format of path_file_server:
% "/E/share2/Projects/1-Bubble/CamConfig_of_07.13.18/Data/07.17.18/Run1_11_of_88_3s"
% format of path_marcc: "~/scratch/07.17.18/"
ep1 = 'e92d0fde-2f33-11e7-bc9e-22000b9a448b'; % endpoint for marcc
ep2 = '70ddb4b8-bd33-11e8-8c1c-0a1d4c5c824a'; % endpoint for file sever

[~,result] = system(['module load globus-cli && ' ...
    'globus login &&' ...
    'globus endpoint activate ' ep1 '&&' ...
    'globus transfer ' ep1 ':' path_marcc ' ' ep2 ':' path_file_server ' --recursive --sync-level size']);
if ~contains(result,'You are already logged in!')
    disp('Login failed!');
else
    disp('Login succeeded!');
end

if ~contains(result, 'Endpoint is already activated.')
    disp('Please activate the endpoint of MARCC first!');
else
    disp('MARCC endpoint is currently activated.');
end

if ~contains(result, 'The transfer has been accepted')
    disp('Transfer request failed!');
    disp(result);
else
    disp('Transfer request accepted!');
    disp('Data is transferring, please wait...');
    %get the task ID
    task_id = extractAfter(result, 'ID: ');
end
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

