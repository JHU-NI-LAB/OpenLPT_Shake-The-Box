function RunVSC(dir, frame_num, calibration_file, save_name)
%dir is the folder for the VSC data
% save_name is the name for the new calibration file
% Read all the tracks
filepath_ALT = [dir 'Tracks/ConvergedTracks/ActiveLongTracks' num2str(frame_num) '.txt'];
active_long_tracks = ReadTracks(filepath_ALT);
filepath_IALT = [dir 'Tracks/ConvergedTracks/InactiveLongTracks'];
if (frame_num > 100)
	filepath_IALT = [dir 'Tracks/ConvergedTracks/InactiveLongTracks'];
	inactive_long_tracks = ReadSeparateTracksFiles(filepath_IALT);
	filepath_ET = [dir 'Tracks/ConvergedTracks/ExitTracks'];
	exit_tracks = ReadSeparateTracksFiles(filepath_ET);
else
	filepath_IALT = [dir 'Tracks/ConvergedTracks/InactiveLongTracks' num2str(frame_num) '.txt'];
	inactive_long_tracks = ReadTracks(filepath_IALT);
	filepath_ET = [dir 'Tracks/ConvergedTracks/ExitTracks' num2str(frame_num) '.txt'];
        exit_tracks = ReadTracks(filepath_ET);
end
tracks = PackTracks(active_long_tracks, inactive_long_tracks, exit_tracks);
VolumeSelfCalib(tracks, 66, dir, 1, calibration_file, save_name); 
% addpath OTF_VSC_Version102717/;
% x_range = [-30, 30];
% y_range = [-30, 30];
% z_range = [-20, 20];
% VSCOTF(tracks, x_range, y_range, z_range, 100, dir, [dir save_name '.mat']);
end


