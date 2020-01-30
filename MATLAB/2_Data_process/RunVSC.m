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
viewsize = max(tracks(:, 3)) - min(tracks(:,3)); 
VolumeSelfCalib(tracks, viewsize, dir, 1, calibration_file, save_name);
end


