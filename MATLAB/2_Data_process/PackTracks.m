function tracks = PackTracks(active_long_tracks, inactive_long_tracks, exit_tracks)
len_activelong = max(active_long_tracks(:, 1));
len_inactivelong = max(inactive_long_tracks(:, 1));
inactive_long_tracks(:, 1) = inactive_long_tracks(:, 1) + len_activelong + 1;
exit_tracks(:, 1) = exit_tracks(:, 1) + len_activelong + len_inactivelong + 2;
tracks = [active_long_tracks; inactive_long_tracks; exit_tracks];
end