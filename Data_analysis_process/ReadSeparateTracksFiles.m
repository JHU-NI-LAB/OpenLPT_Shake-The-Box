function tracks = ReadSeparateTracksFiles(file_path)
%filepath format: ~/InactiveLongTracks
list = dir([file_path '*.txt']); % get the list of file in the same name
num_file = size(list, 1);
tracks = [];
for i = 1 : num_file
    track = ReadTracks([list(i).folder '/' list(i).name]);
    tracks = CombineTracks(tracks, track);
end
end

function tracks = CombineTracks(track1, track2)
if isempty(track1) 
    tracks = track2;
else
    len_1 = max(track1(:, 1));
    track2(:, 1) = track2(:, 1) + len_1 + 1;
    tracks = [track1; track2];
end
end

