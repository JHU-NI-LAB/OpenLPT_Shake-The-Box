function tracks = ReadTracks(filepath)
fileID = fopen(filepath,'r');
formatspec = '%f,';
track_data = fscanf(fileID, formatspec);
total_len = length(track_data);
tracks = zeros(total_len / 5, 5);
for i = 1 :  total_len / 5
    tracks(i, :) = track_data((i - 1) * 5 + 1 : (i - 1) * 5 + 5);
end
end

