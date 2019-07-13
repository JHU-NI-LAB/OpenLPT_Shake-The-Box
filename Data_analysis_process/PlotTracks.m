function PlotTracks(tracks, fig, marker)
if ~exist('fig', 'var') fig = 1;  end
if ~exist('marker', 'var') marker = 'b.';  end
figure(fig);
plot3(tracks(:, 3), tracks(:, 4), tracks(:, 5), marker);
end

