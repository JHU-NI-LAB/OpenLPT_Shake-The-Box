function pos3D = DeleteOverlapParticle(pos3D, image_filepath, nthframe)
num_particle = size(pos3D, 1);
delete_index = zeros(1, num_particle);
ncams = 4;
pixle_error = 3;
yrange = 0; xrange = 0;
for i = 1 : ncams
    img = imread([image_filepath 'cam' num2str(i) '/cam' num2str(i) 'frame' num2str(nthframe,'%05d') '.tif']);
    posistion2D = Get2DPosOnImage(img);
    [yrange, xrange] = size(img);
    for j = 1 : num_particle
%         dist = vecnorm(pos3D(:, (i - 1) * 2 + 4 : (i - 1) * 2 + 5) - pos3D(j, (i - 1) * 2 + 4 : (i - 1) * 2 + 5), 2, 2);
        dist = vecnorm(posistion2D - pos3D(j, (i - 1) * 2 + 4 : (i - 1) * 2 + 5), 2, 2);
        candidate = posistion2D(dist == min(dist), :);
        pos3D(j, (i - 1) * 2 + 4 : (i - 1) * 2 + 5) = candidate(1,:); % use the 2d center to represent the 2D posisiton of the particle
        if sum(dist < pixle_error) > 1 % if the distant between two projection is less than 5 pixels
            delete_index(j) = 1;
        end
    end
end
pos3D(delete_index == 1, :) = [];
pos3D = unique(pos3D, 'rows'); % remove repeated particles
% Also delte particle at the edge
low_range = 2;
high_range = xrange - 2;
out_range = zeros(1, size(pos3D, 1))';
for i = 4 : 2 :11
    out_range = out_range | pos3D(:, i) <= low_range | pos3D(:, i) >= high_range;
end
pos3D(out_range == 1, :) =[];
high_range = yrange - 2;
out_range = zeros(1, size(pos3D, 1))';
for i = 5 : 2 :10
    out_range = out_range | pos3D(:, i) <= low_range | pos3D(:, i) >= high_range;
end
pos3D(out_range == 1, :) =[];
end

