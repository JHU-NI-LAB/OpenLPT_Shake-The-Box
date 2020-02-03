function OTFParam = ExtrapolateOTF(OTFParam)
%there are some grids that no points is found, and therefore needs
%extrapolation
[~, n, ~, ~] = size(OTFParam{1});
%[X,Y,Z] = meshgrid(midpoints,midpoints,midpoints);

% the extrapolation is done from inside to outside
% for l = 1 : 4
    % find the minimun nonzeros cube
%     for ncam = 1 : 4
%         ind = find(reshape(OTFParam{l}(ncam,:,:,:), [n, n, n]) ~= 0);
%         [I1,I2,I3] = ind2sub([n, n, n], ind);
%         n_min = max([I1(1),I2(1),I3(1)]);
%         n_max = min([I1(end),I2(end),I3(end)]);
%         nn = n_max - n_min + 1;
%         new_midpoints = midpoints(n_min : n_max);
%         [X,Y,Z] = meshgrid(new_midpoints,new_midpoints,new_midpoints);
        % expand by 1 until reach the edge
%         while (n_min > 1 || n_max < n)
%             n_min1 = max(1, n_min - 1); n_max1 = min(n, n_max +1);
%             for i = n_min1 : n_max1
%                 for j = n_min1 : n_max1
%                     for k = n_min1 : n_max1
%                         if (OTFParam{l}(ncam, i, j, k) == 0)
%                             OTFParam{l}(ncam, i, j, k) = interp3(X, Y, Z, ...
%                                 reshape(OTFParam{l}(ncam,n_min : n_max,n_min : n_max,n_min : n_max), [nn, nn, nn]), ...
%                                 midpoints(i), midpoints(j), midpoints(k), 'makima');
%                         end
%                     end
%                 end
%             end
%             n_min = n_min1; n_max = n_max1; nn = n_max - n_min + 1;
%         end
        
%     end
% end

% neighbor index
% x = [-1 0 1];
% [X,Y,Z] = meshgrid(x,x,x);
% X = reshape(X, [27,1]);
% Y = reshape(Y, [27,1]);
% Z = reshape(Z, [27,1]);
% offset = [X, Y, Z];
offset = [0 -1 0; -1 0 0; 0 0 -1; 1 0 0;0 0 1; 0 1 0];

zero_flag = 1; 
while (zero_flag)
zero_flag = 0;
for l = 1 : 4
    for ncam = 1 : 4
        for i = 1 : n
            for j = 1 : n
                for k = 1 :n
                    if (OTFParam{l}(ncam, i, j, k) == 0) 
                        neighbors = [i, j, k] + offset;
                        neighbors(neighbors(:, 1) < 1 | neighbors(:, 1) > n | ...
                            neighbors(:, 2) < 1 | neighbors(:, 2) > n | ...
                            neighbors(:, 3) < 1 | neighbors(:, 3) > n, :) = [];
                        
                        num_neighbors = size(neighbors, 1);
                        para = zeros(num_neighbors, 1);
                        for m = 1 : num_neighbors
                            para(m) = OTFParam{l}(ncam, neighbors(m, 1),neighbors(m, 2),neighbors(m, 3));
                        end
                        OTFParam{l}(ncam, i, j, k) =  mean(nonzeros(para));
                        if (isnan(OTFParam{l}(ncam, i, j, k))) 
                            OTFParam{l}(ncam, i, j, k) = 0;
                            zero_flag = zero_flag | 1;
                        end
                    end
                end
            end
        end
    end
end

end

end



