function img = ProjectBy3DPosV2(pos3D, OTFpara, campara, midpoints)
%midpoints = [-20 0 20];
[~, nx, ny, nz] = size(OTFpara{1});
[X,Y,Z] = meshgrid(midpoints,midpoints,midpoints);
a(1) = interp3(X, Y, Z, reshape(OTFpara{1}(1,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
a(2) = interp3(X, Y, Z, reshape(OTFpara{1}(2,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
a(3) = interp3(X, Y, Z, reshape(OTFpara{1}(3,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
a(4) = interp3(X, Y, Z, reshape(OTFpara{1}(4,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));

b(1) = interp3(X, Y, Z, reshape(OTFpara{2}(1,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
b(2) = interp3(X, Y, Z, reshape(OTFpara{2}(2,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
b(3) = interp3(X, Y, Z, reshape(OTFpara{2}(3,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
b(4) = interp3(X, Y, Z, reshape(OTFpara{2}(4,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));

c(1) = interp3(X, Y, Z, reshape(OTFpara{3}(1,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
c(2) = interp3(X, Y, Z, reshape(OTFpara{3}(2,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
c(3) = interp3(X, Y, Z, reshape(OTFpara{3}(3,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
c(4) = interp3(X, Y, Z, reshape(OTFpara{3}(4,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));

alpha(1) = interp3(X, Y, Z, reshape(OTFpara{4}(1,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
alpha(2) = interp3(X, Y, Z, reshape(OTFpara{4}(2,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
alpha(3) = interp3(X, Y, Z, reshape(OTFpara{4}(3,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));
alpha(4) = interp3(X, Y, Z, reshape(OTFpara{4}(4,:,:,:), [nx, ny, nz]), pos3D(2), pos3D(1), pos3D(3));

for i = 1 : 4
    xycenter(1:2) = calibProj(campara(i), pos3D);
    img(i, :, :) = ProjectImage([a(i) alpha(i) b(i) c(i) xycenter(1) xycenter(2)], xycenter);
end

end

