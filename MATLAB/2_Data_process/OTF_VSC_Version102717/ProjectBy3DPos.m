function img = ProjectBy3DPos(pos3D, OTFpara, campara)
[X,Y,Z] = meshgrid(OTFpara{5},OTFpara{6},OTFpara{7});
a(1) = interp3(X, Y, Z, OTFpara{1}.cam1, pos3D(1), pos3D(2), pos3D(3));
a(2) = interp3(X, Y, Z, OTFpara{1}.cam2, pos3D(1), pos3D(2), pos3D(3));
a(3) = interp3(X, Y, Z, OTFpara{1}.cam3, pos3D(1), pos3D(2), pos3D(3));
a(4) = interp3(X, Y, Z, OTFpara{1}.cam4, pos3D(1), pos3D(2), pos3D(3));

b(1) = interp3(X, Y, Z, OTFpara{2}.cam1, pos3D(1), pos3D(2), pos3D(3));
b(2) = interp3(X, Y, Z, OTFpara{2}.cam2, pos3D(1), pos3D(2), pos3D(3));
b(3) = interp3(X, Y, Z, OTFpara{2}.cam3, pos3D(1), pos3D(2), pos3D(3));
b(4) = interp3(X, Y, Z, OTFpara{2}.cam4, pos3D(1), pos3D(2), pos3D(3));

c(1) = interp3(X, Y, Z, OTFpara{3}.cam1, pos3D(1), pos3D(2), pos3D(3));
c(2) = interp3(X, Y, Z, OTFpara{3}.cam2, pos3D(1), pos3D(2), pos3D(3));
c(3) = interp3(X, Y, Z, OTFpara{3}.cam3, pos3D(1), pos3D(2), pos3D(3));
c(4) = interp3(X, Y, Z, OTFpara{3}.cam4, pos3D(1), pos3D(2), pos3D(3));

alpha(1) = interp3(X, Y, Z, OTFpara{4}.cam1, pos3D(1), pos3D(2), pos3D(3));
alpha(2) = interp3(X, Y, Z, OTFpara{4}.cam2, pos3D(1), pos3D(2), pos3D(3));
alpha(3) = interp3(X, Y, Z, OTFpara{4}.cam3, pos3D(1), pos3D(2), pos3D(3));
alpha(4) = interp3(X, Y, Z, OTFpara{4}.cam4, pos3D(1), pos3D(2), pos3D(3));

for i = 1 : 4
    xycenter(1:2) = calibProj(campara(i), pos3D);
    img(i, :, :) = ProjectImage([a(i) alpha(i) b(i) c(i) xycenter(1) xycenter(2)], xycenter);
end

end

