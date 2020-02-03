function img = ImageProjection(pos3D, camPara, OTFPara, midpoints)
img = zeros(4, 1024, 1024);
proj_size = 2;
num_part = size(pos3D, 1);

for j = 1 : num_part
        temp = ProjectBy3DPosV2(pos3D(j, 1:3), OTFPara, camPara, midpoints);
        for i = 1 : 4
            xycenter = round(calibProj(camPara(i), pos3D(j, 1:3)));
            if (xycenter(2) - proj_size < 1 || xycenter(2) + proj_size > 1024 || ...
                   xycenter(1) - proj_size < 1 || xycenter(1) + proj_size > 1024 )
               continue;
            end
            img(i, xycenter(2) - proj_size : xycenter(2) + proj_size, xycenter(1) - proj_size : xycenter(1) + proj_size) = ...
                min(ones(1,5,5) * 255, img(i, xycenter(2) - proj_size : xycenter(2) + proj_size, xycenter(1) - proj_size : xycenter(1) + proj_size)+ ...
                temp(i, :, :));
        end
end

img = uint8(img);
end

