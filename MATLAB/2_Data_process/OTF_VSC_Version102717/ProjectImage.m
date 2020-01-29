function img = ProjectImage(para, xycenter)
xycenter = round(xycenter);
pixel_range = 2;
for i = xycenter(1) - pixel_range : xycenter(1) + pixel_range
    for j = xycenter(2) - pixel_range : xycenter(2) + pixel_range
        img(j - (xycenter(2) - pixel_range) + 1, i - (xycenter(1) - pixel_range) + 1) = Projection(para, [i, j]);
    end
end
end

