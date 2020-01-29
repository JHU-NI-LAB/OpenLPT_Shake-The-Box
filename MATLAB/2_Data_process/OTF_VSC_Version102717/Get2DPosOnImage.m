function position2D = Get2DPosOnImage(img) 
[rows, cols] = size(img);
threshold = 30;
position2D = [];
for i = 2 : rows - 1
    for j = 2 : cols - 1
        if (img(i, j) >= threshold && IsLocalMax(img, i, j))
            x1 = j - 1; x2 = j; x3 = j + 1;
            y1 = i - 1; y2 = i; y3 = i + 1;
            ln1 = NoInfLog(img(i, j - 1));
            ln2 = NoInfLog(img(i, j));
            ln3 = NoInfLog(img(i, j + 1));
        
            xc = -.5 * (ln1 * (x2 ^ 2 - x3 ^ 2) - ln2 * (x1 ^ 2 - x3 ^ 2) + ln3 * (x1 ^ 2 - x2 ^ 2)) / ...
            (ln1 * (x3 - x2) - ln3 * (x1 - x2) + ln2 * (x1 - x3));
            
            ln1 = NoInfLog(img(i - 1, j));
            ln2 = NoInfLog(img(i, j));
            ln3 = NoInfLog(img(i + 1, j));
            yc = -.5 * (ln1 * (y2 ^ 2 - y3 ^ 2) - ln2 * (y1 ^ 2 - y3 ^ 2) + ln3 * (y1 ^ 2 - y2 ^ 2)) / ...
            (ln1 * (y3 - y2) - ln3 * (y1 - y2) + ln2 * (y1 - y3));
            
            if (~isinf(xc) && ~isinf(yc))
                position2D = [position2D; xc, yc];
            end
        end
    end
end

end

function result = IsLocalMax(img, i, j)
   result = 1;
   if (img(i - 1, j) > img(i, j) || img(i + 1, j) > img(i, j) || ...
           img(i, j - 1) > img(i, j) || img(i, j + 1) > img(i, j))
       result = 0;
   end
end

function y = NoInfLog(x)
if x == 0
    x = .0001;
end
y = log(double(x));
end