function finish = CheckSTB(dir)
lastline = GetLastTextLine(dir);
if contains(lastline, 'All of images have been processed!')
    finish = 1;
else
    finish =0;
end

end

function pline = GetLastTextLine(filepath)
fid = fopen(filepath);

while 1
    line = fgetl(fid);

    if ~ischar(line)
        break;
    end

    pline = line;
end
fclose(fid);

end

