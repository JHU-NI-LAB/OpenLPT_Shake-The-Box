function I_subpix = VSCsubdisparpeak(exsize,eysize,stepsize,ncam,nframe,VSC_xrange,VSC_yrange,I_sumframe)
framename = cell(1,nframe);
camname = cell(1,ncam);

for i = 1:ncam
    int=num2str(i);
    camname{i} = strcat('cam',int);
end

for i = 1:nframe
    int=num2str(i);
    framename{i}=strcat('frame',int);
end

[nx,ny,nz]=size(I_sumframe.cam1);

for i = 1:ncam
    type.(camname{i})=cell(nx,ny,nz); 
end

I_subpix =repmat(type,1);

for m = 1:ncam
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                I_subpix.(camname{m}){i,j,k}=[];
            end
        end
    end
end

for m = 1:ncam
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                Isub_sumframe =  I_sumframe.(camname{m})(i,j,k);
                dispar_brightest = VSCpeakfinder(exsize,eysize,stepsize,VSC_xrange,VSC_yrange,Isub_sumframe,m,i,j,k);
                I_subpix.(camname{m}){i,j,k}= dispar_brightest;
            end
        end
    end
end