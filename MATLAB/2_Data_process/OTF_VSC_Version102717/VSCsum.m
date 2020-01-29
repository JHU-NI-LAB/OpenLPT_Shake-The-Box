function I_sumframe = VSCsum(nframe,ncam,I,exsize,eysize)

framename = cell(1,nframe);
camname = cell(1,ncam);
% Define the naming for the cameras based on the number of cameras 
for i = 1:ncam
    int=num2str(i);
    camname{i} = strcat('cam',int);
end
% Define the naming for the frame based on the number of frames
for i = 1:nframe
    int=num2str(i);
    framename{i}=strcat('frame',int);
end

[nx,ny,nz]=size(I.frame1.cam1);

for i = 1:ncam
    type.(camname{i})=cell(nx,ny,nz); 
end
I_sumframe =repmat(type,1);

emptymatrix=repmat(uint8(0),exsize,eysize); %Preallocate

% Preallocate empty matrix for different cameras
for m = 1:ncam
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                I_sumframe.(camname{m}){i,j,k}=emptymatrix;
            end
        end
    end
end


for m = 1:ncam
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                for f = 1:nframe %frame number 
                   % This will sum the intensity map of the disparity
                   % across each subvolume for different frames
                   Itemp = double(cell2mat(I.(framename{f}).(camname{m})(i,j,k)));
                   I_sumframe.(camname{m}){i,j,k}= double(cell2mat(I_sumframe.(camname{m})(i,j,k)))+Itemp;
                   Itemp = zeros(exsize,eysize);
                end
                   % The numbers were in double so that the number can add 
                   % up easily. Codes below will convert double to uint8 by
                   % adjusting the maximum intensity to be 255
                   Isum = I_sumframe.(camname{m}){i,j,k};
                   [M,~] = max(Isum(:)); 
                   I_sumframe.(camname{m}){i,j,k} = uint8(I_sumframe.(camname{m}){i,j,k}*255/M);% 255 is the max intensity value, this will scale the peak of disparity map according to the max intensity 
            end
        end
    end
end
fprintf('VSC disparity summation is completed.\n')

