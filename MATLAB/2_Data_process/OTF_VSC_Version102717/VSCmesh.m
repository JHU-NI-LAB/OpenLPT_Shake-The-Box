function [Isum_mesh] = VSCmesh(n,particlesubvolumedata,Isub,fieldname_cam,ncam,exsize,eysize)

% This will initialize a n*n*n array for each parameter with name of ncam
for i = 1:length(fieldname_cam)
    type.(fieldname_cam{i}) = cell(1,n*n*n);
end
Isum_mesh = repmat(type,1);
emptymatrix=repmat(uint8(0),exsize,eysize); %Preallocate

for i = 1:ncam
    for j = 1:n*n*n
        Isum_mesh.(fieldname_cam{i}){1,j}=emptymatrix;
    end
end

subnum=numel(particlesubvolumedata);

for i =1:ncam
    for j = 1:subnum 
        nth=particlesubvolumedata(j).Xsub+(particlesubvolumedata(j).Ysub-1)*n+(particlesubvolumedata(j).Zsub-1)*n*n; %nth=xsub(i)+(ysub(i)-1)*number of elements in x-direction ...
                                                    ...+(zsub(i)-1)*number of elements in x-direction*number of elements in y-direction;
        Isum_mesh.(fieldname_cam{i}){1,nth}=Isub(j).(fieldname_cam{i});
    end
        Isum_mesh.(fieldname_cam{i}) = permute(reshape(Isum_mesh.(fieldname_cam{i}),[n,n,n]),[2 1 3]);
end
