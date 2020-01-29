function [a,b,c,alpha] = OTFParametersmesh(n,particlesubvolumedata,subvol_OTFParameters,fieldname_cam,ncam)

% This will initialize a n*n*n array for each parameter with name of ncam
for i =1:length(fieldname_cam)
    type.(fieldname_cam{i}) = zeros(1,n*n*n);
end
a = repmat(type,1);
b = repmat(type,1);
c = repmat(type,1);
alpha = repmat(type,1);
subnum=numel(particlesubvolumedata);

for i =1:ncam
    for j = 1:subnum 
        nth=particlesubvolumedata(j).Xsub+(particlesubvolumedata(j).Ysub-1)*n+(particlesubvolumedata(j).Zsub-1)*n*n; %nth=xsub(i)+(ysub(i)-1)*number of elements in x-direction ...
                                                    ...+(zsub(i)-1)*number of elements in x-direction*number of elements in y-direction;
        a.(fieldname_cam{i})(1,nth)=subvol_OTFParameters(j).(fieldname_cam{i})(1);
        alpha.(fieldname_cam{i})(1,nth)=subvol_OTFParameters(j).(fieldname_cam{i})(2);    
        b.(fieldname_cam{i})(1,nth)=subvol_OTFParameters(j).(fieldname_cam{i})(3);
        c.(fieldname_cam{i})(1,nth)=subvol_OTFParameters(j).(fieldname_cam{i})(4);
    end
        a.(fieldname_cam{i}) = permute(reshape(a.(fieldname_cam{i}),[n,n,n]),[2 1 3]);
        alpha.(fieldname_cam{i}) = permute(reshape(alpha.(fieldname_cam{i}),[n,n,n]),[2 1 3]);
        b.(fieldname_cam{i}) = permute(reshape(b.(fieldname_cam{i}),[n,n,n]),[2 1 3]);
        c.(fieldname_cam{i}) = permute(reshape(c.(fieldname_cam{i}),[n,n,n]),[2 1 3]);
end
