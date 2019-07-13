function [Isum_mesh,exsize,eysize] = VSC_disparitymap(subvol_partdispar,ncam,num,n,particlesubvolumedata,temp,epsln_x,epsln_y,stepsize)

% %% Summing disparities in each subvolume 
% subvol_disparity =repmat(type3,num,1);
% dixsum_subvol=0; diysum_subvol=0;
% for i=1:num
%     for j = 1:ncam  
%         for d=1:temp(i) 
%           dix_temp = subvol_partdispar(i).(fieldname_cam{j}){1,d}(1);
%           diy_temp = subvol_partdispar(i).(fieldname_cam{j}){1,d}(2);
%           dixsum_subvol = dixsum_subvol + dix_temp;
%           diysum_subvol = diysum_subvol + diy_temp; 
%           if d == temp(i)
%              Calculating the average of each parameters     
%              subvol_disparity(i).(fieldname_cam{j})(1) = dixsum_subvol;
%              subvol_disparity(i).(fieldname_cam{j})(2) = diysum_subvol;
%              Reinitializing the sum of each parameters to be 0
%              dixsum_subvol = 0;
%              diysum_subvol = 0;
%           end
%         end
%     end
% end

%% Variables that define the disparity map
dx=-epsln_x:stepsize:epsln_x;
dy=-epsln_y:stepsize:epsln_y;
exsize=numel(dx);
eysize=numel(dy);

%% Project disparity map for each subvolume for all the camera
for i =1:ncam
    type.(['cam' num2str(i)]) = zeros(exsize,eysize);
end
Isub =repmat(type,num,1);
Isum=repmat(uint8(0),exsize,eysize); %Preallocate
for i=1:ncam
    for j = 1:num  
        for d=1:temp(j) 
            a = 200;	%Set an appropriate intensity amplitude
            alpha =1;
            b = 1;
            c = 1;
            dxx=subvol_partdispar(j).(['cam' num2str(i)]){1,d}(1);
            dyy=subvol_partdispar(j).(['cam' num2str(i)]){1,d}(2);
            xc = (exsize/2+0.5)+dxx/stepsize; % Matlab Disparity Position
            yc = (eysize/2+0.5)+dyy/stepsize; % Matlab Disparity Position
            Itemp = zeros(exsize,eysize);
            for x = 1:exsize
                   for y = 1:eysize
                        xx = (x-xc)*cos(alpha) + (y-yc)*sin(alpha);
                        yy = -(x-xc)*sin(alpha) + (y-yc)*cos(alpha);
                        Itemp(y,x) = max(Itemp(y,x), a*exp(-(b*(xx)^2 + c*(yy)^2))); % not sure if max is the right thing to use for overlapping particles
                   end
            end
            Itemp = uint8(Itemp);
            Isum=Isum+Itemp;
            if d == temp(j)
                 Isub(j).(['cam' num2str(i)])=Isum;  %Disparity map for each subvolume
                 Isum=repmat(uint8(0),exsize,eysize);
            end
        end
    end
end

%% Create n x n x n matrices for I of all the subvolumes
fieldname_cam = cell(1,ncam);
for k = 1:ncam
    int=num2str(k);
    fieldname_cam{k} = strcat('cam',int);
end
Isum_mesh = VSCmesh(n,particlesubvolumedata,Isub,fieldname_cam,ncam,exsize,eysize);
