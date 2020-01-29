function [subvol_partdispar,num,temp,particlesubvolumedata] = OTFParam(x_size,y_size,ncam,x_range,y_range,z_range,n,pos3D_matrix,imagefilepath,outputfilename,nthframe)

particlespos=[pos3D_matrix(:,1) pos3D_matrix(:,2) pos3D_matrix(:,3)]; % This will take particles 3D position

% Nx=x/n;
% Ny=y/n;
% Nz=z/n;
dx = (x_range(2) - x_range(1)) / (n - 1); % devided by n - 1 in order to make the edge as the center
dy = (y_range(2) - y_range(1)) / (n - 1);
dz = (z_range(2) - z_range(1)) / (n - 1);

% Define the center of each subvolumes
% x_subvolume=x/2;
% y_subvolume=y/2;
% z_subvolume=z/2;
% X = linspace(-x_subvolume,x_subvolume,n);
% Y = linspace(-y_subvolume,y_subvolume,n);
% Z = linspace(-z_subvolume,z_subvolume,n); 
X = x_range(1) : dx : x_range(2);
Y = y_range(1) : dy : y_range(2);
Z = z_range(1) : dz : z_range(2);
% Construct logic lists of particles that are detected in each subvolumes
particlelogicallist=cell(n,n,n); 
S_vertnface=cell(n,n,n);

%% Subvolume checking for particles
% The volume is divided into ixjxk number of subvolumes and each subvolumes is checked one by one from top to the bottom

for k=1:n
    for j=1:n
        for i =1:n
            [in,S] = particledetector(X(i),dx,Y(j),dy,Z(k),dz,particlespos);
            particlelogicallist{i,j,k}=in;
            S_vertnface{i,j,k}=S;
        end 
    end
end



%% Store actual 3D position of particles and xyz coordinate of the subvolume in the volume
% Xsub, Ysub and Zsub are the coordinates for subvolumes that have
% particles 
% type1=struct('Xsub',[],'Ysub',[],'Zsub',[],'X',[],'Y',[],'Z',[]);
fieldname = cell(1,6+2*ncam); % Allocate size for X and Y fieldnames for each of the cameras
fieldname{1}='Xsub'; fieldname{2}='Ysub'; fieldname{3}='Zsub'; 
fieldname{4}='X'; fieldname{5}='Y'; fieldname{6}='Z'; 

% Construct strings for X and Y fieldnames, e.g. X1,Y1
for k = 1:ncam
    int=num2str(k);
    fieldname{6+1+2*(k-1)} = strcat('X',int);
    fieldname{6+2*k} = strcat('Y',int);
end
% Construct a struct with different fieldnames and preallocated to be blank
for ii =1:length(fieldname)
    type1.(fieldname{ii}) = [];
end


% XYZ_particlescoordinates gives the positions for the particles in each of
% the subvolumes
particlesubvolumedata=repmat(type1,n^3,1);
num=0;

% ONLY the subvolumes with particles inside will be stored
for k=1:n
    for j=1:n
        for i =1:n
            if (sum(particlelogicallist{i,j,k})> 0)
                num=num+1;
                particlesubvolumedata(num).Xsub=i;
                particlesubvolumedata(num).Ysub=j;
                particlesubvolumedata(num).Zsub=k;
                particlesubvolumedata(num).X=particlespos(particlelogicallist{i,j,k},1);
                particlesubvolumedata(num).Y=particlespos(particlelogicallist{i,j,k},2);
                particlesubvolumedata(num).Z=particlespos(particlelogicallist{i,j,k},3);
            else 
                continue
            end
        end 
    end 
end 
% XYZ_subvolumecoordinates = XYZ_subvolumecoordinates(1:num, :);
particlesubvolumedata = particlesubvolumedata(1:num);

%% Variable Storing 
for i=1:num
    [a , ~]= size(particlesubvolumedata(i).X); %number of particles in each subvolume
    for d=1:a
        for j = 1:2*ncam 
          % This will find the index for x position from the pos3D_matrix..
           ...so other variables like X1 and Y1 can be called and stored
          % d is number of particles. i is the number of subvolumes
          idx = pos3D_matrix(:,1)== particlesubvolumedata(i).X(d) & ...
              pos3D_matrix(:,2)== particlesubvolumedata(i).Y(d) & ...
              pos3D_matrix(:,3)== particlesubvolumedata(i).Z(d);
%           if ~isempty(idx)
            particlesubvolumedata(i).(fieldname{6+j})(d)= pos3D_matrix(idx,j+3); %'j+3'because X1 is the 4th element in pos3D_matrix
%           end
        end
    end
end
% Example: If you want to call the X1 coordinate of the 1st particle in 2nd
% subvolume on the list, type 'particlesubvolumedata(2).X1(1)'

%% Image Input as A (A is the actual image (temp.tif))
tiff = cell(1,ncam);
imgread = cell(1,ncam);
tempmap = cell(1,ncam);
for k = 1:ncam
%     j = k;
%     if (j == 3)
%         k = 4;
%     end
%     if (j == 4)
%         k = 3;
%     end
    tiff{k} = [imagefilepath 'cam' num2str(k) '/cam' num2str(k) 'frame' num2str(nthframe,'%05d') '.tif'];
    imgread{k}=['cam' num2str(k) 'frame' num2str(nthframe,'%05d')];
end

for ind = 1:length(tempmap)
    [A.(imgread{ind}),map.(imgread{ind})]=imread(tiff{ind});

end
% Example: If you want to call the temp1 image, type 'A.temp1'.

%% 2D Gaussian fitters to get OTF Parameters (amplitude, angle, b, c, xcenter, ycenter)
fieldname_cam = cell(1,ncam);

for k = 1:ncam
    int=num2str(k);
    fieldname_cam{k} = strcat('cam',int);
end

for ii =1:length(fieldname_cam)
    type3.(fieldname_cam{ii}) = [];
end

subvol_particlesOTF =repmat(type3,num,1);
subvol_partdispar=repmat(type3,num,1);

temp=zeros(1,num); % number of particles in each subvolume
tic

% h = waitbar(0, 'Please wait...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); 
% setappdata(h,'canceling',0)
for i=1:num
    % Check for Cancel button press
%     if getappdata(h,'canceling')
%         break
%     end
    % computations take place here
%     waitbar(i/num) 
    [temp(i),~] = size(particlesubvolumedata(i).X);
    for d=1:temp(i)    
        for j = 1:ncam 
          [subvol_particlesOTF(i).(fieldname_cam{j}){d},subvol_partdispar(i).(fieldname_cam{j}){d}] = ...
              matrixstoring_2DG(particlesubvolumedata(i).(fieldname{6+1+2*(j-1)})(d),...
              particlesubvolumedata(i).(fieldname{6+2*j})(d),x_size,y_size,A.(imgread{j}),i);           
        end
    end
end
% delete(h)
fprintf('Calculation time for obtaining OTF Parameters of frame %i: %.4f sec\n',nthframe,toc);


%% Output OTFParameters for comparison
% size = sum(temp);
% OTFParameters_test=zeros(size,8);
% m=0;
% for i = 1:num
%     for p = 1:temp(i)    
% %       for j = 1:4
% %         if j == 1
%          m=m+1;
%          a = subvol_particlesOTF(i).cam1{1,p}(1);
%          alpha = subvol_particlesOTF(i).cam1{1,p}(2);
%          b = subvol_particlesOTF(i).cam1{1,p}(3);
%          c = subvol_particlesOTF(i).cam1{1,p}(4);
%          xo = subvol_particlesOTF(i).cam1{1,p}(5);
%          yo = subvol_particlesOTF(i).cam1{1,p}(6);
%          intensity = double(A.temp1(round(yo),round(xo)));
% %          intensity_test = double(A_test.tmp1(round(yo),round(xo)));
%          OTFParameters_test(m,1)=a;
%          OTFParameters_test(m,2)=alpha;
%          OTFParameters_test(m,3)=b;
%          OTFParameters_test(m,4)=c;
%          OTFParameters_test(m,5)=xo;
%          OTFParameters_test(m,6)=yo;
%          OTFParameters_test(m,7)=intensity;
% %          OTFParameters_test(m,8)=intensity_test;
% %         end
% %         if j == 2
% %         a = subvol_particlesmatrix(n).cam2{1,p}(1);
% %         alpha = subvol_particlesmatrix(n).cam2{1,p}(2);
% %         b = subvol_particlesmatrix(n).cam2{1,p}(3);
% %         c = subvol_particlesmatrix(n).cam2{1,p}(4);
% %         end
% %         if j == 3
% %         a = subvol_particlesmatrix(n).cam3{1,p}(1);
% %         alpha = subvol_particlesmatrix(n).cam3{1,p}(2);
% %         b = subvol_particlesmatrix(n).cam3{1,p}(3);
% %         c = subvol_particlesmatrix(n).cam3{1,p}(4);
% %         end
% %         if j == 4
% %         a = subvol_particlesmatrix(n).cam4{1,p}(1);
% %         alpha = subvol_particlesmatrix(n).cam4{1,p}(2);
% %         b = subvol_particlesmatrix(n).cam4{1,p}(3);
% %         c = subvol_particlesmatrix(n).cam4{1,p}(4);
% %         end
% %       end
%     end
% end


%% Averaging parameters in each subvolume 
subvol_OTFParameters =repmat(type3,num,1);
% amp_subvol=0; ang_subvol=0; b_subvol=0; c_subvol =0;
for i=1:num
    for j = 1:ncam  
        amp_temp = zeros(1,temp(i)); 
        ang_temp = zeros(1,temp(i));
        b_temp = zeros(1,temp(i));
        c_temp = zeros(1,temp(i));
        for d=1:temp(i) 
%           amp_temp = subvol_particlesOTF(i).(fieldname_cam{j}){1,d}(1);
%           ang_temp = subvol_particlesOTF(i).(fieldname_cam{j}){1,d}(2);
%           b_temp   = subvol_particlesOTF(i).(fieldname_cam{j}){1,d}(3);
%           c_temp   = subvol_particlesOTF(i).(fieldname_cam{j}){1,d}(4);
%           amp_subvol = amp_subvol + amp_temp;
%           ang_subvol = ang_subvol + ang_temp;
%           b_subvol = b_subvol + b_temp;
%           c_subvol = c_subvol + c_temp;   
%           if d == temp(i)
%              % Calculating the average of each parameters     
%              subvol_OTFParameters(i).(fieldname_cam{j})(1) = amp_subvol/d;
%              subvol_OTFParameters(i).(fieldname_cam{j})(2) = ang_subvol/d;
%              subvol_OTFParameters(i).(fieldname_cam{j})(3) = b_subvol/d;
%              subvol_OTFParameters(i).(fieldname_cam{j})(4) = c_subvol/d;
%              %Reinitializing the sum of each parameters to be 0
%              amp_subvol = 0;
%              ang_subvol = 0;
%              b_subvol = 0;
%              c_subvol = 0;  
            amp_temp(d) = subvol_particlesOTF(i).(fieldname_cam{j}){1,d}(1);
            ang_temp(d) = subvol_particlesOTF(i).(fieldname_cam{j}){1,d}(2);
            b_temp(d)   = subvol_particlesOTF(i).(fieldname_cam{j}){1,d}(3);
            c_temp(d)   = subvol_particlesOTF(i).(fieldname_cam{j}){1,d}(4);
            if d == temp(i)
                switch sum(amp_temp)
                    case 0
                    subvol_OTFParameters(i).(fieldname_cam{j})(1) = 0;
                    amp_temp = zeros(1,temp(i)); 
                    otherwise
                        outliers = isoutlier(amp_temp,'gesd');
                        amp_temp(outliers) = []; %delete outlier
                    subvol_OTFParameters(i).(fieldname_cam{j})(1) = sum(amp_temp)/nnz(amp_temp);              
                    amp_temp = zeros(1,temp(i));
                end

                switch sum(b_temp)
                    case 0 
                    subvol_OTFParameters(i).(fieldname_cam{j})(3)=0;
                    b_temp = zeros(1,temp(i));
                    otherwise
                        outliers = isoutlier(b_temp,'gesd');
                        b_temp(outliers) = []; %delete outlier
                    subvol_OTFParameters(i).(fieldname_cam{j})(3) = sum(b_temp)/nnz(b_temp);
                    b_temp = zeros(1,temp(i));
                end

                switch sum(c_temp)
                    case 0
                    subvol_OTFParameters(i).(fieldname_cam{j})(4)=0;
                    c_temp = zeros(1,temp(i));
                    otherwise
                        outliers = isoutlier(c_temp,'gesd');
                        c_temp(outliers) = []; %delete outlier
                    subvol_OTFParameters(i).(fieldname_cam{j})(4) = sum(c_temp)/nnz(c_temp);
                    c_temp = zeros(1,temp(i));
                end
                switch sum(ang_temp)
                    case 0
                    subvol_OTFParameters(i).(fieldname_cam{j})(2)=0;
                    ang_temp = zeros(1,temp(i));
                    otherwise
                        outliers = isoutlier(ang_temp,'gesd');
                        ang_temp(outliers) = []; %delete outlier
                    subvol_OTFParameters(i).(fieldname_cam{j})(2) = sum(ang_temp)/nnz(ang_temp);  
                    ang_temp = zeros(1,temp(i));
                end
            end
        end
    end
end




%% Create n x n x n matrices for x, y, z, a, b, c and alpha for the subvolume

[X_map, Y_map, Z_map] = meshgrid(X,Y,Z);
[a,b,c,alpha] = OTFParametersmesh(n,particlesubvolumedata,subvol_OTFParameters,fieldname_cam,ncam);
save(outputfilename,'a','alpha','b','c','X_map','Y_map','Z_map','X','Y','Z','-mat');
end
