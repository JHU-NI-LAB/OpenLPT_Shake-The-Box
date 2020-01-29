function OTFpara = OTFaverage(nframe,ncam,outputfiletitle,averageOTFfilename)


fieldname= cell(1,nframe); % create an array of strings for loading the OTFParameters from different frames
framename = cell(1,nframe);
camname = cell(1,ncam);
for i = 1:ncam
    int=num2str(i);
    camname{i} = strcat('cam',int);
end
alpha=[];%Since you're assigning a value to alpha, MATLAB knows that it should treat alpha as a variable in this function.
for i = 1:nframe
    int=num2str(i);
    fieldname{i} = strcat(outputfiletitle,int,'.mat');
    framename{i}=strcat('frame',int);
    load(fieldname{i});
    A.(framename{i})= a;         
    B.(framename{i})= b;     
    C.(framename{i})= c; 
    Alpha.(framename{i})= alpha; 
end

% To call the 'a' parameter from cam 1 for the second frame, type
% a.frame2.cam1

[nx,ny,nz]=size(X_map);

for i = 1:ncam
    type.(camname{i})=zeros(nx,ny,nz); 
end
a_ave =repmat(type,1);
alpha_ave =repmat(type,1);
b_ave =repmat(type,1);
c_ave =repmat(type,1);
a_temp = zeros(1,nframe);
b_temp = zeros(1,nframe);
c_temp = zeros(1,nframe);
alpha_temp = zeros(1,nframe);
tic
for m = 1:ncam
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                for f = 1:nframe %frame number     
                    a_temp(f) = A.(framename{f}).(camname{m})(i,j,k);
                    b_temp(f) = B.(framename{f}).(camname{m})(i,j,k);
                    c_temp(f) = C.(framename{f}).(camname{m})(i,j,k);
                    alpha_temp(f) = Alpha.(framename{f}).(camname{m})(i,j,k);
                     if f == nframe 
                            switch sum(a_temp)
                                case 0
                                a_ave.(camname{m})(i,j,k)=0;
                                a_temp = zeros(1,nframe); 
                                otherwise
                                    a_temp = nonzeros(a_temp);
                                    outliers = isoutlier(a_temp,'gesd');
                                    a_temp(outliers) = []; %delete outlier
                                a_ave.(camname{m})(i,j,k) = sum(a_temp)/nnz(a_temp);              
                                a_temp = zeros(1,nframe);
                            end

                            switch sum(b_temp)
                                case 0 
                                b_ave.(camname{m})(i,j,k)=0;
                                b_temp = zeros(1,nframe);
                                otherwise
                                    b_temp = nonzeros(b_temp);
                                    outliers = isoutlier(b_temp,'gesd');
                                    b_temp(outliers) = []; %delete outlier
                                b_ave.(camname{m})(i,j,k) = sum(b_temp)/nnz(b_temp);
                                b_temp = zeros(1,nframe);
                            end

                            switch sum(c_temp)
                                case 0
                                c_ave.(camname{m})(i,j,k)=0;
                                c_temp = zeros(1,nframe);
                                otherwise
                                    c_temp = nonzeros(c_temp);
                                    outliers = isoutlier(c_temp,'gesd');
                                    c_temp(outliers) = []; %delete outlier
                                c_ave.(camname{m})(i,j,k) = sum(c_temp)/nnz(c_temp);
                                c_temp = zeros(1,nframe);
                            end
                            switch sum(alpha_temp)
                                case 0
                                alpha_ave.(camname{m})(i,j,k)=0;
                                alpha_temp = zeros(1,nframe);
                                otherwise
                                    alpha_temp = nonzeros(alpha_temp);
                                    outliers = isoutlier(alpha_temp,'gesd');
                                    alpha_temp(outliers) = []; %delete outlier
                                alpha_ave.(camname{m})(i,j,k) = sum(alpha_temp)/nnz(alpha_temp);  
                                alpha_temp = zeros(1,nframe);
                            end    
                    end
                end
            end
        end
    end
end
fprintf('Time taken to average the values: %f\n', toc)

%% Make the variables into 1D for C++
for i = 1:ncam
    type2.(camname{i})=zeros(1,nx*ny*nz); 
end
A =repmat(type2,1);
B =repmat(type2,1);
C =repmat(type2,1);
Alpha =repmat(type2,1);

for i=1:ncam
    A.(camname{i})  = transpose_reshape(a_ave.(camname{i}));
    B.(camname{i})  = transpose_reshape(b_ave.(camname{i}));
    C.(camname{i})  = transpose_reshape(c_ave.(camname{i}));
    Alpha.(camname{i})  =   transpose_reshape(alpha_ave.(camname{i}));
    X_trans = transpose_reshape(X_map);
    Y_trans = transpose_reshape(Y_map);
    Z_trans = transpose_reshape(Z_map);
end
    
    
%% Output text file

% Variables' names
Aname = cell(1,ncam);
Bname = cell(1,ncam);
Cname = cell(1,ncam);
Alphaname = cell(1,ncam);
for i = 1:ncam
    int=num2str(i);
    Aname{i} = strcat('#A',int);
    Bname{i} = strcat('#B',int);
    Cname{i} = strcat('#C',int);
    Alphaname{i} = strcat('#Alpha',int);
end

fileID = fopen(averageOTFfilename,'w');
[~,nelements]=size(X_trans);
fprintf(fileID,'%i # number of elements \r\n', nelements);
% fprintf(fileID,'%i # number of xyz 1D vector elements \r\n', numel(X));
nparameters=7;
for i = 1:nparameters
    switch i
        case 1
            for j = 1:ncam
                fprintf(fileID,'%f\r\n',A.(camname{j}));
                fprintf(fileID,'%s\r\n\r\n',Aname{j});
            end
        case 2
            for j = 1:ncam
                fprintf(fileID,'%f\r\n',B.(camname{j}));
                fprintf(fileID,'%s\r\n\r\n',Bname{j});
            end
        case 3
            for j = 1:ncam
                fprintf(fileID,'%f\r\n',C.(camname{j}));
                fprintf(fileID,'%s\r\n\r\n',Cname{j});
            end
        case 4
            for j = 1:ncam
                fprintf(fileID,'%f\r\n',Alpha.(camname{j}));
                fprintf(fileID,'%s\r\n\r\n',Alphaname{j});
            end
        case 5
            fprintf(fileID,'%g\r\n',X);
            fprintf(fileID,'%s\r\n\r\n','#X');
        case 6
            fprintf(fileID,'%f\r\n',Y);
            fprintf(fileID,'%s\r\n\r\n','#Y');          
        case 7
            fprintf(fileID,'%f\r\n',Z);
            fprintf(fileID,'%s\r\n\r\n','#Z');          
            
    end
        
end
    
fclose(fileID);

OTFpara = {a_ave; b_ave; c_ave; alpha_ave; X; Y; Z};

% % Create text file for XYZ subvolume vectors
% fileID = fopen(XYZsubvolname,'w');
% fprintf(fileID,'%i # number of xyz 1D vector elements \r\n', numel(X));
% for i = 1:3
%     switch i
%         case 1
%             fprintf(fileID,'%g\r\n',X);
%             fprintf(fileID,'%s\r\n\r\n','#X');
%         case 2
%             fprintf(fileID,'%f\r\n',Y);
%             fprintf(fileID,'%s\r\n\r\n','#Y');          
%         case 3
%             fprintf(fileID,'%f\r\n',Z);
%             fprintf(fileID,'%s\r\n\r\n','#Z');
%     end
% end
% fclose(fileID);
%     
