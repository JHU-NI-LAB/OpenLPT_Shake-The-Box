function Target3DTwoLayers(calib_dirr,ncams)
% convert 2D target txt file to 3D
% This code is developed for the two-layer target

for cam=1:ncams
    calib=load([calib_dirr,'cam',num2str(cam),'.txt']);
    calib=[calib,zeros(size(calib,1),1)];
    for i=1:size(calib,1)
        if mod((calib(i,3)+calib(i,4)),2)
            calib(i,5)=-0.1;
        end
    end
    
    temp=calib(:,5);
    calib(:,4:5)=calib(:,3:4);
    calib(:,3)=temp;
    
    if exist([calib_dirr,'cam',num2str(cam),'_3D.txt'],'file')
        delete([calib_dirr,'cam',num2str(cam),'_3D.txt']);
    end
    
    fileid=fopen([calib_dirr,'cam',num2str(cam),'_3D.txt'],'w');
    fprintf(fileid,'%.4f,%.4f,%.4f,%.4f,%.4f\r\n',calib');
    fclose(fileid);
end
end