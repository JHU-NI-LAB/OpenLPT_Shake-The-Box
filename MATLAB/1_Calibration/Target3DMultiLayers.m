function Target3DMultiLayers(calib_dirr,ncams,layer_index,central_index)
% This code is used to convert 2D target data to 3D
% This code is developed for Karuna's target only

for cam=1:ncams
    if exist([calib_dirr,'cam',num2str(cam),'_3D.txt'],'file')
        delete([calib_dirr,'cam',num2str(cam),'_3D.txt']);
    end
    fileid=fopen([calib_dirr,'cam',num2str(cam),'_3D.txt'],'w');
    for i=1:length(layer_index)
        calib=load([calib_dirr,'cam',num2str(cam),'_',num2str(layer_index(i)),'.txt']);
        calib=[calib,zeros(size(calib,1),1)];
        calib(:,4:5)=calib(:,3:4);
        calib(:,3)=(central_index-layer_index(i))*2.54/4;
        if exist([calib_dirr,'cam',num2str(cam),'_',num2str(layer_index(i)),'_3D.txt'],'file')
            delete([calib_dirr,'cam',num2str(cam),'_',num2str(layer_index(i)),'_3D.txt']);
        end
        fprintf(fileid,'%.4f,%.4f,%.4f,%.4f,%.4f\r\n',calib');
    end
    fclose(fileid);
end
