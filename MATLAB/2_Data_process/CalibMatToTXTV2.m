function CalibMatToTXTV2(camParaCalib, n, save_path)
% n is the number of camera

fileID = fopen(save_path,'w');
fprintf(fileID, '# Camera configuration file\n');
fprintf(fileID, '# generated %s\n \n', datetime);

fprintf(fileID, [num2str(n) '    # camera number\n']);

for i = 1 : n
    fprintf(fileID, '\n#camera %d\n', i );
    fprintf(fileID,'%d    #Noffh\n',camParaCalib(i).Noffh);
    fprintf(fileID,'%d    #Noffw\n',camParaCalib(i).Noffw);
    fprintf(fileID,'%d    #Npixw\n',camParaCalib(i).Npixw);
    fprintf(fileID,'%d    #Npixh\n',camParaCalib(i).Npixh);
    fprintf(fileID,'%f    #wpix\n',camParaCalib(i).wpix);
    fprintf(fileID,'%f    #hpix\n',camParaCalib(i).hpix);
    fprintf(fileID,'%f    #f_eff\n',camParaCalib(i).f_eff);
    fprintf(fileID,'%f    #kr\n',camParaCalib(i).k1);
    fprintf(fileID,'%d    #kx\n',1);
    fprintf(fileID,'%f    #R\n',camParaCalib(i).R');
    fprintf(fileID,'%f    #T\n',camParaCalib(i).T);
    fprintf(fileID,'%f    #Rinv\n',camParaCalib(i).Rinv');
    fprintf(fileID,'%f    #Tinv\n',camParaCalib(i).Tinv);
end

fclose(fileID);
end


