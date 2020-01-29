function makeCamConfigFileV2(camParaCalib,fname)
% Expects a 1xN struct coming out of calibTsai.
% Writes a text file that can be read by the 3D tracking code.

fid=fopen(fname,'w');

fprintf(fid,'# Camera configuration file\n');
fprintf(fid,['# Generated ' date '\n\n']);

s=size(camParaCalib);

fprintf(fid,'%d\t\t# number of cameras\n\n',s(2));

for ii=1:s(2)
  fprintf(fid,'######## camera %d ########\n',ii-1);
  fprintf(fid,'%d\t\t\t# Noffh\n',camParaCalib(ii).Noffh);
  fprintf(fid,'%d\t\t\t# Noffw\n',camParaCalib(ii).Noffw);
  fprintf(fid,'%d\t\t\t# Npix_x\n',camParaCalib(ii).Npixw);
  fprintf(fid,'%d\t\t\t# Npix_y\n',camParaCalib(ii).Npixh);
  fprintf(fid,'%g\t\t# pixsize_x\n',camParaCalib(ii).wpix);
  fprintf(fid,'%g\t\t# pixsize_y\n',camParaCalib(ii).hpix);
  fprintf(fid,'%g\t\t# f_eff\n',camParaCalib(ii).f_eff);
  fprintf(fid,'%g\t\t# kr\n',camParaCalib(ii).k1);
  %fprintf(fid,'%g\t\t# kx\n',camParaCalib(ii).k1star);
  fprintf(fid,'%g\t\t# kx\n',1);
  %kx is the translational distortion, which is 1. This is not available
  %for the current calibration. 
  
  % possible that this is backwards (rows v cols)...
  fprintf(fid,'%g\t\t# R\n',camParaCalib(ii).R');
  fprintf(fid,'%g\t\t# T\n',camParaCalib(ii).T);
  fprintf(fid,'%g\t\t# Rinv\n',camParaCalib(ii).Rinv');
  fprintf(fid,'%g\t\t# Tinv\n',camParaCalib(ii).Tinv);
  fprintf(fid,'\n');
end

fprintf(fid,'##### 3D match parameters #####\n\n');
fprintf(fid,'# mindist_pix\n');
fprintf(fid,'# mindist_3D\n');

disp('Don''t forget to specify the 3D matching parameters!');

