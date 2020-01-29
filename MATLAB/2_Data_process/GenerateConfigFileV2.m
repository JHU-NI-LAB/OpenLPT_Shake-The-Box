function GenerateConfigFileV2(dir, startframe, endframe, n, calibration_name)
% n is the number of cameras
GenerateTrackingConfig(dir, n, startframe, endframe);
GenerateIprConfig(dir, calibration_name);
GenerateOTF(dir, n);
GeneratePredictiveField(dir);
% create directory for saving tracks
system(['mkdir ' dir 'Tracks']);
system(['mkdir ' dir 'Tracks/InitialTracks']);
system(['mkdir ' dir 'Tracks/ConvergedTracks']);
system(['mkdir ' dir 'Tracks/BackSTBTracks']);
end

function GenerateTrackingConfig(dir, n, startframe, endframe)
fileID = fopen([dir 'trackingconfig1.txt'],'w');
txt = [ num2str(n) ' #number of cameras \n' ];
for i = 1 : n
    txt = [txt num2str(i) ' # ' num2str(i) 'th camera number \n'];
end
for i = 1 : n
    txt = [txt 'cam' num2str(i) 'ImageNames.txt # text file with' num2str(i) 'th camera image names \n'];
end
% txt = ['4 # Number of cameras \n'...
% '1 # first camera number \n'...
% '2 # second camera number \n'...
% '3 # third camera number \n'...
% '4 # fourth camera number \n'...
% 'cam1ImageNames.txt # text file with first camera image names \n'...
% 'cam2ImageNames.txt # text file with second camera image names \n'...
% 'cam3ImageNames.txt # text file with third camera image names \n'...
% 'cam4ImageNames.txt # text file with fourth camera image names \n'];
fprintf(fileID, txt);
fprintf(fileID, [dir 'iprconfig.txt # Path to ipr configuration file\n']);
fprintf(fileID, [dir 'predictivefield.txt # path to predictive field file\n']);
fprintf(fileID, [ num2str(startframe), ' # first frame\n', num2str(endframe), ' #last frame\n']);
txt =['./matched.gdf # stereomatched 3D positions \n'...
'./tracks.gdf # 3D tracks output filename \n'...
'########### View area limits ############ \n'...
'-30 # xmin \n'...
'30  # xmax \n'...
'-30 # ymin \n'...
'30  # xmax \n'...
'-25 # zmin \n'...
'25  # zmax \n'...
'	######### Initial Phase ############## \n'...
'1 # Flag for using ipr in initialphase (or use .mat files) \n'...
'10 # searchRadius for finding tracks using predictive field \n'...
'	######### Convergence Phase ############# \n'...
'.5 # Shaking range for prediciton (vox)\n' ...
'10 # Avg Interparticle spacing. (vox) to identify neighbouring tracks \n'...
'10 # Largest expected particle shift between frames (vox)for nearest neighbour linking of short tracks \n'...
'6 # Maximum absolute change in particle shift (vox) \n'...
'100  # Maximum relative change in particle shift (percent) \n'...
'2 # A multyplying factor on particle intensity in order to ensure the residual has no traces of tracked particles \n'...
'0.01 # lower intensity threshold (xx*avg. Intensity) to eliminate ghost particles while tracking \n'...
'0 # Back STB is on for 1\n' ...
'3.5 # the distance between two tracks that are supposed to be the same track \n'];
fprintf(fileID, txt);
fclose(fileID);
end

function GenerateIprConfig(dir, calibration_name)
fileID = fopen([dir 'iprconfig.txt'],'w');
txt = ['0 # Triangulation Only? (No IPR or tracking) 1 for yes/ 0 for no \n'...
'0 # Triangulation & IPR only? (No tracking)\n \n'];
fprintf(fileID, txt);
fprintf(fileID, [dir calibration_name ' # Path to camera calibration file\n']);
fprintf(fileID, [dir '# Path to TIFF files \n']);
fprintf(fileID, [dir 'OTFParameters.txt # Path to OTF text data file\n']);
txt = ['4  # Average particle size in pixels\n'...
'4  # No. of outerloop iterations\n'...
'4  # No. of innerloop iterations\n'...
'45 # 2D particle finder threshold\n'...
'8  # number of bits for each pixel\n'...
'0.01 # lower intensity threshold (xx*avg. Intensity) to eliminate ghost particles\n'...
'.75		# mindist_2D \n'...
'.75		# mindist_3D \n'...
'\n'...
'1 # use reduced cams (apply IPR by removing 1 cam each time)? 1 for yes/ 0 for no\n'...
'2 # no. of loops for each reduced camera combination\n'...
'.5 # mindist_2D for reduced cams\n'...
'.5 # mindist_3D for reduced cams\n'];
fprintf(fileID, txt);
fclose(fileID);
end

function GenerateOTF(dir, n)
fileID = fopen([dir 'OTFParameters.txt'],'w');
fprintf(fileID,'27 # number of elements \n');
for i = 1 : n
    A = ones(27, 1) * 255.000000;
    for j = 1 : 27
        fprintf(fileID, [num2str(A(j)) '\n']);
    end
    fprintf(fileID,['#A' num2str(i) '\n\n']);
end
% A1 = ones(27, 1) * 255.000000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(A1(i)) '\n']);
% end
% fprintf(fileID,'#A1\n\n');
% A2 = ones(27, 1) * 255.000000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(A2(i)) '\n']);
% end
% fprintf(fileID,'#A2\n\n');
% A3 = ones(27, 1) * 255.000000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(A3(i)) '\n']);
% end
% fprintf(fileID,'#A3\n\n');
% A4 = ones(27, 1) * 255.000000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(A4(i)) '\n']);
% end
% fprintf(fileID,'#A4\n\n');
for i = 1 : n
    B = ones(27, 1) * 0.400000;
    for j = 1 : 27
        fprintf(fileID, [num2str(B(j)) '\n']);
    end
    fprintf(fileID,['#B' num2str(i) '\n\n']);
end
% B1 = ones(27, 1) * 0.400000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(B1(i)) '\n']);
% end
% fprintf(fileID,'#B1\n\n');
% B2 = ones(27, 1) * 0.400000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(B2(i)) '\n']);
% end
% fprintf(fileID,'#B2\n\n');
% B3 = ones(27, 1) * 0.400000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(B3(i)) '\n']);
% end
% fprintf(fileID,'#B3\n\n');
% B4 = ones(27, 1) * 0.400000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(B4(i)) '\n']);
% end
% fprintf(fileID,'#B4\n\n');
for i = 1 : n
    C = ones(27, 1) * 0.400000;
    for j = 1 : 27
        fprintf(fileID, [num2str(C(j)) '\n']);
    end
    fprintf(fileID,['#C' num2str(i) '\n\n']);
end
% C1 = ones(27, 1) * 0.400000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(B1(i)) '\n']);
% end
% fprintf(fileID,'#C1\n\n');
% C2 = ones(27, 1) * 0.400000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(C2(i)) '\n']);
% end
% fprintf(fileID,'#C2\n\n');
% C3 = ones(27, 1) * 0.400000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(C3(i)) '\n']);
% end
% fprintf(fileID,'#C3\n\n');
% C4 = ones(27, 1) * 0.400000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(C4(i)) '\n']);
% end
% fprintf(fileID,'#C4\n\n');
for i = 1 : n
    alpha = ones(27, 1) * 0.000000;
    for j = 1 : 27
        fprintf(fileID, [num2str(alpha(j)) '\n']);
    end
    fprintf(fileID,['#alpha' num2str(i) '\n\n']);
end
% alpha1 = ones(27, 1) * 0.00000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(alpha1(i)) '\n']);
% end
% fprintf(fileID,'#alpha1\n\n');
% alpha2 = ones(27, 1) * 0.00000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(alpha2(i)) '\n']);
% end
% fprintf(fileID,'#alpha2\n\n');
% alpha3 = ones(27, 1) * 0.00000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(alpha3(i)) '\n']);
% end
% fprintf(fileID,'#alpha3\n\n');
% alpha4 = ones(27, 1) * 0.00000;
% for i = 1 : 27
%     fprintf(fileID, [num2str(alpha4(i)) '\n']);
% end
% fprintf(fileID,'#alpha4\n');
fclose(fileID);
end

function GeneratePredictiveField(dir)
fileID = fopen([dir 'predictivefield.txt'],'w');
txt = [
'########## grid sizes ###############\n'...
'12.5 # xgrid \n'...
'12.5 # ygrid\n'...
'12.5 # zgrid\n'...
'5 # searchRadius\n'... 
'\n'...
'########## mat file ###############\n'...
'0 # Flag for geting predictive field from .mat files\n'...
'../Release/field # path to predictive field matfile\n'];
fprintf(fileID, txt);
fclose(fileID);
end