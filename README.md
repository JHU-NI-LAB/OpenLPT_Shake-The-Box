# ShakeTheBox

ShakeTheBox provides a C++ code for processing images in order to obtain tracks of particles seeded in the flow. Usually at least three cameras are needed to reconstruct 3D tracks. Anyone who uses this code and develop this code for their research should cite the following publication in their work：

Tan, S., Salibindla, A., Masuk, A.U.M. and Ni, R., 2019. An open-source shake-the-box method and its performance evaluation. In 13th International Symposium on Particle Image Velocimetry.

We really welcome any research group that is willing to improve the code with us. Please send us email to join as contributor to the code. 

Instructions on how to use this code are listed as follows:

(1) Preprocessing images, to make particles more bright and images less noise. A sample code of processing images can be found in  ./Data_analysis_process/PreprocessImage.m. Since different camera gives different images, it is not suggested to use this sample code directly, but to use it as a reference.

(2) Setting up a project folder, and put the images under that folder. It is suggested to set up folder like cam1, cam2, cam3 and so on under the project folder, and place the corresponding processed images beneath each camera folder.

(3) Creat a camXImageNames.txt under the project folder to list the relative path to the images. For cam1, cam1ImageNames.txt should be like:

cam1/cam1frame00001.tif

cam1/cam1frame00002.tif

cam1/cam1frame00003.tif

cam1/cam1frame00004.tif

...

(4) Generate the calibration file. It is highly recommended to use ./Calibration/calibGUI.m to get the calibration file. After getting the file, use ./Data_analysis_process/CalibMatToTXTV2.m to change the calibration mat file to txt file.

(5) Generate the configuration file. You can use the matlab file: ./Data_analysis_process/GenerateConfigFile.m to generate the configuration file. There are four arguments, the path of the project folder, start frame NO, end frame NO., and calibration file name. After generation of configuration files, you can open each of the file and modify the parameters to accomodate your experiments.

(6) Compile the C++ code. Use Ecllipse to set up the ShakeTheBox project with the folder ./inc and ./src. Compile the code with g++ compiler, remember to add lable -fopenmp to enable parallelization. The project setting can be found in ./cproject. It is recommended to use ./cproject directly.

(7) Run ShakeTheBox. Go to the compile folder, and run the executable file ShakeTheBox. The argument should be the path of your project folder. After Run the code, you need to key in 0 twice to start the process. There are also other options for development use. 

(8) After data processing, all the tracks are saved in the folder Tracks/ConvergedTracks/ under the project folder. You can use ./Data_analysis_process/ReadAllTracks.m to read all the track data to the workspace in MATLAB. The tracks are save in a format as: X, Y, Z, frame NO, track NO.

A sample dataset used to test our code can be found in the following link:
https://drive.google.com/drive/folders/1Sx1rkMdyt-0LKw-s-LlI__96B4EJ1o5U?usp=sharing

If you have any question on how to use the code or you are interested in improving the code, please contact me: rui.ni@jhu.edu
