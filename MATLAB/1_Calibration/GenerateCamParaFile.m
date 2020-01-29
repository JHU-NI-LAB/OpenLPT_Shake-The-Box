function GenerateCamParaFile(filepath, txtname)
for i = 1 : 4
    campara = load([filepath 'cam' num2str(i) '.mat']);
    camParaCalib(i) = campara.camParaCalib;
end
CameraDistribution(camParaCalib)
save([filepath txtname '.mat'], 'camParaCalib');
CalibMatToTXT(camParaCalib, [filepath txtname '.txt']);
end

