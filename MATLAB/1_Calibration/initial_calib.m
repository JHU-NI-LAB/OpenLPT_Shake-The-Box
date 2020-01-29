function [camParaCalib, varargout] = initial_calib(Ximg, Xp3D, camParaknown, imgsize)
% This function calibrates camera parameters using Tsai's method: i.e., it
% assumes only radial distortion and ignores tangential distortion
% completely.
%
% syntax:
%	camParaCalib = calib_Tsai(Ximg, Xp3D, camParaknown, ...);
%	[camParaCalib, err_img] = calib_Tsai(Ximg, Xp3D, camParaknown, ...);
%
% inputs:
%   Ximg        --  particle coordinates on image plane (in pixels)
%   Xp3D        --  particle coordinates in 3D world system (in mm)
%   camParaknown    --  known camera parameters, including:
%                       Npixh   -- # of pixels along height
%                       Npixw   -- # of pixels along width
%                       hpix    -- pixel size along height (mm)
%                       wpix    -- pixel size along width (mm)
%   imgsize	--  experiement image size, [Npixx, Npixy], optional parameter, 
%		    used if the calibration image size is different from 
%		    experiment image size, default is the same.

%
% outputs:
%   camParaCalib    --  calibrated camera parameters, including 
%	camParaCalib.R = [r1 r2 r3]'
%	camParaCalib.T = [Tx Ty Tz]'
%	camParaCalib.f_eff
%	camParaCalib.hpix = camParaknown.hpix
%	camParaCalib.Npixh = camParaknown.Npixh
%	camParaCalib.Npixw = camParaknown.Npixw
%	camParaCalib.wpix
%	camParaCalib.k1
%	camParaCalib.k1star
%	% In Tsai's model, tangential distortion is ignored. Here p1 and p2 are
%	% kept for compatability with other models.
%	camParaCalib.p1 = 0
%	camParaCalib.p1star = 0
%	camParaCalib.p2 = 0
%	camParaCalib.p2star = 0
%	camParaCalib.err_x
%	camParaCalib.err_y
%	camParaCalib.err_t
%   err_img	--  error in particle center on image plane [err_x err_y]
%

if nargin == 3
	imgsize = [camParaknown.Npixw camParaknown.Npixh];
end



% in Tsai's method, the offset of camera axis to CMOS sensor is assumed to be 0.
camParaCalib.Noffh = 0;
camParaCalib.Noffw = 0;

Ximg(:,1) = Ximg(:,1) - (camParaknown.Npixw - imgsize(1))/2;
Ximg(:,2) = Ximg(:,2) - (camParaknown.Npixh - imgsize(2))/2;

ind = (Ximg(:,1) >= 0) & (Ximg(:,1) <= imgsize(1)) & (Ximg(:,2) >= 0) & (Ximg(:,2) <= imgsize(2));
Ximg = Ximg(ind, :);
Xp3D = Xp3D(ind, :);

Np = size(Ximg,1);
X = zeros(Np,1);
Y = zeros(Np,1);


camParaCalib.Noffw = 0;
camParaCalib.Noffh = 0;
X = Ximg(:,1) - camParaCalib.Noffw - imgsize(1)/2;
Y = imgsize(2)/2 - camParaCalib.Noffh - Ximg(:,2);  % image y-axis is pointing downward

% check if calibration points are co-planar or not
if sum(abs(Xp3D(:,1))) == 0
    % co-planar case: (assume the calibration mask lies on plane xw = 0)
    % the pixel size is assumed known
    X = X*camParaknown.wpix;
    Y = Y*camParaknown.hpix;
    % calculate radial distance from the image center, in pixel units, to be used
    % later for compute distortion coefficient
    rdsq = X.^2 + Y.^2;     
    % solve the radial alignment constraint by SVD for 5 parameters
    A = [Y.*Xp3D(:,2) Y.*Xp3D(:,3) Y -X.*Xp3D(:,2) -X.*Xp3D(:,3)];
    [U w V] = svd(A,0);
    x = U'*X;
    minsv = 1.E-10;
    a = zeros(5,1);
    for i=1:5
        if w(i,i) > minsv
            a = a + x(i)/w(i,i) * V(:,i);
        else
            warning('calib_Tsai: singular value of A1 < minsv, probably A1 is singular');
        end
    end

    Sr = a(1)^2 + a(2)^2 + a(4)^2 + a(5)^2;
    Ty = sqrt(2./(Sr + sqrt(Sr^2 - 4*(a(1)*a(5) - a(2)*a(4))^2))); % this is the magnitude of Ty

    % determine the sign of Ty
    % pick up a point that is far away from the center (to avoid offset of axis)
    i = 1;
    r2 = X(i)*X(i) + Y(i)*Y(i);
    rcr = (imgsize(1)/2/2*camParaknown.wpix)^2 + (imgsize(2)/2/2*camParaknown.hpix)^2;
    while r2 < rcr && i < numel(X);                
        i = i+1;
        r2 = X(i)*X(i) + Y(i)*Y(i);
    end
    xi = Xp3D(i,2:3)*a(1:2) + a(3);
    yi = Xp3D(i,2:3)*a(4:5) + 1;
    if xi*X(i)*yi*Y(i) > 0  % if the signs of xi and yi are consistent with X(i) and Y(i)
        if xi*X(i) < 0      % if xi and X(i) have different sign, then Ty is negative
            Ty = -Ty;
        end                 % if xi and X(i) have the same sign, then Ty is positive, do nothing
    else    % if the signs of xi, yi and X(i), Y(i) are not consistent, display a warning
        warning('calib_Tsai: inconsistent signs of xi, yi and X(i), Y(i)');
    end

    % now determine the rotation matrix R and Tx and Ty
    Tx = a(3)*Ty;
    r1 = zeros(3,1);
    r2 = zeros(3,1);
    r3 = zeros(3,1);
    r1(2:3) = a(1:2)*Ty;
    r2(2:3) = a(4:5)*Ty;
    r3(1) = r1(2)*r2(3) - r1(3)*r2(2);
    % The following are only known of their magnitude
    r1(1) = sqrt(1-r1(2)^2-r1(3)^2);
    r2(1) = sqrt(1-r2(2)^2-r2(3)^2);
    if (r1(2)*r2(2) + r1(3)*r2(3)) > 0
        r2(1) = -r2(1);
    end
    r3(2) = r1(3)*r2(1) - r1(1)*r2(3);
    r3(3) = r1(1)*r2(2) - r2(1)*r1(2);

    % determine the sign of the elements in R
    Xc = Xp3D * [r1 r2 r3];
    A = [Xc(:,1)+Tx  -X];
    [U w V] = svd(A,0);
    x = U' * (Xc(:,3).*X);
    a = zeros(2,1);
    for i = 1:2
        if w(i,i) > minsv
            a = a + x(i)/w(i,i) *V(:,i);
        else
            warning('calib_Tsai: singular value of A2 < minsv, probably A2 is singular');
        end
    end
    f_eff = a(1);
    Tz = a(2);
    if f_eff < 0
        r1(1) = -r1(1);
        r2(1) = -r2(1);
        r3(2) = -r3(2);
        r3(3) = -r3(3);
    end
    Xc = Xp3D*[r1 r2 r3];
    Xc(:,1) = Xc(:,1) + Tx;
    Xc(:,2) = Xc(:,2) + Ty;

    % now compute f_eff, Tz, and k by iteration
    max_iter = 50;
    iter = 0;
    k1 = 0;
    f_err = 1.;
    Tz_err = 1.;
    f_tol = 1.E-6;
    Tz_tol = 1.E-6;
    while ((iter< max_iter) && ((f_err > f_tol) || (Tz_err > Tz_tol)) )
        A = [Xc(:,1) -X.*(1+k1*rdsq)];
        b = X.*(1+k1*rdsq).*Xc(:,3);
        [U w V] = svd(A,0);
        x = U' * b;
        a = zeros(2,1);
        for i = 1:2
            if w(i,i) > minsv
                a = a + x(i)/w(i,i) *V(:,i);
            else
                txt = strcat('calib_Tsai: iter = ', num2str(iter), ': singular value of A3 < minsv, probably A2 is singular');
                warning(txt);
            end
        end
        f_err = abs(f_eff - a(1))/f_eff;
        Tz_err = abs((Tz - a(2))/Tz);
        f_eff = a(1);
        Tz = a(2);
        k1 = sum((Xc(:,1)*f_eff - X.*(Xc(:,3)+Tz)).*(X.*(Xc(:,3)+Tz).*rdsq))/sum((X.*(Xc(:,3)+Tz).*rdsq).^2);
        iter = iter+1;
    end

    dx_eff = camParaknown.wpix;
    k1star = k1*dx_eff^2*(imgsize(1)/2)^2 * (1+(imgsize(2) * camParaknown.hpix/(imgsize(1)*dx_eff))^2);

    objval = Tz_err;
else
    % if not co-planar, then use the most general method

    % solve the radial alignment constraint by SVD
    A = [Y.*Xp3D(:,1) Y.*Xp3D(:,2) Y.*Xp3D(:,3) Y -X.*Xp3D(:,1) -X.*Xp3D(:,2) -X.*Xp3D(:,3)];
    [U w V] = svd(A,0);
    x = U' * X;
    minsv = 1.E-10;     % minimum singular value, below this matrix A is considered as singular
    a = zeros(7,1);
    for i=1:7
        if w(i,i) > minsv
            a = a + x(i)/w(i,i) * V(:,i);
        else
            warning('calib_Tsai: singular value of A < minsv, probably A is singular');
        end
    end

    Ty = 1/sqrt(sum(a(5:7).*a(5:7)));   % magnitude of Ty
    fxfy = sqrt(sum(a(1:3).*a(1:3))) * Ty;  % fx/fy = sx*dy/dx

    % Now determine the sign of Ty
    % pick up a point that is far away from the center (to avoid offset of axis)
    i = 1;
    r2 = X(i)*X(i) + Y(i)*Y(i);
    rcr = (imgsize(1)/2/2)^2 + (imgsize(2)/2/2)^2;
    while r2 < rcr
        i = i+1;
        r2 = X(i)*X(i) + Y(i)*Y(i);
    end
    xi = Xp3D(i,:)*a(1:3) + a(4);
    yi = Xp3D(i,:)*a(5:7) + 1;
    if xi*X(i)*yi*Y(i) > 0  % if the signs of xi and yi are consistent with X(i) and Y(i)
        if xi*X(i) < 0      % if xi and X(i) have different sign, then Ty is negative
            Ty = -Ty;
        end                 % if xi and X(i) have the same sign, then Ty is positive, do nothing
    else    % if the signs of xi, yi and X(i), Y(i) are not consistent, display a warning
        warning('calib_Tsai: inconsistent signs of xi, yi and X(i), Y(i)');
    end

    % solve for the rotation matrix R and the translation along xc-axis
    r1 = a(1:3)*Ty/fxfy;
    r2 = a(5:7)*Ty;
    % r3 is defined such that the camera coordinate system is left-hand
%     r3 = [-r1(2)*r2(3)+r1(3)*r2(2);
%           -r1(3)*r2(1)+r1(1)*r2(3);
%           -r1(1)*r2(2)+r1(2)*r2(1)];
    r3=cross(r1,r2);
    Tx = a(4)*Ty/fxfy;

    % particle coordinates in camera system, since Tz is not known yet, we set Tz = 0
    Xc = Xp3D * [r1 r2 r3];
    Xc(:,1) = Xc(:,1) + Tx;
    Xc(:,2) = Xc(:,2) + Ty;
    rc = sqrt(Xc(:,1).*Xc(:,1) + Xc(:,2).*Xc(:,2));

    % solve for fx, Tz and k1 by iteration
    rs2 = X.*X + (fxfy*fxfy)*(Y.*Y);
    rs = sqrt(rs2);
    rs3 = rs2.*rs;
    sumrs6 = sum(rs3.*rs3);
    k1s = 0;    % initial guess of k1

%             
%             A = [Xc(:,1) -X];
%             A1 = [Xc(:,2) -Y]
%             B = Xc(:,3).*X;
%             B1 = Xc(:,3).*Y;
%             A = [A;A1];
%             B = [B;B1];
%             
% 
%             bb=A'*A;
% 
%             bc=A'*B;
% 
%             result=bb\bc;
%             
%             dummy = result(1)./(Xc(:,3)+result(2));
%             X_approx = dummy.*Xc(:,1);
%             Y_approx = dummy.*Xc(:,2);
%             err_t = sqrt(sum((X_approx-X).^2 + (Y_approx-Y).^2)/Np);

    Maxiter = 1;
    for niter = 1:Maxiter
        d = rs - k1s*rs3;
        A = [rc -d];
        [U w V] = svd(A,0);
        x = U' * (Xc(:,3).*d);
        a = x(1)/w(1,1)*V(:,1) + x(2)/w(2,2)*V(:,2);
        fx = a(1);
        Tz = a(2);
        dummy = fx./(Xc(:,3)+Tz);
        k1s = sum(rs3.*(rs - dummy.*rc))/sumrs6;
        % calculate residual errors in image coordinates
        dummy = dummy./(1-k1s*rs2);
        err_r = sqrt(sum((dummy.*rc-rs).^2)/Np);    % radial error
        X_approx = dummy.*Xc(:,1);
        Y_approx = dummy.*Xc(:,2);
        err_t = sqrt(sum((X_approx-X).^2 + (Y_approx-Y).^2)/Np);    % total distance error
    end

    objval = err_t;

    % with known camera parameters, rescale k1, f, and dx/Sx
    dx_eff = camParaknown.hpix/fxfy;    % effective pixel size in x-direction on CMOS sensor (mm)
    f_eff = fx*dx_eff;      % effective focal length (mm)
    k1 = k1s/(dx_eff^2);    % distortion coefficient, dimensional (mm^-2)
    % distortion coefficient, nondimensional
    k1star = k1s*(imgsize(1)/2)^2 * (1+(imgsize(2) * camParaknown.hpix/(imgsize(1)*dx_eff))^2);
end


% camParaCalib.Noffw = offset(1);
% camParaCalib.Noffh = offset(2);

% find the inverse matrix for transferring back from camera coords to world coords
T = [[[r1 r2 r3]' [Tx; Ty; Tz]]; [0 0 0 1]];
Tinv = T^(-1);
camParaCalib.R = [r1 r2 r3]';
camParaCalib.T = [Tx Ty Tz]';
camParaCalib.Rinv = Tinv(1:3, 1:3);
camParaCalib.Tinv = Tinv(1:3, 4);
camParaCalib.f_eff = f_eff;
camParaCalib.hpix = camParaknown.hpix;
camParaCalib.wpix = dx_eff;
camParaCalib.Npixh = imgsize(2);
camParaCalib.Npixw = imgsize(1);
camParaCalib.k1 = k1;
camParaCalib.k1star = k1star;
% In Tsai's model, tangential distortion is ignored. Here p1 and p2 are
% kept for compatability with other models.
camParaCalib.p1 = 0;
camParaCalib.p1star = 0;
camParaCalib.p2 = 0;
camParaCalib.p2star = 0;
Xproj = calibProj_Tsai(camParaCalib, Xp3D);
% camParaCalib.err_t = sqrt(sum(sum((Ximg - calibProj_Tsai(camParaCalib, Xp3D)).^2))/Np);
err_img = Xproj - Ximg;
camParaCalib.err_x = sqrt(mean((err_img(:,1)).^2));
camParaCalib.err_y = sqrt(mean((err_img(:,2)).^2));
camParaCalib.err_t = sqrt((camParaCalib.err_x)^2 + (camParaCalib.err_y)^2);
[m, I] = max(sqrt((Ximg(:,1)-Xproj(:,1)).^2 + (Ximg(:,2)-Xproj(:,2)).^2));
[m I Ximg(I,:) Xp3D(I,:)];
if nargout >= 2
	varargout(1) = {err_img};
	if nargout == 3
		varargout(2) = {[Ximg Xp3D]};
	end
end

end