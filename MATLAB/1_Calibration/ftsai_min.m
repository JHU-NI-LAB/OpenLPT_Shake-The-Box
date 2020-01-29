% f = Tsai_8b(params, params_const, xw, yw, zw, X, Y) 
% 
% ********************************************************************************************** 
% *******         Calibrating a Camera Using a Monoview Coplanar Set of Points           ******* 
% ********************************************************************************************** 
%                              6/2004   Simon Wan  
%                              simonwan@hit.edu.cn 
% 
% Note:    This is not called directly but as a function handle from the "fminsearch " 
% 
function f = ftsai_min(params, params_const, Xtest3D, Ximg) 
% unpack the params 
    angles=params(1:3);
    
    if size(angles,1)~=3
        angles=angles';
    end
    rotm = eul2rotm(angles);
    

    R = rotm;
    
    T = params(4:6);
    f = params(7);
    k1 = params(8);
    offw = params(9);
    offh = params(10);
    
    wpix = params_const(1);
    hpix = params_const(2);
    Npixw = params_const(3);
    Npixh = params_const(4);
    
    
    Xc = Xtest3D * R';
    Xc(:,1) = Xc(:,1) + T(1);
    Xc(:,2) = Xc(:,2) + T(2);
    Xc(:,3) = Xc(:,3) + T(3);
    dummy = f./Xc(:,3);
    Xu = Xc(:,1).*dummy;  % undistorted image coordinates
    Yu = Xc(:,2).*dummy;
    ru2 = Xu.*Xu + Yu.*Yu;
    dummy = 1+k1*ru2;
    Xd = Xu./dummy;
    Yd = Yu./dummy;
    %iterate once
    %dummy = 1 + k1*(Xd.*Xd + Yd.*Yd);
     
    %Xd = Xu.*dummy;
    %Yd = Yu.*dummy;
    Np = size(Xtest3D,1);
    Xtest_proj = zeros(Np,2);
    Xtest_proj(:,1) = Xd/wpix + offw + Npixw/2;
    Xtest_proj(:,2) = Npixh/2 - offh - Yd/hpix;
    f = sqrt(mean(sum((Xtest_proj-Ximg).^2,2)));

