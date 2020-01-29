function rotm = eul2rotm(eul, sequence)
    if (size(eul,1) ~= 3)
        error('eul2rotm: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end

    if ~exist('sequence', 'var')
        % use the default axis sequence ...
        sequence = 'ZYX';
    end
    rotm = zeros(3,3);

    s_1 = sin(eul(1,1)); % theta_z or theta_z1
    c_1 = cos(eul(1,1));
    s_2 = sin(eul(2,1)); % theta_y
    c_2 = cos(eul(2,1));
    s_3 = sin(eul(3,1)); % theta_x or theta_z2
    c_3 = cos(eul(3,1));

    %% Convert the given Euler angles theta for the x, y and z-axis into the corresponding
    %  direction cosine rotation matrix R, in dependency of the axis rotation sequence for
    %  the multiplication order of the rotation factors:
    % For further details see:
    %   [1] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/EulerAngles.pdf>, p. 9 & 16.
    %   [2] MATLAB TOOLBOX FOR RIGID BODY KINEMATICS, Hanspeter Schaub & John L. Junkins,
    %       9th AAS/AIAA Astrodynamics Specialist Conference, AAS 99-139, 1999, <http://hanspeterschaub.info/Papers/mtb1.pdf>, p. 4.
    %   [3] GitHub: ShoolCode/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox,
    %       <https://github.com/ZachDischner/SchoolCode/tree/master/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox/>
    %   [4] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
    %       pp. 31-32, formulas (2.18) and (2.20).
    switch sequence
        case 'ZYX'
            %            |c_1*c_2    c_1*s_2*s_3 - s_1*c_3    c_1*s_2*c_3 + s_1*s_3|
            % R(Theta) = |s_1*c_2    s_1*s_2*s_3 + c_1*c_3    s_1*s_2*c_3 - c_1*s_3|
            %            |   -s_2                  c_2*s_3                  c_2*c_3|
            rotm(1,1) =  c_1*c_2;
            rotm(1,2) =  c_1*s_2*s_3 - s_1*c_3;
            rotm(1,3) =  c_1*s_2*c_3 + s_1*s_3;

            rotm(2,1) =  s_1*c_2;
            rotm(2,2) =  s_1*s_2*s_3 + c_1*c_3;
            rotm(2,3) =  s_1*s_2*c_3 - c_1*s_3;

            rotm(3,1) = -s_2;
            rotm(3,2) =  c_2*s_3;
            rotm(3,3) =  c_2*c_3;
        case 'ZYZ'
            %            |c_1*c_2*c_3 - s_1*s_3   -c_1*c_2*s_3 - s_1*c_3    c_1*s_2|
            % R(Theta) = |s_1*c_2*c_3 + c_1*s_3   -s_1*c_2*s_3 + c_1*c_3    s_1*s_2|
            %            |             -s_2*c_3                  s_2*s_3        c_2|
            rotm(1,1) =  c_1*c_2*c_3 - s_1*s_3;
            rotm(1,2) = -c_1*c_2*s_3 - s_1*c_3;
            rotm(1,3) =  c_1*s_2;

            rotm(2,1) =  s_1*c_2*c_3 + c_1*s_3;
            rotm(2,2) = -s_1*c_2*s_3 + c_1*c_3;
            rotm(2,3) =  s_1*s_2;

            rotm(3,1) = -s_2*c_3;
            rotm(3,2) =  s_2*s_3;
            rotm(3,3) =  c_2;
        otherwise
            error('eul2rotm: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end