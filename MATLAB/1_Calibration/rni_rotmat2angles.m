function [angles rotout] = rni_rotmat2angles(rot)
%This function calculates the angles that correspond to a rotation matrix.
%It is needed so that dynamic calibration can be performed on the minimum
%number of parameters rather than on the full rotation matrix.  

%inputs: 
%       rot:  a rotation matrix 
%outputs:
%       angles: three angles that define the rotation matrix--theta, phi,
%       and mu where (starting at Tinv) mu is first rotation about z axis, phi is second
%       rotation around y axis and theta is third rotation around z axis.
%
%       rotout:  best match rotation matrix
%
%some checks 
%determinant = det(rot) %should be -1 (for some reason this calibration uses reflection and rotation matrix
%check1=(rot(1,3)^2 + rot(2,3)^2+rot(3,3)^2) % should be 1
%check2 = (rot(3,1)^2 + rot(3,2)^2+rot(3,3)^2) %should be 1

% Rui changed it to exclude the ambiguity of angles
% new algorithm based on
% https://truesculpt.googlecode.com/hg-history/38000e9dfece971460473d5788c235fbbe82f31b/Doc/rotation_matrix_to_euler.pdf


phi = -1*acos(-1*rot(3,3));  % As I first did this, phi is always between -pi and 0--probably should be changed

mu= atan2(-1*rot(3,2)/sin(phi),rot(3,1)/sin(phi)); 
theta = atan2(rot(2,3)/sin(phi),rot(1,3)/sin(phi));
%theta=mod(theta,pi)+pi;  %this just seems to work -DB 2/2010

if det(rot) > 0  
   mu = -1*mu ;    %These corrections change the handedness of the coord system.
   theta=mod(theta,pi)+pi;
end % if initial coordinate system is left handed

rotztheta=[cos(theta), -sin(theta), 0; +sin(theta), cos(theta), 0; 0 0 1];
rotzmu=[cos(mu), -sin(mu), 0; +sin(mu), cos(mu), 0; 0 0 1];
rotyphi=[cos(phi), 0, sin(phi); 0 1 0; -sin(phi), 0, cos(phi)];
reflect = [1 0 0;0 1 0;0 0 -1];

rotout = reflect*rotztheta*rotyphi*rotzmu;

if ~isempty(nonzeros(rotout==rot))
    phi=-phi;
    mu= atan2(-1*rot(3,2)/sin(phi),rot(3,1)/sin(phi)); 
    theta = atan2(rot(2,3)/sin(phi),rot(1,3)/sin(phi));
end
%some checks on the output
% diff = rot - rotout
% fracdiff = diff ./ rot

angles=[theta, phi, mu];