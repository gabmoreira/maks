function R = quatToRot(q)
% QUATTOROT Unit quaternion to rotation matrix in SO(3).
%
% q = [qx, qy, qz, qw], with qw the real part.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

q = q / sqrt(q' * q);

qr = q(4);
qi = q(1);
qj = q(2);
qk = q(3);

% First row of the rotation matrix
r00 = 1 - 2 * (qj^2 + qk^2);
r01 = 2 * (qi * qj - qk * qr);
r02 = 2 * (qi * qk + qj * qr);

% Second row of the rotation matrix
r10 = 2 * (qi * qj + qk * qr);
r11 = 1 - 2 * (qi^2 + qk^2);
r12 = 2 * (qj * qk - qi * qr);

% Third row of the rotation matrix
r20 = 2 * (qi * qk - qj * qr);
r21 = 2 * (qj * qk + qi * qr);
r22 = 1 - 2 * (qi^2 + qj^2);

R = [r00, r01, r02; r10, r11, r12; r20, r21, r22];
end

