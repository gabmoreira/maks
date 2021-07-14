function R = angleAxisToRot(u, theta)
% ANGLEAXISTOROT Angle-Axis to matrix representation.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

c = cos(theta);
s = sin(theta);
cc = 1 - c;
R = [    c + (u(1)^2)*cc,  u(1)*u(2)*cc-u(3)*s,  u(1)*u(3)*cc+u(2)*s;
     u(1)*u(2)*cc+s*u(3),        c+cc*(u(2)^2),  u(2)*u(3)*cc-u(1)*s;
     u(1)*u(3)*cc-u(2)*s,  u(3)*u(2)*cc+u(1)*s,       c+(u(3)^2)*cc];
end

