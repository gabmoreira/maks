function R = rotz(theta)
% ROTZ Rotation matrix around z-axis.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];

end

