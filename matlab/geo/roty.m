function R = roty(theta)
% ROTY Rotation matrix around y-axis.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];

end

