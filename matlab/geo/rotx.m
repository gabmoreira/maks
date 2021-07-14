function R = rotx(theta)
% ROTX Rotation matrix around x-axis.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

R = [1,0,0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)];

end

