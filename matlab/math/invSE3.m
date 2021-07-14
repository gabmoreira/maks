function dst = invSE3(src)
% INVSE3 Inverse of a matrix belonging to the special Eucliden group.
%
% Syntax:  dst = invSE3(src)
%
% Inputs:
%    src - 4 x 4 matrix in SE(3)
%
% Outputs:
%    dst - Inverse of mat
%    
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:  Gabriel Moreira
% email:   gmoreira (at) isr.tecnico.ulisboa.pt
% Website: https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

dst          = eye(4,4);
dst(1:3,1:3) = src(1:3,1:3)';
dst(1:3,4)   = -dst(1:3,1:3) * src(1:3,4);

end