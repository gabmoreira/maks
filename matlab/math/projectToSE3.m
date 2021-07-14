function dst = projectToSE3(src)
% PROJECTTOSE3 Computes the matrix in SE(3) which is the closest to mat.
%
% Syntax: dst = projectToSE3(src)
%
% Inputs:
%    src - 4 x 4 matrix
%
% Outputs:
%    dst - Rigid transformation matrix closest to src (min ||dst - src||_F)
%    
% Other m-files required: projectToSO3
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

dst          = src ./ src(4,4);
dst(1:3,1:3) = projectToSO3(dst(1:3,1:3));
dst(4,:)     = [0, 0, 0, 1];

end

