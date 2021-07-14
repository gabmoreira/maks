function dst = projectToSO3(src)
% PROJECTTOSO3 Computes the matrix in SO(3) which is the closest to mat.
%
% Syntax: dst = projectToSO3(src)
%
% Inputs:
%    src - 3 x 3 matrix
%
% Outputs:
%    dst - 3D rotation matrix closest to mat (min ||dst - src||_F)
%    
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

[U, ~, V]     = svd(real(src));
detCorrection = diag([1,1,det(U*V')]);
dst           = U * detCorrection * V';

end
