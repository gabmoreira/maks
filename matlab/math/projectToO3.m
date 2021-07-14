function dst = projectToO3(src)
% PROJECTTOO3 Computes the matrix in O(3) which is the closest to mat.
%
% Syntax: dst = projectToSO3(src)
%
% Inputs:
%    src - 3 x 3 matrix
%
% Outputs:
%    dst - 3D rotation matrix closest to src (min ||dst - src||_F)
%    
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 25-March-2020

[U, ~, V] = svd(real(src));
dst       = U * V';

end
