function K = intrinsicsMatrix(fx, fy, cx, cy, s)
% INTRINSICSMATRIX Create camera intrinsics matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

K = [fx,  s, cx;
      0, fy, cy;
      0,  0,  1];

end

