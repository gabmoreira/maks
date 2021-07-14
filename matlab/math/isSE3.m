function v = isSE3(mat)
% ISSE3 Checks if matrix is in SE(3)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

eps = 1e-6;

if size(mat,1) ~= 4 || size(mat,2) ~= 4
    v = false;
    
elseif ~isSO3(mat(1:3,1:3))
    v = false;
    
elseif (abs(mat(4,1)) > eps) || (abs(mat(4,2)) > eps) || (abs(mat(4,3)) > eps) || (abs(mat(4,4)-1) > eps)
    v = false;
    
else
    v = true;
end

end

