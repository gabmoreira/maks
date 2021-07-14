function q = isSO3(mat)
% ISSO3 Checks if a matrix is in SO(3)
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

if (size(mat,1) ~= 3) || (size(mat,2) ~= 3)
    q = false;
elseif abs(det(mat)-1) > eps
    q = false;
elseif norm(mat * mat' - eye(3), 'fro') > eps
    q = false;
else
    q = true;
end

end

