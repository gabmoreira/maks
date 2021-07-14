function dst = rankReduction(src, rank)
% RANKREDUCTION Projects src matrix onto lower rank
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

[U,S,V] = svd(src);
dst = U(:,1:rank) * S(1:rank,1:rank) * V(:,1:rank)';

end

