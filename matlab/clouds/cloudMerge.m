function dst = cloudMerge(srcCloud1, srcCloud2)

% CLOUDMERGE Sets all points to same color
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

dst.pts =   [srcCloud1.pts; srcCloud2.pts];
dst.color = [srcCloud1.color; srcCloud2.color];

end

