function dst = cloudSetColor(src, color)
% CLOUDSETCOLOR Sets all points to same the color
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

dst = src;
dst.color = ones(size(dst.pts,1), 3) .* color;

end

