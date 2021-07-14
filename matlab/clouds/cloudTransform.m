function dst = cloudTransform(src, R, t)
% CLOUDTRANSFORM Applies rigid transformation to point cloud
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

dst.pts = dst.pts(:,1:3) * R' + t';

end

