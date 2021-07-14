function dst = cloudDownsample(src, factor)
% CLOUDDOWNSAMPLE Random point cloud downsampling
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

idx     = randperm(size(src.pts,1));
newSize = round(factor * size(src.pts,1));

% Get fields in the src point clouds structure
fields = fieldnames(src);

for i=1:size(fields,1)
    aux = src.(fields{i});
    dst.(fields{i}) = aux(idx(1:newSize), :);
end

end

