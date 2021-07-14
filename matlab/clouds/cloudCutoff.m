function dst = cloudCutoff(src, cutoff, dim)
% CLOUDCUTOFF Deletes points outside cutoff interval along specified dim
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

idx = (src.pts(:,dim) > cutoff(1)) & (src.pts(:,dim) < cutoff(2));

% Get fields in the src point clouds structure
fields = fieldnames(src);

for i=1:size(fields,1)
    aux = src.(fields{i});
    dst.(fields{i}) = aux(idx, :);
end

end
