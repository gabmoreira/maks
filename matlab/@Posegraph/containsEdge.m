function idx = containsEdge(obj, id)
% CONTAINSEDGE Checks if posegraph has edge
%
% Other m-files required: Posegraph.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

assert(isnumeric(id) && (numel(id) == 2) && (floor(id(1)) == id(1)) ...
       && (floor(id(2)) == id(2)), 'Invalid edge');
   
if obj.numEdges > 0
    ids = obj.edges.id;
    [q, idx] = ismember(reshape(id,1,2), [ids{:}]', 'rows');
    
    if (~q)
        idx = obj.numEdges + 1;
    end
else
    idx = 1;
end

end


