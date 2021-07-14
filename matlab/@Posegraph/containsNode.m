function idx = containsNode(obj, id)
% CONTAINSNODE Check if specified node is in the posegraph.
%
% Other m-files required: Posegraph.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

q = find(id == [obj.nodes.id{:}], 1);

if isempty(q)
    idx = obj.numNodes + 1;
else
    idx = q;
end

end

