function nodeOp(obj, op, field)
% NODEOP Applies function to all nodes of the graph on a specified field
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

obj.nodes.(field) = cellfun(op, obj.nodes.(field), 'UniformOutput', false);
    
end

