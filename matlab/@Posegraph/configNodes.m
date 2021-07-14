function configNodes(obj, varargin)
% CONFIGNODES Sets attributes to graph nodes
%
% Other m-files required: Posegraph.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

for i=1:nargin-1
    assert(isstring(varargin{i}), "Node fields must be strings");
end

for i=1:nargin-1
    obj.nodeFields = [obj.nodeFields, varargin{i}];
end

fprintf("Graph nodes have the following fields:\n");
for i=1:numel(obj.nodeFields)
    obj.nodes.(obj.nodeFields{i}) = {};
    fprintf(strcat(" o ", obj.nodeFields{i}, "\n"));
end

end

