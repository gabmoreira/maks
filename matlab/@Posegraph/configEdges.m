function configEdges(obj, varargin)
% CONFIGEDGES Sets attributes to graph edges
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
    assert(isstring(varargin{i}), "Edge fields must be strings");
end

for i=1:nargin-1
    obj.edgeFields = [obj.edgeFields, varargin{i}];
end

fprintf("Graph edges have the following fields:\n");
for i=1:numel(obj.edgeFields)
    obj.edges.(obj.edgeFields{i}) = {};
    fprintf(strcat(" o ", obj.edgeFields{i}, "\n"));
end

end

