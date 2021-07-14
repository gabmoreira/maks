function addEdge(obj, id, varargin)
% ADDEDGE Adds edge to the posegraph.
%
% Other m-files required: Posegraph.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

p = inputParser;

validId = @(x) all(numel(id) == 2);

% Node id from / node id to
addRequired(p, 'id', validId);

% Any other parameter must have been specified in configEdges
defaultParam = [];
validParam   = @(x) true;

for field=obj.nodeFields(2:end)
    addParameter(p, field, defaultParam, validParam);
end

parse(p,id,varargin{:});

% Check if graph already contains this edge
idx = obj.containsEdge(p.Results.id);

% Add the new information.
for field=obj.edgeFields(1:end)
    obj.edges.(field){idx} = p.Results.(field);
end

if (idx > obj.numEdges)
    obj.numEdges = obj.numEdges + 1;
end

end

