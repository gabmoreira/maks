function addNode(obj, id, varargin)
% ADDNODE Adds node to the posegraph.
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

% Field id is always required to identify the node
validId      = @(x) isscalar(x);
addRequired(p, 'id', validId);

% Any other parameter must have been specified in configNodes
defaultParam = [];
validParam   = @(x) true;

for field=obj.nodeFields(2:end)
    addParameter(p, field, defaultParam, validParam);
end

parse(p, id, varargin{:});

idx = obj.containsNode(id);

% Create node struct
for field=obj.nodeFields(1:end)
    obj.nodes.(field){idx} = p.Results.(field);
end

% Update no. nodes
if (idx > obj.numNodes)
    obj.numNodes = obj.numNodes + 1;
end
end

