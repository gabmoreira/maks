function A = adjacency(obj, varargin)
% ADJACENCY Builds sparse adjacency matrix from specified field
%
% Other m-files required: Posegraph.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

% Parse input arguments
validField   = @(x) isstring(x);
defaultField = 'default';

p = inputParser;
addParameter(p, 'field', defaultField, validField);
parse(p, varargin{:});

% Sort node ids and create respective indices
nodes = sort([obj.nodes.id{:}]); 
nodesIdx(nodes) = 1:obj.numNodes;

edgeIds = [obj.edges.id{:}];

% Edge indices
ij = [nodesIdx(edgeIds(1,:)); nodesIdx(edgeIds(2,:))];

switch(p.Results.field)
    
    % Default case returns binary adjacency matrix
    case 'default'
        A = blockMatrix(ones(1,obj.numEdges), 1, ij, 'symmetric', true);
        
    % Otherwise returns block matrix
    otherwise
        edgeData = obj.edges.(p.Results.field);
        edgeData = [edgeData{:}];
        stride   = size(obj.edges.(p.Results.field){1}, 2);
        A = blockMatrix(edgeData, stride, ij, 'symmetric', true);
end

end

