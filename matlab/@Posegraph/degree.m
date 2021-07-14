function D = degree(obj, varargin)
% DEGREE Builds degree matrix or vector (if 'asmatrix' is true/false)
%
% Other m-files required: Posegraph.m, adjacency.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

p = inputParser;

validAsmatrix = @(x) islogical(x);
defaultAsmatrix = false;

addParameter(p,'asmatrix', defaultAsmatrix, validAsmatrix);

parse(p,varargin{:});

A = obj.adjacency();
D = full(sum(A,2)); 
D = D(:);

if (p.Results.asmatrix)
    D = spdiags(reshape(D, obj.numNodes,1), 0, obj.numNodes, obj.numNodes);
end
end

