function  L = laplacian(obj, varargin)
% LAPLACIAN Builds sparse graph Laplacian matrix
%
% Other m-files required: Posegraph.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

validField   = @(x) isstring(x);
defaultField = 'default';

p = inputParser;
addParameter(p, 'field', defaultField, validField);
parse(p, varargin{:});

A = obj.adjacency('field', string(p.Results.field));
n = size(A,1) / obj.numNodes;

L = kron(obj.degree('asmatrix', true), speye(n));
L = L - A;

end

