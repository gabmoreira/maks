function A = cycleGraph(n)
% CYCLEGRAPH Creates sparse adjacency matrix for a cycle graph
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

A = spdiags(ones(n,2), [-1, 1], n, n);
A(n,1) = 1;
A(1,n) = 1;
    
end

