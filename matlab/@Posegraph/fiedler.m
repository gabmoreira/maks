function f = fiedler(obj)
% FIEDLER Computes the graph's Fiedler value.
%
% Other m-files required: Posegraph.m, laplacian.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         g.antunes.moreira (at) gmail.com
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

L = obj.laplacian();
sigma = -1e-6;
lambdas = eigs(L, 2, sigma);
f = max(lambdas);

end

