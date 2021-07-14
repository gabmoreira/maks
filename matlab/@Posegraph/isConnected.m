function [q, num_components] = isConnected(obj)
% ISCONNECTED Checks if graph is connected and returns number of connected
% components. This is done by checking the spectrum of the graph laplacian.
%
% Other m-files required: Posegraph.m, laplacian.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         g.antunes.moreira (at) gmail.com
% Website:       https://www.github.com/gabmoreira/computervision
% Last revision: 12-June-2020

%------------------------------- BEGIN CODE -------------------------------
eps = 1e-10;
L = obj.laplacian;
lambdas = eig(L);
num_components = sum(abs(lambdas) < eps);

% Number of connected components of a graph is given by the number of
% null eigenvalues (counting repeated eigenvalues multiple times)
q = (num_components == 1);
%------------------------------ END OF CODE -------------------------------
end

