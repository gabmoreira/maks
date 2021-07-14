function [R1,t1,R2,t2] = computePose(E)
% COMPUTEPOSE Compute relative pose from Essential matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

[U, S, V] = svd(E);

assert(S(3) < 1e-10, 'Essential matrix must be rank 2');

S = [1,0,0;0,1,0;0,0,0];
E = U * S * V';

R1 = U * [0, 1, 0; -1, 0, 0; 0, 0, 1] * V';
R2 = U * [0, 1, 0; -1, 0, 0; 0, 0, 1]' * V';
t1 = U(:,3);
t2 = -U(:,3);

res1 = norm(crossMat(t1)*R1 / norm(crossMat(t1)*R1, 'fro') - (E / norm(E, 'fro')), 'fro');
res2 = norm(crossMat(t2)*R2 / norm(crossMat(t2)*R2, 'fro') - (E / norm(E, 'fro')), 'fro');

assert(res1 < 1e-10, 'Error with 1');
assert(res2 < 1e-10, 'Error with 2');

end

