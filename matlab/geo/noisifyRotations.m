function Rtilde = noisifyRotations(R, A, sigma)
% NOISIFYROTATIONS Takes absolute rotations and computes noisy pairwise
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 28-February-2021

[ii,jj] = find(triu(A));
m = numel(ii);
n = size(A,1);

% Generate noise rotations
noise = rotrnd(m, sigma);

Rtilde = speye(3*n, 3*n) + kron(A, ones(3,3)) .* (R * R');

% Add noise
for k=1:m
    i = ii(k);
    j = jj(k);
    Rtilde(i*3-2:i*3,j*3-2:j*3) = noise(k*3-2:k*3,:) * Rtilde(i*3-2:i*3,j*3-2:j*3);
    Rtilde(j*3-2:j*3,i*3-2:i*3) = Rtilde(i*3-2:i*3,j*3-2:j*3)';
end

end

