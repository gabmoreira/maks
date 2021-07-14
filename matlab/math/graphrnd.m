function A = graphrnd(n, m, diagpull)
% GRAPHRND Creates adjacency matrix of a random graph
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

A = spdiags(ones(n,2), [1,-1], n, n);

k = 1;
p = 1;
T = zeros(m*2,3);

while (k <= m)
    i = randi(n);
    j = randi(n);
    if (i ~= j)
        if (abs(i-j) < (diagpull + randn*0.1) * n)
            T(p,:) = [i,j,1]; p = p + 1;
            T(p,:) = [j,i,1]; p = p + 1;
            k = k+1;
        end
    end
end

A = A + sparse(T(:,1), T(:,2), T(:,3), n, n);
A = A > 0;

end

