function dst = bilinearInterp(src, distance)
% BILINEARINTERP Performs bilinear interpolation of an image
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

ex = @(k,l) min(mod(k,2*l), 2*l-1-mod(k,2*l));

[n,m] = size(src);
nnew = floor(n / distance);
mnew = floor(m / distance);

dst = zeros(nnew, mnew);

for inew=0:nnew-1
    for jnew=0:mnew-1
        i = distance * inew;
        j = distance * jnew;
        ii = floor(i);
        ji = floor(j);
        
        dst(inew+1,jnew+1) = (i-ii) * (j-ji) * src(ex(ii+1,n)+1, ex(ji+1,m)+1) + ...
                            (1+ii-i) * (j-ji) * src(ex(ii,n)+1, ex(ji+1,m)+1) + ...
                            (i-ii) * (1+ji-j) * src(ex(ii+1,n)+1, ex(ji,m)+1) + ...
                            (1+ii-i) * (1+ji-j) * src(ex(ii,n)+1, ex(ji,m)+1); 
    end
end

dst = 255*(dst - min(dst, [], 'all')) / (max(dst, [], 'all') - min(dst, [], 'all'));
end

