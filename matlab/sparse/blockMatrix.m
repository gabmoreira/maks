function dst = blockMatrix(src, srcStride, ij, varargin)
% BLOCKMATRIX Build sparse block matrix from blocks in src and edges in ij
%
%       'src' - Sparse p x qn matrix containing n blocks of size p x q
% 'srcStride' - scalar specifying the row stride between blocks
%        'ij' - 2 x n matrix whose edges contain the nonzero entry indices
% 'symmetric' - Flag to symmetrize the output
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

% Input parser
validSrc         = @(x) isnumeric(x) && ismatrix(x);
validSrcStride   = @(x) isscalar(x);
validij          = @(x) isnumeric(x) && ismatrix(x);

validSymmetric   = @(x) islogical(x);
defaultSymmetric = false;

p = inputParser;
addRequired(p, 'src',       validSrc);
addRequired(p, 'srcStride', validSrcStride);
addRequired(p, 'ij',        validij);
addOptional(p, 'symmetric', defaultSymmetric, validSymmetric);

parse(p, src, srcStride, ij, varargin{:});

% Size of the output matrix
if (p.Results.symmetric)
    dstH = max(ij(:)) * srcStride;
    dstW = dstH;
else
    dstH = max(ij(1,:)) * srcStride;
    dstW = max(ij(2,:)) * srcStride;
end

% Size of the blocks
blockH = size(src,1);
blockW = srcStride;

numBlocks = size(src,2) / srcStride;

% Allocate vectors
% Twice the amount of entries for symmetric matrix
if (p.Results.symmetric)
    entries = zeros(numel(src)*2, 1);
    ei      = zeros(numel(src)*2, 1);
    ej      = zeros(numel(src)*2, 1);
else
    entries = zeros(numel(src), 1);
    ei      = zeros(numel(src), 1);
    ej      = zeros(numel(src), 1);
end

k=1;
for i=1:blockH
    for j=1:blockW
        entries(k:k+numBlocks-1) = src(i,j:srcStride:end)';
        ei(k:k+numBlocks-1)      = ij(1,:)' * blockW - blockW + i;
        ej(k:k+numBlocks-1)      = ij(2,:)' * blockW - blockW + j;
        k = k + numBlocks;
    end
end

% Duplicate entries if symmetric flag is true
if (p.Results.symmetric)
    entries(numel(src)+1:end) = entries(1:numel(src));
    ei(numel(src)+1:end)      = ej(1:numel(src));
    ej(numel(src)+1:end)      = ei(1:numel(src));
end

% Create sparse block matrix from triplets
dst = sparse(ei, ej, double(entries), dstH, dstW);

end

