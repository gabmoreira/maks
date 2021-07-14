function F = epipolarGeometry(pts1, pts2)
% EPIPOLARGEOMETRY Least-squares fundamental matrix from stereo matches.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

npts = size(pts1,1);
if (npts < 9)
    error("Insufficient number of points");
end

W = zeros(npts,9);

for i=1:npts
    W(i,1) = pts1(i,1) * pts2(i,1);
    W(i,2) = pts1(i,1) * pts2(i,2);
    W(i,3) = pts1(i,1);
    W(i,4) = pts1(i,2) * pts2(i,1);
    W(i,5) = pts1(i,2) * pts2(i,2);
    W(i,6) = pts1(i,2);
    W(i,7) = pts2(i,1);
    W(i,8) = pts2(i,2);
    W(i,9) = 1;
end

[~,~,V] = svd(W);

F = reshape(V(:,end), 3, 3)';
F = rankReduction(F,2);

end

