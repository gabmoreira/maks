function plotCam(R, t, scale, color)
% PLOTCAM Plots n cameras on a circle, given absolute rotations.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

n = size(R,1) / 3;

campts = [ 1, -1, 1;
           0,  0, 0;
           1,  1, 1;
          -1, -1, 1;
          -1,  1, 1;
           1,  1, 1;
           1, -1, 1;
          -1, -1, 1;
           0,  0, 0;
          -1,  1, 1;
           1, -1, 1];

for i=1:n
    pts = scale*campts*R(i*3-2:i*3,:)';
    
    plot3(pts(:,1)+t(i*3-2),pts(:,2)+t(i*3-1),pts(:,3)+t(i*3), ...
        '-', 'Color', color); hold on;
    
end
axis equal;

end

