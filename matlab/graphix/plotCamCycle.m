function plotCamCycle(R, color, labels)
% CAMPLOTCYCLE Plots n cameras on a circle, given absolute rotations.
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

campts=[ 1, -1, 0;
         0,  0, 2;
         1,  1, 0;
        -1, -1, 0;
        -1,  1, 0;
         1,  1, 0;
         1, -1, 0;
        -1, -1, 0;
         0,  0, 2;
        -1,  1, 0;
         1, -1, 0];
     
radius = n;

for i=1:n
    z = ones(11,1)*radius*cos((i-1)*2*pi/n);
    y = ones(11,1)*radius*sin((i-1)*2*pi/n);
    pts = campts*R(i*3-2:i*3,:);
    plot3(pts(:,1),pts(:,2)+y,pts(:,3)+z, 'color', color); hold on;
    
    if(labels)
        text(pts(1,1), pts(1,2)+y(1), pts(1,3)+z(1), ...
            "cam-" + num2str(i), 'color', 'w');
    end
    
end

axis equal;
axis off;

end

