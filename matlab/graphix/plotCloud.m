function plotCloud(cloud, s)
% PLOTCLOUD Plots point cloud
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

figure();

set(gcf,'Color',[0 0 0]);

scatter3(cloud.pts(:,1), cloud.pts(:,2), cloud.pts(:,3), ...
    s*ones(size(cloud.pts,1),1), double(cloud.color) ./ 255.0, 'filled');
hold on;

xlabel('x');
ylabel('y');
zlabel('z');

grid on;

ax = gca;
set(gca,'color',[0 0 0]);
ax.GridColor = [1, 1, 1];  % [R, G, B]
ax.GridAlpha = 0.3;
axis equal;

end

