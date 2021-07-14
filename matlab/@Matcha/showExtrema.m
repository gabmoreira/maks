function showExtrema(obj)
% SHOWEXTREMA Shows DoG extrema
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira @ isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

imagesc(obj.inputImage); 
colormap gray; 
hold on; 

plot(obj.extrema(:,5), obj.extrema(:,6), '.', 'color', [0,1,0]);
axis equal;

for i=1:size(obj.keys, 1)
    x1 = obj.extrema(i,5);
    x2 = x1 + 10 * sin(obj.extrema(i,9));

    y1 = obj.extrema(i,6);
    y2 = y1 + 10 * cos(obj.extrema(i,9));
    plot([x1,x2], [y1,y2], '-', 'color', [0,1,0]);
end

xlabel("x (px) | j (px)");
ylabel("y (px) | i (px)");

title("Matcha: Keypoints and orientations");

end

