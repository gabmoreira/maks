function h = showSideBySide(im1, im2, alpha)
% SHOWSIDEBYSIDE Shows two images side by side
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

im = [im2gray(im1), im2gray(im2)];

h = imshow(uint8(im)); hold on;

axis equal;

set(h, 'AlphaData', alpha);

end

