function h = showAnaglyph(im1, im2, alpha)
% SHOWANAGLYPH Displays anaglyph of two images
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

anaglyph = zeros(size(im1,1), size(im1,2), 3);

anaglyph(:,:,1) = im2gray(im1);
anaglyph(:,:,2) = im2gray(im2);
anaglyph(:,:,3) = im2gray(im2);

anaglyph = uint8(anaglyph);

h = imshow(anaglyph); hold on;

axis equal;

set(h, 'AlphaData', alpha);

xlabel("x (px) | j (px)");
ylabel("y (px) | i (px)");

end

