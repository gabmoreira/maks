function showMatches(im1, pts1, im2, pts2, type, alpha)
% SHOWMATCHES Displays matches between two images
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

im1 = im2gray(im1);
im2 = im2gray(im2);

figure();

switch type
    case 'anaglyph'
        showAnaglyph(im1, im2, alpha);

        for i=1:size(pts1,1)
            plot(pts1(i,1), pts1(i,2), 'o', 'color', [0,1,0], 'linewidth', 1);
        end

        for i=1:size(pts1,1)
            plot(pts2(i,1), pts2(i,2), 'o', 'color', [1,0,0], 'linewidth', 1);
        end

        for i=1:size(pts1,1)
            s = plot([pts1(i,1),pts2(i,1)], [pts1(i,2),pts2(i,2)], '-', 'color', [1,1,0], 'linewidth', 1.2);
            s.Color(4) = 0.25;
        end

    case 'side'
        showSideBySide(im1, im2, alpha);

        for i=1:size(pts1,1)
            plot(pts1(i,1), pts1(i,2), 'o', 'color', [0,1,0], 'linewidth', 1);
        end

        for i=1:size(pts2,1)
            plot(pts2(i,1)+size(im1,2), pts2(i,2), 'o', 'color', [1,0,0], 'linewidth', 1);
        end

        for i=1:size(pts1,1)
            s = plot([pts1(i,1),pts2(i,1)+size(im1,2)], [pts1(i,2),pts2(i,2)], '-', 'color', [1,1,0], 'linewidth', 1.2);
            s.Color(4) = 0.25;
        end
        
    otherwise
        error('Unknown display type');
end

