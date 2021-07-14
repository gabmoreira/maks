function load(obj, im)
% LOAD Load RGB or grayscale image
%
% Other m-files required: Matcha.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

obj.rgb  = im;
obj.inputImage = double(im2gray(obj.rgb));

end

