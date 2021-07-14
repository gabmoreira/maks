function computePyramid(obj)
% COMPUTEPYRAMID Compute Gaussian and DoG pyramid
%
% Other m-files required: Matcha.m
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

% Create seed image by interpolating the input image with
% a pixel distance of 0.5 and then blurring it to bring the 
% blur to obj.params.seedSigma
if (obj.params.inputSigma < obj.params.seedSigma)
    fprintf("Matcha: Interpolating image to double its size\n");
    seed   = bilinearInterp(obj.inputImage, 0.5);
    kernel = obj.gaussianKernel(false, (1.0/obj.params.seedPixelDelta)*sqrt(obj.params.seedSigma^2 - obj.params.inputSigma^2), -1);
    seed   = conv2(seed, kernel, 'same');
else
    fprintf("Matcha: Image interpolation not performed\n");
    seed = obj.inputImage;
    obj.params.seedPixelDelta = 2 * obj.params.seedPixelDelta;
end

[n,m] = size(seed);

% Each octave contains nspo scales plus three other blurs
obj.octave{1}        = zeros(n, m, obj.params.nspo + 3);
obj.octave{1}(:,:,1) = seed;

obj.pixelDelta{1,1} = obj.params.seedPixelDelta;
obj.sigma{1,1}      = obj.params.seedSigma;

obj.dog{1} = zeros(n, m, obj.params.nspo + 2);


% Octave index. Starting at the first octave
o = 1;

% Iterate until we hit maximum number of octaves or images get too small
while (o < obj.params.maxOctave)
    fprintf("Matcha: Creating octave %d\n", o);
    fprintf("  Scale 1 | kernel size -- | sigma %2.2e | pixel delta %2.2e\n", ...
        obj.sigma{o,1}, obj.pixelDelta{o,1});
                        
    % Form each level by adding incremental blur from previous level.
    for s=2:size(obj.octave{o}, 3)
        
        obj.pixelDelta{o,s} = obj.params.seedPixelDelta * power(2.0, o-1);
        obj.sigma{o,s}      = (obj.pixelDelta{o,s} / obj.params.seedPixelDelta) * ...
                                obj.params.seedSigma * power(2, (s-1) / obj.params.nspo);
    
        sigmaIncrease = (obj.params.seedSigma / obj.params.seedPixelDelta) * ...
            sqrt(power(2.0, 2.0 * (s-1) / obj.params.nspo) - power(2.0, 2.0 * (s-2) / obj.params.nspo));

        kernel = obj.gaussianKernel(false, sigmaIncrease, -1);

        obj.octave{o}(:,:,s) = conv2(obj.octave{o}(:,:,s-1), kernel, 'same');

        fprintf("  Scale %d | kernel size %2d | sigma %2.2e | pixel delta %2.2e\n", ...
            s, size(kernel,1), obj.sigma{o,s}, obj.pixelDelta{o,s});
    end

    % Compute an array, dogs, of difference-of-Gaussian images by
    % subtracting each image from its next blurred version
    for s=1:size(obj.dog{o}, 3)
        obj.dog{o}(:,:,s) = obj.octave{o}(:,:,s+1) - obj.octave{o}(:,:,s);
    end

    if ( (floor(size(obj.octave{o}(:,:,1),1)) / 2 > obj.params.minsize) && ...
         (floor(size(obj.octave{o}(:,:,1),2)) / 2 > obj.params.minsize) && (o+1 < obj.params.maxOctave) )

        % Allocate next octave
        obj.octave{o+1} = zeros(round(size(obj.octave{o},1) / 2), round(size(obj.octave{o},2) / 2), obj.params.nspo+3);
        % Set seed for the next octave as the frame nspo of current octave
        obj.octave{o+1}(:,:,1) = obj.octave{o}(1:2:end,1:2:end,obj.params.nspo+1);

        % Allocate next DoG
        obj.dog{o+1} = zeros(round(size(obj.octave{o},1) / 2), round(size(obj.octave{o},2) / 2), obj.params.nspo+2);
  
        % Store inter pixel distance and sigma
        obj.pixelDelta{o+1,1} = obj.params.seedPixelDelta * power(2.0, o);
        obj.sigma{o+1,1}      = (obj.pixelDelta{o+1,1} / obj.params.seedPixelDelta) * obj.params.seedSigma;
    
        % Move to next octave
        o = o + 1;
    else
        break;
    end
end
obj.COMPUTED_DOG = true;
fprintf("Matcha: Pyramid complete.\n");
end

