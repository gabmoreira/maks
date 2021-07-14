function setDefaultParams(obj)
% SETDEFAULTPARAMS Plots n cameras on a circle, given absolute rotations.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira @ isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

% Set border padding
obj.params.padding = 5;

% Minimum image size when generating octaves with subsampling
obj.params.minsize = 2 * obj.params.padding + 2;

% We set the inter pixel distance of the first image (seed) to this
obj.params.seedPixelDelta = 0.5;
% We set the blur sigma for the first image (seed) to this
obj.params.seedSigma = 0.8;

% Maximum number of octaves
obj.params.maxOctave = 5;
            
% Number of scales per octave (nspo)
obj.params.nspo = 3;

% Discard DoG extrema whose abs value is smaller than this
obj.params.peakThreshold  = 255 * 0.04 / 3.0;

% Thresholding for edge identification and discarding
obj.params.edgeThreshold  = 0.08;
obj.params.edgeThreshold1 = 0.08;

% Orientation assignament params: Number of bins / histogram
obj.params.orientationBins       = 36;
obj.params.orientationSigma      = 1.5;
obj.params.orientationHistThresh = 0.8;
obj.params.orientationLambda     = 1.5;

% Descriptor parameters: Patch size
obj.params.descriptorLambda = 6;

% Descriptor parameters: nist x nist histograms per keypoint
obj.params.nhist = 4; %

% Descriptor parameters: Number of bins per each of the histograms
obj.params.nori = 8; 
end

