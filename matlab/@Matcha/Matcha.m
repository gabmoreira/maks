classdef Matcha < handle
% MATCHA Feature extraction using HoG
%
% Other m-files required: @Matcha/*
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

    properties (Access = public)
        params;
        
        % Hold original uint8 RGB
        rgb;
        
        % Holds the grayscale input image
        inputImage;
        
        % Array holding octaves each with nspo+3 scales
        octave;
        
        % Array holding octaves each with nspo+2 DoG
        dog;  
        
        % Camera
        camera
        exif;
        
        pixelDelta;
        sigma;
                
        keys;
        extrema;
                
        % Status flags
        COMPUTED_DOG;
        COMPUTED_DOG_EXTREMA;
        ASSIGNED_ORIENTATIONS;
        COMPUTED_DESCRIPTIONS;
        COMPUTED_MATCHES;
    end
    
    
    methods (Access = public)
        function obj = Matcha(im, varargin)
            
            validIm         = @(x) isnumeric(x);
            validInputSigma = @(x) isscalar(x);
            
            defaultInputSigma = 0.8;
            
            p = inputParser;
            addRequired(p,'im', validIm);
            addOptional(p,'inputSigma', defaultInputSigma, validInputSigma);
            
            parse(p,im,varargin{:});

            obj.params.inputSigma = p.Results.inputSigma;
            
            % Assumed blur level of input image
            if ( abs(obj.params.inputSigma - defaultInputSigma) > 1e-2 )
                fprintf("Matcha: Input image will be resized\n");
            else
                fprintf("Matcha: Input image will not be resized\n");
            end

            obj.setDefaultParams();
            
            obj.load(p.Results.im);

            obj.COMPUTED_DOG          = false;
            obj.COMPUTED_DOG_EXTREMA  = false;
            obj.ASSIGNED_ORIENTATIONS = false;
            obj.COMPUTED_DESCRIPTIONS = false;
            obj.COMPUTED_MATCHES      = false;
        end
                
        load(obj, filename);
        
        setDefaultParams(obj);
        
        im = bilinearInterp(obj, seed, interpixelDistance);

        computePyramid(obj);
        
        computeFeatures(obj);
        
        computeOrientations(obj);
        
        computeDogExtrema(obj);
          
        [inlierPts1, inlierPts2, F, pts1, pts2] = match(obj, im2, ratio, pnorm, maxiter, errorThresh, inlierThresh)

        [h, size] = gaussianKernel(obj, sflag, std, size);
                 
        [success, s, i, j, x, y, sigma, peakval] = interpolatePoint(obj, o, s, i, j);
        
        interpolateExtrema(obj);
        
        [offset, peakval, info] = fitQuadratic(obj, dogs, s, i, j);
        
        inEdge = isEdge(obj, o, s, i, j);
        
        filterEdges(obj);
        
        showExtrema(obj);
        
        varargout = epipolarRansac(obj, pts1, pts2, varargin);
        
        F = epipolarGeometry(obj, pts1, pts2);
    end
end

