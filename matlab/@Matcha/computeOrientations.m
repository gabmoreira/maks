function computeOrientations(obj)

assert(obj.COMPUTED_DOG_EXTREMA, 'Matcha: Must compute Difference-of-Gaussian (dog) extrema first.');

if (size(obj.extrema, 2) ~= 9)
    obj.extrema = [obj.extrema(:,1:8), nan * ones(size(obj.extrema,1),1)];
end

consoleCounter = 0;
fprintf("Matcha: Assigning orientations  - [")
for k=1:size(obj.extrema,1)
    if (round(100*k/size(obj.extrema, 1)) > consoleCounter)
        fprintf("=");
        consoleCounter = consoleCounter + 10;
    end
    
    % Get keypoint information
    o = obj.extrema(k,1);       % Octave (o)
    s = obj.extrema(k,2);       % Scale (s)
    i = obj.extrema(k,3);       % Row (i or y)
    j = obj.extrema(k,4);       % Column (j or x)
    sigma = obj.extrema(k,7);   % Blur level
    
    h = size(obj.dog{o}, 1); % No. rows (height)
    w = size(obj.dog{o}, 2); % No. cols (width)

    % Allocation histogram vector
    histbins = zeros(obj.params.orientationBins,1);

    % Look at pixels within 3 sigma around the point and sum their
    % Gaussian weighted gradient magnitudes into the histogram.
    radius  = round(sigma * 3.0 * obj.params.orientationLambda);
    radius2 = radius * radius;
    sigma2  = sigma*sigma;
    
    % Check if keypoint is distant enough from borders
    if ( (i >= 3*obj.params.orientationLambda*sigma) && (i <= h - 3*obj.params.orientationLambda*sigma) && ...
         (j >= 3*obj.params.orientationLambda*sigma) && (j <= w - 3*obj.params.orientationLambda*sigma)  )
        % Patch around the keypoint
        patch_imin = max(2,i-radius);
        patch_imax = min(i+radius,h-2);
        patch_jmin = max(2,j-radius);
        patch_jmax = min(j+radius,w-2);

        % Iterative search inside patch
        for ii=patch_imin:patch_imax
            for jj=patch_jmin:patch_jmax
                % Gradient magnitude  
                ddx = 0.5 * (obj.octave{o}(ii,jj+1,s) - obj.octave{o}(ii,jj-1,s));
                ddy = 0.5 * (obj.octave{o}(ii+1,jj,s) - obj.octave{o}(ii-1,jj,s));

                gval = sqrt(ddx^2 + ddy^2);

                % Squared distance to the center of the patch
                distsq = (ii - i)^2 + (jj - j)^2;

                if ( (gval > 0.0)  &&  (distsq < radius2 + 0.5) )
                    % Sample contribution: gradient magnitude weighted
                    contribution = exp(- distsq / (2 * sigma2) ) * gval;

                    % Bin index
                    bin = round( ((obj.params.orientationBins-1) / (2*pi)) * wrapTo2Pi(atan2(-ddy, ddx)) );

                    % Add to histogram
                    histbins(bin+1) = histbins(bin+1) + contribution;
                end
            end
        end

        % Apply smoothing 6 times for accurate Gaussian approximation.
        for iter=1:6
            histbins = conv(histbins, [1.0/3.0; 1.0/3.0; 1.0/3.0;], 'same');
        end

        % Find maximum value in histogram
        maxval = max(histbins);

        % Get indices of peaks of histogram 
        idx_peaks = find(islocalmax(histbins));

        % If value is within
        % par.OriHistThresh of maximum value, then generate a keypoint.
        for ii=1:numel(idx_peaks)
            i = idx_peaks(ii);
            % Use parabolic fit to interpolate peak location from 3 samples.
            % Set angle in [-PI, PI]

            if (histbins(i) >= obj.params.orientationHistThresh * maxval)
                prev = i - 1;
                next = i + 1;

                if (i == 1)
                    prev = numel(idx_peaks);
                end

                if (i == numel(idx_peaks))
                    next = 1;
                end

                interp = interpPeak(histbins(prev), histbins(i), histbins(next));
                orientation = 2.0 * pi * (i + 0.5 + interp) / obj.params.orientationBins - pi;

                obj.extrema(k,end) = orientation;
            end
        end
    end
end

% Remove keypoints without orientation
obj.extrema(isnan(obj.extrema(:,end)),:) = [];
fprintf("] - %d orientations assigned.\n", size(obj.extrema,1));
obj.ASSIGNED_ORIENTATIONS = true;
end

function peak = interpPeak(a, b, c)
if (b < 0.0)
    a = -a;
    b = -b;
    c = -c;
end
    
peak = 0.5 * (a - c) / (a - 2.0 * b + c);
end