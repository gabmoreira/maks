function computeFeatures(obj)

assert(obj.COMPUTED_DOG_EXTREMA, 'Matcha: Must compute Difference-of-Gaussian (dog) extrema first.');
assert(obj.ASSIGNED_ORIENTATIONS, 'Matcha: Must assign orientations first.');

fprintf("Matcha: Building descriptions   - [")
obj.keys.x        = zeros(size(obj.extrema,1), 1);
obj.keys.y        = zeros(size(obj.extrema,1), 1);
obj.keys.sigma    = zeros(size(obj.extrema,1), 1);
obj.keys.theta    = zeros(size(obj.extrema,1), 1);
obj.keys.features = zeros(size(obj.extrema,1), obj.params.nhist * obj.params.nhist * obj.params.nori);

counter = 1;

consoleCounter = 0;

% Iterate over keypoints
for k=1:size(obj.extrema, 1)
    if (round(100*k/size(obj.extrema, 1)) > consoleCounter)
        fprintf("=");
        consoleCounter = consoleCounter + 10;
    end
    
    % Placeholder for nhist x nhist histograms each with nori bins;
    histdata = zeros(obj.params.nhist, obj.params.nhist, obj.params.nori);
    
    % Keypoint information
    okey     = obj.extrema(k,1);              % Octave (o)
    skey     = obj.extrema(k,2);              % Scale (s)
    xkey     = obj.extrema(k,5);              % Absolute coordinates in reference frame of input
    ykey     = obj.extrema(k,6);              % Absolute coordinates in reference frame of input
    sigmakey = obj.extrema(k,7);              % Blur level
    thetakey = obj.extrema(k,9);              % Keypoint orientation (rad)
    deltakey = obj.pixelDelta{okey, skey};    % Inter pixel distance in this octave

    % Look at pixels within radius around the point
    radius = 1.4142 * obj.params.descriptorLambda * sigmakey * (obj.params.nhist+1.0)/obj.params.nhist;
    sigmakey2 = sigmakey * sigmakey;
    
    % Boundary in the coordinate system of the octave 
    mmin = round((xkey - radius) / deltakey);
    mmax = round((xkey + radius) / deltakey);
    nmin = round((ykey - radius) / deltakey);
    nmax = round((ykey + radius) / deltakey);
    
    % Check if patch is within image
    if ( (mmin > 1) && (mmax < size(obj.octave{okey},2)-1) && (nmin > 1) && (nmax < size(obj.octave{okey},1)-1))
        % Search inside patch. m and n are in the octave reference frame
        for m = mmin:mmax
            for n = nmin:nmax
                % Normalized coordinates of the keypoint in the reference 
                % frame of the input image
                xhat = ((m*deltakey - xkey)*cos(thetakey) + (n*deltakey - ykey)*sin(thetakey)) / sigmakey;
                yhat = ((xkey - m*deltakey)*sin(thetakey) + (n*deltakey - ykey)*cos(thetakey)) / sigmakey;

                if ( max(abs(xhat), abs(yhat)) < obj.params.descriptorLambda*(obj.params.nhist-1)/obj.params.nhist )
                    % Compute horizontal derivative @ (okey, skey, n, m)
                    ddx = 0.5 * (obj.octave{okey}(n,m+1,skey) - obj.octave{okey}(n,m-1,skey));

                    % Vertical derivative @ (okey, skey, n, m)
                    ddy = 0.5 * (obj.octave{okey}(n+1,m,skey) - obj.octave{okey}(n-1,m,skey));

                    % Gradient norm @ (okey, skey, n, m)
                    gradnorm = sqrt(ddx^2 + ddy^2);

                    % Normalize gradient orientation;
                    thetahat = mod(atan2(-ddy, ddx) - thetakey, 2*pi);

                    % Squared distance from (n,m) to keypoint in input image
                    % coordinates
                    distsq = (m*deltakey - xkey)^2 + (n*deltakey - ykey)^2;

                    % Weight of this pixel in particular
                    c = exp(-distsq/(2*sigmakey2*obj.params.descriptorLambda^2) ) * gradnorm;

                    aux    = 0.5 * obj.params.nhist / obj.params.descriptorLambda;
                    auxinv = 1.0 / aux;

                    % Update nhist x nhist histograms
                    for i=1:obj.params.nhist
                        % Update nhist x nhist histograms
                        for j=1:obj.params.nhist
                            
                            xhati = (j-(1+obj.params.nhist)/2) * auxinv;
                            yhati = (i-(1+obj.params.nhist)/2) * auxinv;
                            if ( (abs(xhati-xhat)<=auxinv) && (abs(yhati-yhat)<=auxinv) )
                                % Iterate over #nori bins
                                for bin=1:obj.params.nori
                                    thetabin = 2*pi*bin/obj.params.nori;
                                    if ( abs(mod(thetahat-thetabin, 2*pi)) < (2*pi / obj.params.nori) )
                                        factor1 = 1 - aux * abs(xhat-xhati);
                                        factor2 = 1 - aux * abs(yhat-yhati);
                                        factor3 = 1 - (obj.params.nori/(2*pi)) * abs(mod(thetahat-thetabin, 2*pi));
                                        histdata(i,j,bin) = histdata(i,j,bin) + factor1 * factor2 * factor3 * c;
                                    end
                                end
                            end
                        end
                    end
                   
                end
            end
        end
    end
    % Threshold
    normf = norm(histdata(:));
    if (normf > 1e-6)
        features = histdata(:);
        features(features > 0.2*normf) = 0.2*normf;
        features = round(512 * features / norm(features));
        obj.keys.features(counter,:) = features;
        obj.keys.x(counter)          = xkey;
        obj.keys.y(counter)          = ykey;
        obj.keys.sigma(counter)      = sigmakey;
        obj.keys.theta(counter)      = thetakey;
        counter = counter + 1;
    end
end
fprintf("] - Complete.\n");

% Trim vectors
obj.keys.x        = obj.keys.x(1:counter-1);
obj.keys.y        = obj.keys.y(1:counter-1);
obj.keys.sigma    = obj.keys.sigma(1:counter-1);
obj.keys.theta    = obj.keys.theta(1:counter-1);
obj.keys.features = obj.keys.features(1:counter-1,:);

obj.COMPUTED_DESCRIPTIONS = true;
end