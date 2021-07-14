function [bestModel,bestError,inliers] = ransac3D(xyz1, xyz2, params)
% RANSAC3D Implements RANSAC to find a rigid transformation.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% email:         gmoreira (at) isr.tecnico.ulisboa.pt
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

inliers           = [];
minMatches        = 4;
bestError         = realmax;
invalidMatchesIdx = union(find(xyz1(:,3)==0), find(xyz2(:,3)==0));
validMatchesIdx   = setdiff(1:size(xyz1,1), invalidMatchesIdx);

% Default best transformation estimate
bestModel.R = eye(3);
bestModel.t = zeros(1,3);

% RANSAC main loop
for i = 1:params.maxiter
    
    % Selects 4 random matches 
    maybeInliersIdx = validMatchesIdx(randperm(numel(validMatchesIdx), minMatches));
    
    % Finds rigid transformation from pts2 to pts1
    model = generalProcrustes(xyz1(maybeInliersIdx,:), xyz2(maybeInliersIdx,:)); 

    % Indices which were not selected previously
    otherMatchesIdx = setdiff(validMatchesIdx, maybeInliersIdx);

    % See if these other matches verify the model 
    delta = xyz2(otherMatchesIdx,:) - (xyz1(otherMatchesIdx,:) * model.R + model.t');
                                   
    errors = vecnorm(delta,2,2);
    alsoInliersIdx = otherMatchesIdx(errors < params.errorThreshold);

    % If they are in enough numbers we accept the matches 
    if (numel(alsoInliersIdx) > params.minMatches)
        
        % Compute model with maybe inliers and also inliers
        currentModelIdx  = union(alsoInliersIdx, maybeInliersIdx);
        
        [betterModel, e] = generalProcrustes(xyz1(currentModelIdx,:), xyz2(currentModelIdx,:)); 

        % Update the best model and error
        if e < bestError
            bestModel = betterModel;
            bestError = e;
            inliers   = currentModelIdx;   
            % Verbose
            if params.verbose
                fprintf("RMSE: %.5e | No. of inliers: %d\n", bestError, numel(currentModelIdx));
            end
        end
    end    
end

end

