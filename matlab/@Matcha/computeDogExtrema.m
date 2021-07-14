function computeDogExtrema(obj)
% Scale-space extrema detection in this octave

assert(obj.COMPUTED_DOG, 'Must compute Difference-of-Gaussian (DoG) first.');

obj.extrema = [];

fprintf("Matcha: Searching octave [")

% Iterate over all the octaves
for o=1:numel(obj.octave)
    fprintf(" %d ", o);
    % Find DoG points with magnitude larger than threshold
    [ii,jj,ss] = ind2sub(size(obj.dog{o}), find(obj.dog{o} > obj.params.peakThreshold));
    
    % Filter subscripts according to padding
    padding_ii_map = ((ii > obj.params.padding) & (ii < size(obj.dog{o},1) - obj.params.padding));
    padding_jj_map = ((jj > obj.params.padding) & (jj < size(obj.dog{o},2) - obj.params.padding));
    padding_ss_map = ((ss > 1) & (ss < size(obj.dog{o}, 3)));
    
    % Apply filter
    padding_filter = (padding_ii_map & padding_jj_map & padding_ss_map);
    ii = ii(padding_filter);
    jj = jj(padding_filter);
    ss = ss(padding_filter);

    % Iterate over valid points
    for idx=1:numel(ii)
        islocalmin = [true, true, true];
        islocalmax = [true, true, true];
        center     = obj.dog{o}(ii(idx),jj(idx),ss(idx));

        % Look for extrema in the scales ss(idx)-1, ss(idx), ss(idx)+1
        for iter=-1:1:1
            neighbors = obj.dog{o}(ii(idx)-1:ii(idx)+1,jj(idx)-1:jj(idx)+1,ss(idx)+iter);

            neighbors(2,2) = neighbors(2,2) + 0.1;
            if (any(neighbors <= center, 'all'))
                islocalmin(iter+2) = false;
            end

            neighbors(2,2) = neighbors(2,2) - 0.2;
            if (any(neighbors >= center, 'all'))
                islocalmax(iter+2) = false;
            end
        end

        % Store extremum
        if ((all(islocalmax) || all(islocalmin)))
            obj.extrema = [obj.extrema; o, ss(idx), ii(idx), jj(idx)];
        end
    end

end
fprintf("] - %d extrema found.\n", size(obj.extrema, 1));
obj.COMPUTED_DOG_EXTREMA = true;
end

