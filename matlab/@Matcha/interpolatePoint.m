function [success, news, newi, newj, newx, newy, sigma, peakval] = interpolatePoint(obj, o, s, i, j)

news    = s;
newi    = i;
newj    = j;
newx    = 0;
newy    = 0;
sigma   = 0;
peakval = 0;

success = false;

tries = 5;
accept_point = false;

while ((tries >= 0) && (~accept_point))
    % Fit quadratic @ (octave=o, scale=s, row=i, column=j)
    [offset, peakval, info] = obj.fitQuadratic(o, news, newi, newj);

    % Check if interpolation was successful
    if (info)            
        % Compute coordinates in the local frame (in respective DoG)
        if ( (offset(1) > 0.6) && (news < size(obj.dog{o},3) - 2) )
            news = news + 1;
        elseif ( (offset(1) < -0.6) && (news > 3) )
            news = news - 1;
        end
        
        if ( (offset(2) > 0.6) && (newj < size(obj.dog{o},2) - 2) )
            newj = newj + 1;
        elseif ( (offset(2) < -0.6) && (newj > 3) )
            newj = newj - 1;
        end
        
        if ( (offset(3) > 0.6) && (newi < size(obj.dog{o},1) - 2) ) 
            newi = newi + 1;
        elseif ( (offset(3) < -0.6) && (newi > 3) )
            newi = newi - 1;
        end
            
        if ( max(abs(offset)) < 0.6 )
            % Compute global coordinates (in original image)
            sigma = (obj.pixelDelta{o,news} / obj.params.seedPixelDelta) * ...
                obj.params.seedSigma * power(2.0, (offset(1) + (news-1))/obj.params.nspo);
            newx = obj.pixelDelta{o,news} * (offset(2) + newj-1) + 1;
            newy = obj.pixelDelta{o,news} * (offset(3) + newi-1) + 1;
            accept_point = true;
            success = true;
            break;
        end
        
        tries = tries - 1;
    else
        break;
    end
end
