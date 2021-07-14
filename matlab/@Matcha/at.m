function pixels = at(obj, location)

assert(ismatrix(location), 'Location must be a matrix');
assert(size(location,2) > 1, 'Location must be n x 2');

pixels = zeros(size(location,1),3);

for i=1:size(location,1)
    pixels(i,:) = [obj.rgb(floor(location(i,2)), floor(location(i,1)),1), ...
                   obj.rgb(floor(location(i,2)), floor(location(i,1)),2), ...
                   obj.rgb(floor(location(i,2)), floor(location(i,1)),3)];
               
end
end

