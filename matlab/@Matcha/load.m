function load(obj, im)

obj.rgb  = im;
obj.inputImage = double(im2gray(obj.rgb));

end

