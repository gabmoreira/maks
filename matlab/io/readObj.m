function dst = readObj(filename)
% READOBJ Reads object (.obj) file
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:        Gabriel Moreira
% Website:       https://www.github.com/gabmoreira/maks
% Last revision: 14-July-2021

v    = []; 
vt   = []; 
vn   = []; 
f.v  = [];
f.vt = [];
f.vn = [];

fid = fopen(filename);

while 1    
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    
    % Line type 
    ln = sscanf(tline,'%s',1);

    switch ln
        % Mesh vertexs
        case 'v'   
            v = [v; sscanf(tline(2:end),'%f')'];
        
        % Texture
        case 'vt'  
            vt = [vt; sscanf(tline(3:end),'%f')'];
            
        % Normal
        case 'vn' 
            vn = [vn; sscanf(tline(3:end),'%f')'];
            
        % Face   
        case 'f'  
            fv = []; fvt = []; fvn = [];
            str = textscan(tline(2:end),'%s'); str = str{1};
       
            nf = length(strfind(str{1},'/')); 
            [tok, str] = strtok(str,'//');     
            for k = 1:length(tok) 
                fv = [fv str2double(tok{k})]; 
            end
           
            if (nf > 0) 
                [tok, str] = strtok(str,'//');  
                for k = 1:length(tok) 
                    fvt = [fvt str2double(tok{k})]; 
                end
            end
            
            if (nf > 1) 
                [tok, str] = strtok(str,'//');
                for k = 1:length(tok)
                    fvn = [fvn str2double(tok{k})]; 
                end
            end
            f.v  = [f.v; fv];
            f.vt = [f.vt; fvt];
            f.vn = [f.vn; fvn];
    end
end

fclose(fid);

dst.v  = v;
dst.vt = vt;
dst.vn = vn; 
dst.f  = f;

end
