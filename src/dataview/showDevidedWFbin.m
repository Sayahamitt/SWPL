nrows = 256;
numSubmat = 36;
sizeCompDouble = 16;

dirpath = '\\filesrv\Personal\3B\jt3b1215\working\SWPL\test\parares\';
fnamesuffix = 'th_obs.bin';

entireResult = [];

for ii=0:(numSubmat-1)
    fpath = [dirpath,num2str(ii),fnamesuffix];
    fileinfo = dir(fpath);
    fileID = fopen(fpath);
    ncols = fileinfo.bytes/(sizeCompDouble*nrows);
    
    same_real = fread(fileID, [nrows, ncols], 'double', 8);
    fseek(fileID, 8, 'bof');
    same_imag = fread(fileID, [nrows, ncols], 'double', 8);
    fclose(fileID);
    same_z = complex(same_real, same_imag);
    
    entireResult = horzcat(entireResult,same_z);
end

figure;surf(abs(entireResult),'linestyle','none');