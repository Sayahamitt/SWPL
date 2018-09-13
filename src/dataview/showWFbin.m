nrows = 128;
ncols = 32;
%fileID = fopen('..\lastr01.bin');
%fileID = fopen('..\ei_distination.bin');
%fileID = fopen('source.bin');
fileID = fopen('\\filesrv\Personal\3B\jt3b1215\working\SWPL\test\0th_obs.bin');

same_real = fread(fileID, [nrows, ncols], 'double', 8);
fseek(fileID, 8, 'bof');
same_imag = fread(fileID, [nrows, ncols], 'double', 8);
fclose(fileID);

same_z = complex(same_real, same_imag);

figure;surf(abs(same_z),'linestyle','none');