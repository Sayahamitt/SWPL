fid = fopen('\\filesrv\Personal\3B\jt3b1215\working\SWPL\test\0th_axis.bin');

vec = fread(fid, 256, 'double');
fclose(fid);

plot(vec,"-*");