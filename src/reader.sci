function mat = readComplexDouble(fname,nrows,ncols)
	[fd,err] = mopen(fname,'rb');
	if err~=0 then
		disp("Faild to open the file.")
		return;
	end
	fsize = fileinfo(fname);
	fsize = fsize(1);

	beforConstuctData = mget(fsize,'dl',fd);
	mclose(fd);

	reals = beforConstuctData(modulo((1:size(beforConstuctData,2)),2)==1);
	imgs = beforConstuctData(modulo((1:size(beforConstuctData,2)),2)==0);
	comps = reals+imgs*%i; 

	mat = matrix(comps,nrows,ncols)
endfunction