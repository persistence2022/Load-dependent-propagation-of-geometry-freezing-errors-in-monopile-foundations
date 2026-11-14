function h=hash_file(path)
fid=fopen(path,'rb');assert(fid>=0,'stage:MissingFile','Cannot open %s',path);cl=onCleanup(@()fclose(fid));b=fread(fid,inf,'*uint8');h=stage.hash_bytes(b);
end
