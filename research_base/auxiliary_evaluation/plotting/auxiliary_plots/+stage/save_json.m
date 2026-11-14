function save_json(path,s)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);cl=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(s,'PrettyPrint',true),'char');
end
