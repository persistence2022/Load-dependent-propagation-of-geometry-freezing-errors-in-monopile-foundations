function snapshot_execution(root,out)

dest=fullfile(out,'source_snapshot');mkdir(dest);mkdir(fullfile(dest,'pisa_model'));
files=[dir(fullfile(root,'pisa_model','*.m'));dir(fullfile(root,'code','*.m'))];rows={};
inputnames={fullfile(root,'sources','design_benchmarks','published_design_curves.csv'),...
 fullfile(root,'soil_profiles','offshore_profiles_0p05m.csv'),...
 fullfile(root,'pisa_model','core_manifest.json')};
for j=1:numel(files)
 f=fullfile(files(j).folder,files(j).name);if contains(files(j).folder,[filesep 'pisa_model']),d=fullfile(dest,'pisa_model',files(j).name);else,d=fullfile(dest,files(j).name);end
 copyfile(f,d);rows(end+1,:)={f,sha(f),files(j).bytes};
end
for j=1:numel(inputnames),f=inputnames{j};copyfile(f,fullfile(dest,extractAfter(f,find(f==filesep,1,'last'))));info=dir(f);rows(end+1,:)={f,sha(f),info.bytes};end
writetable(cell2table(rows,'VariableNames',{'source_path','sha256','bytes'}),fullfile(out,'input_code_manifest.csv'));
addpath(fullfile(dest,'pisa_model'),'-begin');
end
function s=sha(f),fid=fopen(f,'r');c=onCleanup(@()fclose(fid));b=fread(fid,inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);s=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));end
