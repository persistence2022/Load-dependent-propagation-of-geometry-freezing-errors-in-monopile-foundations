function [runDir,meta,resolution]=resolve_current_run(pointer,projectRoot)

if nargin<2,projectRoot=fileparts(fileparts(mfilename('fullpath')));end
pointer=char(pointer);projectRoot=char(projectRoot);
if ~isfile(pointer),pointer=fullfile(projectRoot,pointer);end
assert(isfile(pointer),'Model:MissingPointer','Run pointer does not exist: %s',pointer);
meta=jsondecode(fileread(pointer));parent=fileparts(pointer);
candidates={};labels={};
for key={'run_relative','run_dir_relative','relative_run'}
 if isfield(meta,key{1}) && ~isempty(meta.(key{1}))
  value=char(meta.(key{1}));value=strrep(value,char(92),'/');
  candidates{end+1}=fullfile(projectRoot,value);labels{end+1}=['relative:' key{1}];
 end
end
recorded='';
for key={'run','run_dir','path'}
 if isfield(meta,key{1}) && ~isempty(meta.(key{1}))
  recorded=char(meta.(key{1}));break;
 end
end
if ~isempty(recorded)
 normalized=strrep(recorded,char(92),'/');
 absolute=~isempty(regexp(normalized,'^(?:[A-Za-z]:/|/|//)','once'));
 if ~absolute
  candidates=[candidates,{fullfile(parent,normalized),fullfile(projectRoot,normalized)}];
  labels=[labels,{'pointer-relative legacy path','project-relative legacy path'}];
 end

 token=regexp(normalized,'(?:^|/)((?:field_validation|design_baselines)/runs/[^/]+)$','tokens','once');
 if ~isempty(token)
  candidates{end+1}=fullfile(projectRoot,token{1});labels{end+1}='relocated project subdirectory';
 end
 token=regexp(normalized,'(?:^|/)(results/[^/]+)$','tokens','once');
 if ~isempty(token)
  candidates{end+1}=fullfile(projectRoot,token{1});labels{end+1}='relocated results suffix';
 end
 token=regexp(normalized,'(?:^|/)(runs/[^/]+)$','tokens','once');
 if ~isempty(token)
  candidates{end+1}=fullfile(parent,token{1});labels{end+1}='relocated pointer-folder runs suffix';
 end
 candidates{end+1}=recorded;labels{end+1}='recorded path fallback';
end
for k=1:numel(candidates)
 if isfolder(candidates{k})
  runDir=char(java.io.File(candidates{k}).getCanonicalPath());
  resolution=labels{k};meta.recorded_run=recorded;meta.run=runDir;
  return
 end
end
error('Model:MissingRun','No accessible run for %s. Keep its results/runs folder or set run_relative relative to the project root. Recorded path: %s',pointer,recorded);
end
