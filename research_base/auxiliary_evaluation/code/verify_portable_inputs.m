function verify_portable_inputs
here=fileparts(mfilename('fullpath'));review=fileparts(here);addpath(fullfile(here,'prediction'));
manifest=readtable(fullfile(review,'inputs','input_manifest.csv'),'TextType','string');
for j=1:height(manifest)
 target=fullfile(review,'inputs',manifest.relative_path(j));
 assert(isfile(target)&&strcmp(stage.hash_file(target),manifest.sha256(j)),'review:InputIntegrity','Copied input hash mismatch');
end
resolved=which('run_evaluation_audit');assert(strcmp(resolved,fullfile(here,'run_evaluation_audit.m')));
code=fileread(resolved);
assert(~contains(code,"fullfile(p09"),'review:UpstreamDependency','Legacy cache dependency remains');
assert(contains(code,"fullfile(review,'inputs','formal')"));
fprintf('Verified %d local input SHA256 values. Entry resolves only local package paths.\n',height(manifest));
end
