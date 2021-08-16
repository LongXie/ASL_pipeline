function Convert3Dto4D_session(session_dir)

%% Find bold run directories
d = listdir(fullfile(session_dir,'*ASL*'),'dirs');
if isempty(d) %MV
    d = listdir(fullfile(session_dir,'*asl*'),'dirs');
end
nruns = length(d);

if nruns == 0
    fprintf('No ASL directories found in %s.\n',session_dir);
    return;
%else
%    fprintf('The following ASL directories were found in %s:\n',session_dir);
%    d
end


%% Process each run
for r = 1:nruns    
    fprintf('....Convert 3D data to 4D (nifti) for run %s:\n',d{r});
        
    % convert 3D data to 4D (nifti)
    files = spm_select('FPlist', fullfile(session_dir,d{r}), '^raw\w*.*img');
    
    if isempty(files)
        fprintf('........Run %s has been processed.\n',d{r})
    else
        spm_file_merge(cellstr(files), 'raw.nii.gz', 4);
        delete(fullfile(session_dir, d{r}, 'raw*.img'));
        delete(fullfile(session_dir, d{r}, 'raw*.hdr'));
    end
    
end

