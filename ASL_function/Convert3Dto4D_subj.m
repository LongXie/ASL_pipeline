function Convert3Dto4D_subj(subj_dir)

%% Find bold run directories
d = listdir(fullfile(subj_dir,'*'),'dirs');
ntps = length(d);

if ntps == 0
    fprintf('No time point directories found in %s.\n',subj_dir);
    return;
end


%% Process each time point
for r = 1:ntps    
    fprintf('..Convert 3D data to 4D (nifti) for time point %s:\n', d{r});

    % process each time point
    session_dir = fullfile(subj_dir, d{r});
    Convert3Dto4D_session(session_dir);

end