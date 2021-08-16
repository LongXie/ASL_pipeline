function makeimages(h, infile, outfile, range, saveimage, pausetime, maskfile_gz)

%% check input
if nargin < 5
    saveimage = 1;
end

if nargin < 6
    pausetime = 0;
end

if nargin < 7
    maskflag = 0;
else
    maskflag = 1;
end

if isempty(range)
    range_type = 0;
else
    range_type = 1;
end

%% check existance of input file
if ~exist(infile, 'file')
    error('%s does not exist.\n', infile)
else
    [in_dir, infile_filename, in_ext] = fileparts(infile);
    if strcmp(in_ext, '.gz')
        infile_gz = infile;
        infile = fullfile(in_dir, infile_filename);
        system(sprintf('gunzip -c %s > %s', infile_gz, infile));
    elseif strcmp(in_ext, '.nii')
        infile = infile;
    else
        error('Extension %s of input file is not supported.\n', in_ext)
    end
end

%% check existance of mask file
if maskflag == 1
    if ~exist(maskfile_gz, 'file')
        error('%s does not exist.\n', maskfile_gz)
    else
        [mask_dir, mask_filename, mask_ext] = fileparts(maskfile_gz);
        if strcmp(mask_ext, '.gz')
            maskfile = fullfile(mask_dir, mask_filename);
            system(sprintf('gunzip -c %s > %s', maskfile_gz, maskfile));
        elseif strcmp(mask_ext, '.nii')
            maskfile = maskfile_gz;
        else
            error('Extension %s of input file is not supported.\n', mask_ext)
        end
    end
    mask = spm_read_vols(spm_vol(maskfile));
end

%% check and make output dir
[out_dir, out_filename, ext] = fileparts(outfile);
if ~strcmp(ext, '.jpg') || strcmp(ext, '.jepg')
    error('Extension %s of output file is not supported.\n', ext);
end
    
if ~exist(out_dir, 'dir')
    mkdir(out_dir)
end

%% make images
%scrsz=get(0,'Screensize');
%h = figure('Position',[1 1 scrsz(4) scrsz(4)],'Numbertitle','off');

clf(h)

Y = spm_read_vols(spm_vol(infile));
%Y=permute(Y,[1 3 2]);
if ~range_type
    range=[min(Y(:)),max(Y(:))];
end
%figure(h)
if maskflag == 1
    Y(mask==0) = mean(range);
end
M = createmontage(Y);

imshow(M,range,'InitialMagnification','fit')
title(sprintf('%s', out_filename));
colorbar;
drawnow;
if saveimage
    print('-djpeg', '-r0', outfile);
end
if pausetime>0
    pause(pausetime)
else
    pause
end
%close(h)

%% delete unzip file
if strcmp(in_ext, '.gz')
    delete(infile)
end

