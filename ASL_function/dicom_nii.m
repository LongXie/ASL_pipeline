function [] = dicom_nii(dicomDir,outputDir,outputFile)
%%
%
%   Usage: 
%   dicom_nii(dicomDir,outputDir,outputFile)
%
% Inputs:
%   -dicomDir = path to dicom directory (assumes single series)
%       -if directory has multiple serires, use dicom_sort.m first.
%   -outputDir = path to output directory for nifti file
%   -outputFile = name of output file, e.g. f.nii.gz
%
% Written by Andrew S Bock Feb 2014

fprintf('Converting dicoms to nifti\n');
fprintf(['Input ' dicomDir '\n']);
fprintf(['Output ' fullfile(outputDir,outputFile) '\n']);
dicomList = listdir(dicomDir,'files');
dicom_file = strrep(fullfile(dicomDir,dicomList{1}), ' ', '\ ');
out_file = strrep(fullfile(outputDir,outputFile), ' ', '\ ');

[~, series] = fileparts(dicomDir);

if ~isempty(dicomList)
    if isempty(strfind(series,'asl')) && isempty(strfind(series,'ASL')) 
        [~,~] = system(['mri_convert -odt float ' dicom_file ' ' out_file]);
    elseif  (~isempty(strfind(series,'2D')) || ~isempty(strfind(series,'2d')))
        %2d asl
        [~,~] = system(['mri_convert -odt float ' dicom_file ' ' out_file]);
    else
        %correction for 3D ASL coming out upside down in slice direction
        [~,~] = system(['mri_convert -odt float ' '--in_orientation LPI ' dicom_file ' ' out_file]);
    end
else
    fprintf('Skipping Dir %s - no dicoms found inside.\n',dicomDir);

end