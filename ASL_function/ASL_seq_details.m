function ASL_seq_details(dcmDir,outDir)

%   Gets the echo-spacing (msec) and EPI TE (msec) from an EPI dicom file.
%   Requires AFNI (dicom_hdr -sexinfo Siemens EXtra INFO text).
%
%   Outputs the Echo Spacing and EPI TE as text files.
%
%   Written by Andrew S Bock Feb 2014
%   Modified my M Vidorreta Jan 2015

%% Get Dicoms
dicomlist = listdir(dcmDir,'files');

dictionary = 'Mydicom-dict.txt';
dcmfilename = fullfile(dcmDir,dicomlist{end});
dcmhdr = ReadSiemensExtraInfoHeaderQuickFix( dcmfilename, dictionary );

%% Find TE
TE = dcmhdr.EchoTime;

%% Find if it's 3D sequence
if isempty(strfind(dcmhdr.tSequenceFileName, '3D'))
    is3D = 0;
    if (isempty(strfind(dcmhdr.tSequenceFileName, 'PASL')) && isempty(strfind(dcmhdr.tSequenceFileName, 'FAIR')) && isempty(strfind(dcmhdr.tSequenceFileName, 'pasl')) && isempty(strfind(dcmhdr.tSequenceFileName, 'fair')))
        isPCASL = 1;
    end
else
    is3D = 1;
    isPCASL = 1;
end



%% Save echo spacing and EPI TE as text files
system(['echo ' num2str(TE) ' > ' fullfile(outDir,'TE')]);
system(['echo ' num2str(is3D) ' > ' fullfile(outDir,'is3D')]);
system(['echo ' num2str(isPCASL) ' > ' fullfile(outDir,'isPCASL')]);

end
