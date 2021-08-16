function sort_nifti(session_dir,dicom_dir,series,outFN)

%   Sorts dicoms into series directories, converts to nifti files based on
%   series type (e.g. MPRAGE, BOLD, DTI). Also ACPC aligns anatomical
%   images. Assumes dicoms are found in a "DICOMS" directory, in the
%   session directory.
%
%   Usage: sort_nifti(session_dir,dicom_dir)
%
%   note: for the ACPC step to work, you need to have python installed.
%   Also, by default 'ACPC.m' uses an FSL anatomical atlas, so it is
%   recommended that FSL also be installed.
%
%   note: for the B0 step to work, be sure that AFNI is installed, and that
%   you have the line below in your startup file:
%
%   setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
%
%   Written by Andrew S Bock Feb 2014
%
% 12/15/14      spitschan       Output to be run in the shell also saved in
%                               file.
%
%   Modified my M Vidorreta Jan 2015

%% Add to log
%SaveLogInfo(session_dir, mfilename) %MV error in SaveLogInfo - GetGITInfo
%function not found

%% Set default parameters
%if ~exist('dicom_dir','var')
%    dicom_dir = fullfile(session_dir,'DICOMS');
%end
%% sort dicoms within this directory
%dicom_sort(dicom_dir);
%series = listdir(dicom_dir,'dirs');
%series
% process series types
if ~isempty(series)
    mpragect = 0;
    mp2ragect = 0;
    PDct = 0;
    boldct = 0;
    aslct = 0;
    DTIct = 0;
    %progBar = ProgressBar(length(series),'Converting dicoms'); %function not found
    for s = 1:length(series)
        fprintf(['\nProcessing ' series{s} ' series ' num2str(s) ' of ' ...
            num2str(length(series)) '\n'])
        % Anatomical image
        if ( ~isempty(strfind(series{s},'MPRAGE')) && isempty(strfind(series{s},'MPRAGE_NAV')) );
            
            mpragect = mpragect + 1;
            fprintf(['PROCESSING ANATOMICAL IMAGE ' num2str(mpragect) '\n']);
            
            % Convert dicoms to nifti
            outputDir = fullfile(session_dir,outFN{s}); 
            mkdir(outputDir);
            outFile = 'MPRAGE.nii.gz';
            dicom_nii(fullfile(dicom_dir,series{s}),outputDir,outFile);
            
%             betname = 'MPRAGE';
%             system(['bet ' fullfile(outputDir,outFile) ' ' ...
%                 fullfile(outputDir,[betname '_brain.nii.gz']) ' -f 0.1']);
%             %ACPC(outputDir,betname);

            system(['echo ' series{s} ' > ' fullfile(outputDir,'series_name')]);
            disp('done.')
            
        elseif ~isempty(strfind(series{s},'mp2rage'));
            mp2ragect = mp2ragect + 1;
            fprintf(['\nPROCESSING ANATOMICAL IMAGE ' num2str(mp2ragect) '\n']);
            % Convert dicoms to nifti
            outputDir = fullfile(session_dir,'MP2RAGE',['00' num2str(mp2ragect)]);
            mkdir(outputDir);
            outFile = 'MP2RAGE.nii.gz';
            dicom_nii(fullfile(dicom_dir,series{s}),outputDir,outFile)
            betname = 'MP2RAGE';
            system(['bet ' fullfile(outputDir,outFile) ' ' ...
                fullfile(outputDir,[betname '_brain.nii.gz']) ' -f 0.1']);
            ACPC(outputDir,betname);
            system(['echo ' series{s} ' > ' fullfile(outputDir,'series_name')]);
            disp('done.')
        elseif ~isempty(strfind(series{s},'PD'));
            PDct = PDct + 1;
            fprintf(['\nPROCESSING PROTON DENSITY IMAGE ' num2str(PDct) '\n']);
            % Convert dicoms to nifti
            outputDir = fullfile(session_dir,'PD',['00' num2str(PDct)]);
            mkdir(outputDir);
            outFile = 'PD.nii.gz';
            dicom_nii(fullfile(dicom_dir,series{s}),outputDir,outFile)
            system(['echo ' series{s} ' > ' fullfile(outputDir,'series_name')]);
            disp('done.')
%         elseif ~isempty(strfind(series{s},'ep2d')) || ~isempty(strfind(series{s},'BOLD')) ...
%                 || ~isempty(strfind(series{s},'bold')) || ~isempty(strfind(series{s},'EPI')) ...
%                 || ~isempty(strfind(series{s},'RUN'));
        elseif ~isempty(strfind(series{s},'BOLD')) ...
                || ~isempty(strfind(series{s},'bold')) ;
            boldct = boldct + 1;
            fprintf(['\nPROCESSING BOLD IMAGE ' num2str(boldct) '\n']);
            % Convert dicoms to nifti
            %outputDir = fullfile(session_dir,series{s});
            outputDir = fullfile(session_dir, sprintf('BOLD_%03d',boldct));
            mkdir(outputDir);
            outFile = 'raw_f.nii.gz';
            dicom_nii(fullfile(dicom_dir,series{s}),outputDir,outFile)
            system(['echo ' series{s} ' > ' fullfile(outputDir,'series_name')]);
            echo_spacing(fullfile(dicom_dir,series{s}),outputDir);
            slice_timing(fullfile(dicom_dir,series{s}),outputDir);
            disp('done.')       

        elseif ~isempty(strfind(series{s},'PCASL')) ...
                || ~isempty(strfind(series{s},'pcasl')) ;
            aslct = aslct + 1;
            fprintf(['PROCESSING ASL IMAGE ' num2str(aslct) '\n']);
            
            % Convert dicoms to nifti
            %outputDir = fullfile(session_dir,series{s});
            if isempty(strfind(series{s},'M0'))
                outputDir = fullfile(session_dir, outFN{s});
            else
                outputDir = fullfile(session_dir, [outFN{s}, '_M0']);
            end
            mkdir(outputDir);
            outFile = 'RAWPCASL.nii.gz';
            dicom_nii(fullfile(dicom_dir,series{s}),outputDir,outFile)
            
            system(['echo ' series{s} ' > ' fullfile(outputDir,'series_name')]);
            ASL_seq_details(fullfile(dicom_dir,series{s}),outputDir);
            disp('done.')   
        elseif ~isempty(strfind(series{s},'DTI'));
            DTIct = DTIct + 1;
            disp('done.')
        end
        %progBar(s);
    end
end
% fprintf('\n');
% disp('Preprocess finished');
% fprintf('\n');
% disp('to run recon-all, run the following in a terminal (assuming using 3T and MPRAGE):')
% fprintf('\n');
% commandc = ['recon-all -i ' fullfile(session_dir,'MPRAGE','001','ACPC','MPRAGE.ACPC.nii.gz') ' -s <subject_name> -all'];
% disp(commandc);

% Also save this out in a file.
%fid = fopen(fullfile(session_dir, 'recon_all_scripts'), 'a');
%fprintf(fid, commandc);
%fclose(fid);

end
