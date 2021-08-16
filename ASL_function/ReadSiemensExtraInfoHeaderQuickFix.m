%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Siemens CSA Header info from dicom header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% dcmfilename = fullfile('C:', '20141216', 'dicom', '001_000005_000001.dcm');
% matdocpath = userpath;
% dictionary  = fullfile(matdocpath(1:end-1), 'Mydicom-dict.txt');
%
% OUTPUT:
% dcmhdr (struct with multiple fields)
% Look for:
% dcmhdr.ucReadOutMode
% dcmhdr.tSequenceFileName
% dcmhdr.TProtocolName
% dcmhdr.TE_0
% dcmhdr.TE_1
% ...
% other product attributes are:
% dcmhdr.RepetitionTime
% dcmhdr.EchoTime
% dcmhdr.SeriesNumber
% ...

function dcmhdr = ReadSiemensExtraInfoHeaderQuickFix( dcmfilename, dictionary )

    dcmhdr = dicominfo(dcmfilename, 'dictionary', dictionary);
    
    %quick fix
    fid = fopen(dcmfilename);
    A=fread(fid,'char');
    fclose(fid);
    Achar = char(A');
    idx1=strfind(Achar,'ASCCONV BEGIN');
    idx2=strfind(Achar,'ASCCONV END');
   
    SiemensHeader = Achar(idx1(1):idx2(end));
        
%     SiemensHeader = char(dcmhdr.SiemensCSAHeader');
    dcmhdr.SiemensCSAHeader = SiemensHeader;
    
    %tags.attributes = {'ucReadOutMode', 'tSequenceFileName', 'tProtocolName','alTE[0]','alTE[1]','sPat.lAccelFactPE','sSliceArray.ucMode','sSliceArray.lSize'};
    %tags.outnames   = {'ucReadOutMode', 'tSequenceFileName', 'tProtocolName','TE_0','TE_1','AccelFactPE','AcqDirMode','NSlices'};
    
    tags.attributes = {'ucReadOutMode', 'tSequenceFileName', 'tProtocolName','sPat.lAccelFactPE','sSliceArray.ucMode','sSliceArray.lSize'};
    tags.outnames   = {'ucReadOutMode', 'tSequenceFileName', 'tProtocolName','AccelFactPE','AcqDirMode','NSlices'};
    
    
    n_tags = length(tags.attributes);
    
    tags.values     = cell(1, n_tags);
    
    for i = 1:n_tags
        idx = strfind(SiemensHeader, tags.attributes{i});
        if isempty(idx)
            fprintf('==== WARNING: Tag %d = %s not found in dicom header. ====\n',i, tags.attributes{i});
        else
            remain = SiemensHeader( (idx + length(tags.attributes{i}):end));
            [token, remain] = strtok(remain, [' =["]', char(32),char(10),char(13)]); 
    
            %Just to polish output xd
            if ~strcmp(tags.attributes{i},'tSequenceFileName')
                tags.values{i} = token;
            else
                token2   = token;
                remain2  = token;
                finished = 0;
                while (finished==0)
                    [token2,remain2]=strtok(remain2,'%\/');
                    if isempty(token2)
                        finished=1;
                    elseif isempty(remain2)
                        finished=2;
                    end
                end
                if finished==1
                    tags.values{i} = remain2;
                else
                    tags.values{i} = token2;
                end
            end
            
            if ~isempty(strfind(tags.values{i},'+AF8-'))
                idx2 = strfind(tags.values{i},'+AF8-');
                tags.values{i} = [tags.values{i}(1:(idx2-1)) '_' tags.values{i}(idx2+5:length(tags.values{i}))];
            end
        end
        
        eval(sprintf('dcmhdr.%s = ''%s'';\n',tags.outnames{i},tags.values{i}));
    end        

%     if isempty(strfind(dcmhdr.tProtocolName,'B0'))
%         fid=fopen(dcmfilename);
%         A=fread(fid,'char');
%         Achar = char(A');
%         idx=strfind(Achar,'MosaicRefAcqTimes');
%         remain = Achar( (idx(1) + length('MosaicRefAcqTimes'):end));
%         startstring = char([77    0    0    0   11]); %This corresponds to 'M   ' string preceding the slice timings and ensures nothing in between is read by mistake
%         idx = strfind(remain,startstring);
%         remain = remain(idx(1):end);
%         SliceTiming = zeros(str2double(dcmhdr.NSlices),1);
%         islice = 1;
%         while islice <= str2double(dcmhdr.NSlices)     
%             [token,remain] = strtok(remain,'0123456789.');
%             idx=1;
%             while (int8(remain(idx)) > 45 && int8(remain(idx)) < 58 && int8(remain(idx)) ~= 47) %ASCII indeces for 0123456789 are 48-57 and . is 46
%                 idx = idx + 1;
%             end
%             currtiming = str2double(remain(1:idx-1));
%             if currtiming > dcmhdr.RepetitionTime
%                 ; %skip this number, is not a slice timing index
%             else
%                 SliceTiming(islice,1) = str2double(remain(1:idx-1));
%                 islice = islice + 1;
%             end                
%         end    
%         fclose(fid);
%            
%         dcmhdr.SliceTimings = SliceTiming;
%     end
    

end