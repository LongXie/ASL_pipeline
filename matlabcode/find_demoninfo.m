%% init
clear
clc
cd /Users/long/Desktop

%% read dataset
INCLUDED = dataset('File', 'id_age_sex_includedsubjects.csv', 'Delimiter', ',');
DEMOG = dataset('File', 'Demographics_Combined_10112017.csv', 'Delimiter', ',');

%% select the data
idxout = 1;
OUTDEMOG = DEMOG(1,:);

for i = 1:size(INCLUDED, 1)
    
    % find the index
    id = INCLUDED.ID(i);
    age = INCLUDED.Age(i);
    idx = find(DEMOG.ID == id & DEMOG.Age == age);
    
    % decide where to go
    if size(idx,1) == 0 
        fprintf('Could not find data for %d (age: %d). \n', id, age);
    elseif size(idx,1) > 1
        fprintf('Found multiple input for %d (age: %d). \n', id, age);
    else
        fprintf('Found data for %d (age: %d). \n', id, age);
        
        % combine and get the spreadsheet for the included subjects
        OUTDEMOG(idxout, :) = DEMOG(idx, :);
        idxout = idxout + 1;
        
    end
    
    
end

%% save output
export(OUTDEMOG, 'File', 'Demog_included.csv', 'Delimiter', ',');


