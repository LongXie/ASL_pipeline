function b = createmontage(I_in)

I=flipdim(permute(I_in,[2 1 4 3]),1);
[nRows,nCols,~,nFrames]=size(I);

aspectRatio = 1;
montageCols = sqrt(aspectRatio * nRows * nFrames / nCols);

% Make sure montage rows and columns are integers. The order in
% the adjustment matters because the montage image is created
% horizontally across columns.
montageCols = ceil(montageCols);
montageRows = ceil(nFrames / montageCols);
montageSize = [montageRows montageCols];

nMontageRows = montageSize(1);
nMontageCols = montageSize(2);
nBands = size(I, 3);

sizeOfBigImage = [nMontageRows*nRows nMontageCols*nCols nBands 1];
rows = 1 : nRows;
cols = 1 : nCols;
k = 1;

for i = 0 : nMontageRows-1
    for j = 0 : nMontageCols-1,
        if k <= nFrames
            b(rows + i * nRows, cols + j * nCols, :) = ...
                I(:,:,:,k);
        else
            return;
        end
        k = k + 1;
    end
end


