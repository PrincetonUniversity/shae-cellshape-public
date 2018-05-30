function outputArray = normalizeRows(inputArray);

try
    if issparse(inputArray)
        inputArray = full(inputArray);
    end
    
    rowVals = sqrt(nansum(inputArray.^2,2));
    outputArray = inputArray./repmat(rowVals,[1,size(inputArray,2),1,1,1,1,1,1]);
    
catch ME % use the builtin function (part of neural net toolbox)
   outputArray = normr(inputArray); 
end