function buildProgressBarString(inputMat)
% takes inputs similar to what should go into progressbar, 
% but as a matrix instead of individual input args

nDimBar = numel(inputMat);

% build string for eval
evalString = 'progressbar(';
for iiDimBar = 1:nDimBar
    evalString = [evalString,...
        num2str(inputMat(iiDimBar)),' , '];
end

evalString(end-1:end) = [];
evalString = [evalString,');'];
eval(evalString);