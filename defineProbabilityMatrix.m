% Define Probability Matrix

% This function is part of the Statistical Learning of Visual Working
% Memory experiment. This defines the probabilities of the colour pairs in
% each block for each participant

% WN started writing this 1/4/16

% -------------------------------------------------------------------------

% Create a matrix of probabilities for each colour pair.

stimulus.probMatrix = NaN(stimulus.nColours);
stimulus.unifMatrix = NaN(stimulus.nColours);

% Set the diagonal of the matrix to be zero

for thisDiagProb = 1:stimulus.nColours

    stimulus.probMatrix(thisDiagProb,thisDiagProb) = 0;
    stimulus.unifMatrix(thisDiagProb,thisDiagProb) = 0;

end

% Select which colour pairs will be the high-probability pair

randColour = randperm(stimulus.nColours); % Create a random order of 1 to the number of colours
stimulus.highProbPairs_Code = NaN(stimulus.nPairs,2);

for thisPair = 1:stimulus.nPairs

    stimulus.highProbPairs_Code(thisPair,:) = [randColour(thisPair*2-1) randColour(thisPair*2)]; % Creates a list of four pairs of 1 to 8

end

for thisPair = 1:stimulus.nPairs

    for thisColour = 1:2

    stimulus.highProbPairs_Name(thisPair, thisColour) = stimulus.colourList(stimulus.highProbPairs_Code(thisPair,thisColour)); % Creates a list of four pairs of colours from code

    end

end

% Assign high probability to selected colour pairs in the prob matrix

for thisPair = 1:stimulus.nPairs

    thisSelectedPair = stimulus.highProbPairs_Code(thisPair,:);
    stimulus.probMatrix(thisSelectedPair(1),thisSelectedPair(2)) = stimulus.highProb;

end

% Assign low probability to rest of probability matrix

for thisRow = 1:stimulus.nColours

    for thisColumn = 1:stimulus.nColours

        if isnan(stimulus.probMatrix(thisRow,thisColumn))

            stimulus.probMatrix(thisRow,thisColumn) = stimulus.lowProb;

        end

        if isnan(stimulus.unifMatrix(thisRow,thisColumn))

            stimulus.unifMatrix(thisRow,thisColumn) = stimulus.unifProb;

        end

    end

end

% Turn probability matrix into a vector for sampling purposes

stimulus.probVector = [];
stimulus.unifVector = [];

for thisRow = 1:stimulus.nColours
    
    stimulus.probVector = [stimulus.probVector stimulus.probMatrix(thisRow,:)];
    stimulus.unifVector = [stimulus.unifVector stimulus.unifMatrix(thisRow,:)];
    
end

for thisPair = 1:stimulus.nColours^2
    
    if mod(thisPair,8) == 0
        
        stimulus.pairVector(thisPair,:) = [ceil(thisPair/8) 8];
    
    else
        
        stimulus.pairVector(thisPair,:) = [ceil(thisPair/8) mod(thisPair,8)];
            
    end

end

% For convenience, save the high-probability pair numbers

stimulus.highProbPairs = find(stimulus.probVector == stimulus.highProb);

% Save to participant user file
