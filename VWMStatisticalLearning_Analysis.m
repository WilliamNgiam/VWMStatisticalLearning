% VWM Statistical Learning Analysis

cd('/Users/alexh/Documents/MATLAB/William/VWMStatisticalLearning/Data/');

theseFiles = what;
theseFiles = theseFiles.mat;
nFiles = numel(theseFiles);
load(theseFiles{1});

participantID = {'02'};
numParticipants = numel(participantID);

for thisNum = 1:numParticipants
    
     % Compile responses
    nCorr = NaN(1,20);
    pCorr = NaN(1,20);   
    
    for thisFile = 1:nFiles

        % Load the file
        load(theseFiles{thisFile});
        participantID = {'02'};
        thisParticipantID = participantID{thisNum};
        thisID = participant.ID;
        thisBlock = block.thisBlock; 
        
        if strcmp(thisID,thisParticipantID)
        
            nCorr(thisBlock) = sum(block.allCorrect);
            pCorr(thisBlock) = sum(block.allCorrect)/experiment.nTrialsPerBlock;

        end
        
    end
     
    disp(nCorr);
    nCorr_unif = NaN(1,10);
    pCorr_unif = NaN(1,10);
    nCorr_patt = NaN(1,10);
    pCorr_patt = NaN(1,10);

    thisParticipant = participantID{thisNum};
    firstFileName = [thisParticipant '_VWMStatLearning_1.mat'];
    load(firstFileName);
    
    if block.thisCondition == 1  % Uniform first
        
        for thisBlock = 1:experiment.nBlocks
            
            if mod(thisBlock,2) == 1
                
                thisConditionBlock = (thisBlock+1)/2;
                nCorr_unif(thisConditionBlock) = nCorr(thisBlock);
                pCorr_unif(thisConditionBlock) = pCorr(thisBlock);
                
            elseif mod(thisBlock,2) == 0
                
                thisConditionBock = thisBlock/2;
                nCorr_patt(thisConditionBlock) = nCorr(thisBlock);
                pCorr_patt(thisConditionBlock) = pCorr(thisBlock);
                
            end
            
        end
        
    elseif block.thisCondition == 2 % Patterned first
                       
        for thisBlock = 1:experiment.nBlocks
            
            if mod(thisBlock,2) == 0
                
                thisConditionBlock = (thisBlock)/2;
                nCorr_unif(thisConditionBlock) = nCorr(thisBlock);
                pCorr_unif(thisConditionBlock) = pCorr(thisBlock);
                
            elseif mod(thisBlock,2) == 1
                
                thisConditionBlock = (thisBlock+1)/2;
                nCorr_patt(thisConditionBlock) = nCorr(thisBlock);
                pCorr_patt(thisConditionBlock) = pCorr(thisBlock);
                
            end
            
        end

    end
    
    figure('Color','white','Name',thisParticipantID);
    subplot(2,1,1);
    plot(nCorr_unif,'bx-');
    hold on;
    plot(nCorr_patt,'rx-');

    subplot(2,1,2);
    plot(pCorr_unif,'bx-');
    hold on;
    plot(pCorr_patt,'rx-');
    hold off;
    
end

cd('/Users/alexh/Documents/MATLAB/William/VWMStatisticalLearning/PilotData/');

theseFiles = what;
theseFiles = theseFiles.mat;
nFiles = numel(theseFiles);
load(theseFiles{1});

participantID = {'02'};
numParticipants = numel(participantID);

for thisParticipantNum = 1:numParticipants
        
    data.allLowProbPairsNum = [];
    data.allHighProbPairTrialsNum = [];
    
    for thisFile = 1:nFiles

        % Load the file
        load(theseFiles{thisFile});
        thisParticipantID = participantID{thisParticipantNum};
        thisID = participant.ID;
        thisBlock = block.thisBlock;
        thisCondition = block.thisCondition;
        
        a = zeros(experiment.nTrialsPerBlock,stimulus.nPairs);
        
        if exist('stimulus.highProbPairs') == 0

            stimulus.highProbPairs = find(stimulus.probVector == stimulus.highProb);

        end
        
        if thisCondition == 2
        
            for thisPair = 1:stimulus.nPairs

                thisOne = stimulus.highProbPairs(thisPair);

                a = a + (block.allPairs == thisOne);

            end
            disp('a');
            disp(a)

            disp('nLowProbPairs');
            nLowProbPairs = numel(find(a == 0));
            disp(nLowProbPairs);
            data.allLowProbPairsNum = [data.allLowProbPairsNum nLowProbPairs];

            disp('nHighProbPairTrials');
            nHighProbPairTrials = sum(sum(a')==4);
            disp(nHighProbPairTrials);
            data.allHighProbPairTrialsNum = [data.allHighProbPairTrialsNum nHighProbPairTrials];

        end

    end
    
end
