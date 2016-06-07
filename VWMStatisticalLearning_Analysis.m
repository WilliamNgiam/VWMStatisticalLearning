% VWM Statistical Learning Analysis
%
% WXQN started writing this on 13th June 2016
% This code is for data analysis of VWMStatisticalLearning.m
%
% Currently this code:
%
% * Compiles the number of correct trials in each block for each participant
% * Calculates the percentage correct and VWM capacity (K) in each block
% * Plots percentage correct and K over blocks for each participant
% * Counts the number of trials including low probability pairs in each
%   pattern condition block for each participant

% ----------------------------------------------------------------------- %

cd('/Users/alexh/Documents/MATLAB/William/VWMStatisticalLearning/Data/');

theseFiles = what;
theseFiles = theseFiles.mat;
nFiles = numel(theseFiles);
load(theseFiles{1});
participantID = [];

for thisFile = 1:nFiles

load(theseFiles{thisFile});
thisID = participant.ID;

participantID = [participantID; thisID];

end

participantID = unique(cellstr(participantID))';
numParticipants = numel(participantID);
    
% Compile responses
nCorr = NaN(numParticipants,20);
pCorr = NaN(numParticipants,20);   
K = NaN(numParticipants,20);

nCorr_unif = NaN(numParticipants,10);
pCorr_unif = NaN(numParticipants,10);
K_unif = NaN(numParticipants,10);

nCorr_patt = NaN(numParticipants,10);
pCorr_patt = NaN(numParticipants,10);
K_patt = NaN(numParticipants,10);

for thisNum = 1:numParticipants

    for thisFile = 1:nFiles

        % Load the file
        load(theseFiles{thisFile});
        thisParticipantID = participantID{thisNum};
        thisID = participant.ID;
        thisBlock = block.thisBlock; 
        
        if strcmp(thisID,thisParticipantID)
        
            nCorr(thisNum,thisBlock) = sum(block.allCorrect);
            pCorr(thisNum,thisBlock) = sum(block.allCorrect)/experiment.nTrialsPerBlock;
            K(thisNum,thisBlock) = (64*pCorr(thisNum,thisBlock)-8)/7;
        end
        
    end
     
    disp(nCorr);
   
    thisParticipant = participantID{thisNum};
    firstFileName = [thisParticipant '_VWMStatLearning_1.mat'];
    load(firstFileName);
    
    if block.thisCondition == 1  % Uniform first
        
        for thisBlock = 1:experiment.nBlocks
            
            if mod(thisBlock,2) == 1
                
                thisConditionBlock = (thisBlock+1)/2;
                nCorr_unif(thisNum,thisConditionBlock) = nCorr(thisNum,thisBlock);
                pCorr_unif(thisNum,thisConditionBlock) = pCorr(thisNum,thisBlock);
                K_unif(thisNum,thisConditionBlock) = K(thisNum,thisBlock);
                
            elseif mod(thisBlock,2) == 0
                
                thisConditionBock = thisBlock/2;
                nCorr_patt(thisNum,thisConditionBlock) = nCorr(thisNum,thisBlock);
                pCorr_patt(thisNum,thisConditionBlock) = pCorr(thisNum,thisBlock);
                K_patt(thisNum,thisConditionBlock) = K(thisNum,thisBlock);
                
            end
            
        end
        
    elseif block.thisCondition == 2 % Patterned first
                       
        for thisBlock = 1:experiment.nBlocks
            
            if mod(thisBlock,2) == 0
                
                thisConditionBlock = (thisBlock)/2;
                nCorr_unif(thisNum,thisConditionBlock) = nCorr(thisNum,thisBlock);
                pCorr_unif(thisNum,thisConditionBlock) = pCorr(thisNum,thisBlock);
                K_unif(thisNum,thisConditionBlock) = K(thisNum,thisBlock);
                
            elseif mod(thisBlock,2) == 1
                
                thisConditionBlock = (thisBlock+1)/2;
                nCorr_patt(thisNum,thisConditionBlock) = nCorr(thisNum,thisBlock);
                pCorr_patt(thisNum,thisConditionBlock) = pCorr(thisNum,thisBlock);
                K_patt(thisNum,thisConditionBlock) = K(thisNum,thisBlock);
                
            end
            
        end

    end
    
    Conditions = {'Uniform', 'Pattern'};
    figure('Color','white','Name',thisParticipantID);
    subplot(2,1,1);
    plot(nCorr_unif(thisNum,:),'bx-');
    xlabel('Block Number');
    ylabel('Number of Correct Trials');
    set(gca,'TickDir', 'out')
    set(gca,'YMinorTick', 'on')
    hold on;
    plot(nCorr_patt(thisNum,:),'rx-');
    legend(Conditions,'Location','NorthOutside');
    
    subplot(2,1,2);
    plot(K_unif(thisNum,:),'bx-');
    xlabel('Block Number');
    ylabel('Capacity (K)');
    set(gca,'TickDir', 'out')
    set(gca,'YMinorTick', 'on')
    hold on;
    plot(K_patt(thisNum,:),'rx-');
    hold off;
    legend(Conditions,'Location','NorthOutside');
       
end

% Combined figure of all participants (better for low n)

figure('Color','white','Name','Mean Data');

for thisParticipant = 1:numParticipants
    
    subplot(numParticipants,2,thisParticipant*2-1);
    plot(nCorr_unif(thisParticipant,:),'bx-');
    xlabel('Block Number');
    ylabel('Number of Correct Trials');
    set(gca,'TickDir', 'out')
    set(gca,'YMinorTick', 'on')
    hold on;
    plot(nCorr_patt(thisParticipant,:),'rx-');
    
%     if thisParticipant == 1
%     
%         legend(Conditions,'Location','NorthOutside');
%         
%     end
    
    subplot(numParticipants,2,thisParticipant*2);
    plot(K_unif(thisParticipant,:),'bx-');
    xlabel('Block Number');
    ylabel('Capacity (K)');
    set(gca,'TickDir', 'out')
    set(gca,'YMinorTick', 'on')
    hold on;
    plot(K_patt(thisParticipant,:),'rx-');
    hold off;
%     
%     if thisParticipant == 1
%         
%         legend(Conditions,'Location','NorthOutside');
%         
%     end
    
end

% Figure of mean data

figure;
subplot(2,1,1);
plot(mean(nCorr_unif),'bx-');
xlabel('Block Number');
ylabel('Mean Number of Correct Trials');
set(gca,'TickDir','out');
set(gca,'YMinorTick','on');
hold on;
plot(mean(nCorr_patt),'rx-');
legend(Conditions,'Location','NorthOutside');

subplot(2,1,2);
plot(mean(K_unif),'bx-');
xlabel('Block Number');
ylabel('Mean K');
set(gca,'TickDir','out');
set(gca,'YMinorTick','on');
hold on;
plot(mean(K_patt),'rx-');
legend(Conditions,'Location','NorthOutside');

% Checking the number of trials of low probability pairs and number of
% trials that contained only high probability pairs in each block for each
% participant.

data.allLowProbPairsNum = NaN(numParticipants,9);
data.allHighProbPairTrialsNum = NaN(numParticipants,9);

for thisParticipantNum = 1:numParticipants
        
    allLowProbPairsNum = [];
    allHighProbPairTrialsNum = [];
    
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
        
        if strcmp(thisID,thisParticipantID)
        
            if thisCondition == 2

                for thisPair = 1:stimulus.nPairs

                    thisOne = stimulus.highProbPairs(thisPair);

                    a = a + (block.allPairs == thisOne);

                end
    %           disp('a');
    %           disp(a);

                disp('nLowProbPairs');
                nLowProbPairs = numel(find(a == 0));
                disp(nLowProbPairs);
                allLowProbPairsNum = [allLowProbPairsNum nLowProbPairs];

                disp('nHighProbPairTrials');
                nHighProbPairTrials = sum(sum(a')==4);
                disp(nHighProbPairTrials);
                allHighProbPairTrialsNum = [allHighProbPairTrialsNum nHighProbPairTrials];

            end

        end
        
    end
    
    data.allLowProbPairsNum(thisParticipantNum,:) = allLowProbPairsNum;
    data.allHighProbPairTrialsNum(thisParticipantNum,:) = allHighProbPairTrialsNum;
    
end

% Testing awareness

cd('/Users/alexh/Documents/MATLAB/William/VWMStatisticalLearning/UserData/');

theseFiles = what;
theseFiles = theseFiles.mat;
nFiles = numel(theseFiles);
load(theseFiles{1});
participantID = [];

for thisFile = 1:nFiles

load(theseFiles{thisFile});
thisID = participant.ID;

participantID = [participantID; thisID];

end

participantID = unique(cellstr(participantID))';
numParticipants = numel(participantID);

for thisParticipantNum = 1:numParticipants
    
    awareness.ifCorrect = [];
    
    % Retrieve the high probability pairs
    
    load(theseFiles{thisParticipantNum});
    
    % Retrieve their awareness test responses
    
    for thisColour = 1:stimulus.nColours
    
        theirResponse = awareness.responseColours(thisColour);          % Retrieve their response to what was paired to this Colour
        thisPosition = find(stimulus.highProbPairs_Code == thisColour); % Retrieve the position of the colour in the 4x2 matrix
        
        if thisPosition > 4
            
            thisColourPairPosition = thisPosition - 4;
            
        elseif thisPosition < 5
            
            thisColourPairPosition = thisPosition + 4;
            
        end
        
        theColourPair = stimulus.highProbPairs_Code(thisColourPairPosition);
        
        if theirResponse == theColourPair       % Correctly identified the pair
            
            awareness.ifCorrect = [awareness.ifCorrect 1];
            
        elseif theirResponse ~= theColourPair   % Did not identify the pair
            
            awareness.ifCorrect = [awareness.ifCorrect 0];
            
        end
        
    end
    
    % Save and attach awareness to participant
    
    userFileName = [participant.ID '_VWMStatLearning_Exp1.mat'];
    save(userFileName, 'participant', 'experiment', 'equipment', 'colour', 'stimulus', 'timing', 'awareness');

end        
        
            