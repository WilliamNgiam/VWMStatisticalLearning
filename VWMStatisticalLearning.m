% Statistical Learning in Visual Working Memory

% This experiment is designed to examine the influence of statistical
% learning on visual working memory capacity. This experiment is a
% replication of the retracted Anderson paper

% WN started writing this on 1/4/16

% -------------------------------------------------------------------------

clear all;
Screen('CloseAll');

% Set up experiment parameters

experiment.nBlocks = 20;            % Number of blocks
experiment.nTrialsPerBlock = 60;    % Number of trials per block

% Set up participant parameters

userDirectory = '/Users/wngi5916/Documents/MATLAB/VWMStatisticalLearning/UserData/';
saveDirectory = '/Users/wngi5916/Documents/MATLAB/VWMStatisticalLearning/Data/';
experimentDirectory = '/Users/wngi5916/Documents/MATLAB/VWMStatisticalLearning/';    

% Set up equipment parameters

equipment.viewDist = 770;           % Viewing distance in mm
equipment.ppm = 2.7;                % Pixels per mm
equipment.gammaVals = 1.0./[3.0902 2.4049 2.3194];      % Gamma values for CRT in GT519 (recalibrated 25/9/15)

equipment.greyVal = .5;             % Grey background
equipment.blackVal = 0;
equipment.whiteVal = 1;

equipment.responseKeys = 30:37;     % Responses key for 1 to 8 on Macbook Pro keyboard

% Set up colour parameters

colour.yellow = [255, 255, 0];
colour.blue = [0, 0, 255];
colour.red = [255, 0, 0];
colour.green = [0, 255, 0];
colour.magenta = [255, 0, 255];
colour.cyan = [0, 255, 255];
colour.white = [255, 255, 255];
colour.black = [0, 0, 0];

colour.textVal = 0;
colour.fixVal = 1;
colour.probeVal = 0;

% Set up stimulus parameters

    % Colours
    
stimulus.nColours = 8;
stimulus.nPairs = stimulus.nColours/2;
stimulus.colourList = {'Yellow', 'Blue', 'Red', 'Green', 'Magenta', 'Cyan', 'White', 'Black'};
stimulus.colours = [colour.yellow; colour.blue; colour.red; colour.green; colour.magenta; colour.cyan; colour.white; colour.black];

stimulus.size_dva = .9;                     % Colour size (diameter of circle, width of square)
stimulus.fixationEccentricity_dva = 3;    % Eccentricity of stimulus from fixation point
stimulus.pairEccentricity_dva = 1;          % Eccentricity between (center of) colour pairs
stimulus.refEccentricity_dva = 1.5;         % Eccentricity between the colour references

    % Probabilities

stimulus.highProb = 80/372;
stimulus.lowProb = 1/372;
stimulus.unifProb = 1/(stimulus.nColours^2-stimulus.nColours);  % 1/56 if stimulus.nColours = 8

    % Fixation
    
stimulus.fixationOn = 1;            % Fixation on (1) or off (0)
stimulus.fixationSize_dva = .3;     % Fixation size in degrees of visual angle

    % Probe Rectangle
    
stimulus.thinPenWidth = 1;          % Thin line for memoranda rect that aren't the probe
stimulus.thickPenWidth = 5;         % Thick line for the memoranda rect that is the probe

% Set up temporal parameters

timing.ITI = .75;       % Inter-trial interval
timing.memory = 1;
timing.delay = 1;

% Set up Psychtoolbox Pipeline

AssertOpenGL;

    % Imaging set-up

screenID = max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
Screen('Preference','SkipSyncTests',2);

    % Window set-up
    
[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.greyVal);
PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);

    % Text set-up
    
Screen('TextFont',ptbWindow,'Courier New');
Screen('TextSize',ptbWindow,16);
Screen('TextStyle',ptbWindow,1);        % Bold text

% Enable alpha blending for typical drawing of masked textures
Screen('BlendFunction', ptbWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Calculate equipment parameters

equipment.mpd = (equipment.viewDist)*tan(deg2rad(2*stimulus.fixationEccentricity_dva))/stimulus.fixationEccentricity_dva; % Calculate mm per degree of visual angle to the ecccentricity of the stimuli
equipment.ppd = equipment.ppm*equipment.mpd;        % Pixels per degree

% Calculate spatial parameters

stimulus.size_pix = round(stimulus.size_dva*equipment.ppd);                                     % Item size in pixels
stimulus.fixationEccentricity_pix = round(stimulus.fixationEccentricity_dva*equipment.ppd);     % Eccentricity of stimulus in pixels
stimulus.fixationSize_pix = stimulus.fixationSize_dva*equipment.ppd;                            % Fixation dot size in pixels
stimulus.pairEccentricity_pix = round(stimulus.pairEccentricity_dva*equipment.ppd);             % Eccentricity between colour pairs in pixels
stimulus.refEccentricity_pix = round(stimulus.refEccentricity_dva*equipment.ppd);               % Eccentricity between colour references in response screen

% Get participant ID

while true
    participant.ID = upper(input('Enter two-digit participant number: ', 's'));
    if length(participant.ID)==2
        break
    end
end

cd(userDirectory);
participant.userFile = [userDirectory participant.ID '_VWMStatLearning_Exp1.mat'];
newUser = ~exist(participant.userFile,'file');

if newUser
    
    % Create participant parameters
    
        % Run function to designate probabilities to colour pairs for each
        % block. Remember to change to uniform probabilities for last block
    
        defineProbabilityMatrix;
       
        % Designate conditions for each block
        
            % Determine if rectangles first or circles first
            
            whichShapeFirst = mod(str2num(participant.ID),2)+1;
        
            if whichShapeFirst == 1
                
                participant.shapeOrder = repmat([1,2],1,10);
            
            elseif whichShapeFirst == 2
                
                participant.shapeOrder = repmat([2,1],1,10);
                
            end
            
            % Determine if pattern first or uniform first
        
            whichCondnFirst = ceil((mod(str2num(participant.ID),4)+1)/2);
            
            if whichCondnFirst == 1
                
                participant.condnOrder = repmat([1,2],1,9);
                participant.condnOrder = [participant.condnOrder 1 1];  % Make the final two blocks uniform conditions
                
            elseif whichCondnFirst == 2
                
                participant.condnOrder = repmat([2,1],1,9);
                participant.condnOrder = [participant.condnOrder 1 1];  % Make the final two blocks uniform conditions
                
            end     
        
    participant.thisBlock = 1;
    save(participant.userFile, 'participant');
    
else
    
    load(participant.userFile);
    participant.thisBlock = participant.thisBlock+1;
    
end

% Instruction Text



% Set up location rects

fixRect = [0 0 stimulus.fixationSize_pix stimulus.fixationSize_pix];    % Fixation rect
fixRect = CenterRectOnPoint(fixRect, screenCentreX, screenCentreY);     % Centred in the middle of the screen

colourRect = [0 0 stimulus.size_pix stimulus.size_pix];                 % Colour stimulus rect

% Build location rects of the colours going clockwise starting from north

colourRects = NaN(4,stimulus.nColours);

colourTheta_deg = linspace(0,270,stimulus.nPairs)-90;       %[-90 0 90 180] - Rotation of colour pairs in degrees (North, East, South, West)
colourTheta_rad = deg2rad(colourTheta_deg)';                % Rotation of the colour pairs in radians (North, East, South, West)
[theseColourX, theseColourY] = pol2cart(colourTheta_rad,stimulus.fixationEccentricity_pix*ones(stimulus.nPairs,1));

for thisPair = 1:stimulus.nPairs
    
    for thisColour = 1:2        % Create rect for each colour block in pair, separate them vertically by adding eccentricity to y-value of rect
        
        if thisColour == 1      % Subtracting half the eccentricity
            
            colourRects(:,thisPair*2-1) = CenterRectOnPoint(colourRect, theseColourX(thisPair)+screenCentreX, theseColourY(thisPair)+screenCentreY-stimulus.pairEccentricity_pix);
            
        elseif thisColour == 2  % Adding half the eccentricity
            
            colourRects(:,thisPair*2) = CenterRectOnPoint(colourRect, theseColourX(thisPair)+screenCentreX, theseColourY(thisPair)+screenCentreY+stimulus.pairEccentricity_pix);

        end
        
    end
    
end

% Build probe rects of the response screen

probeRect = [0 0 stimulus.size_pix stimulus.size_pix];
probeRects = colourRects;           % Should be the same as the colour rect screen

% Build location rects of the reference colours

refRect = [0 0 stimulus.size_pix stimulus.size_pix];        % Same size as colour rects
refRects = NaN(4,stimulus.nColours);
refJitter = [-3.5 -2.5 -1.5 -.5 +.5 +1.5 +2.5 +3.5];        % Reference of the colour blocks

for thisItem = 1:stimulus.nColours
    
    refRects(:,thisItem) = CenterRectOnPoint(refRect, screenCentreX+(refJitter(thisItem)*stimulus.refEccentricity_pix),screenCentreY+4*stimulus.refEccentricity_pix);
   
end

% Build location rects of the reference numbers

numRect = [0 0 stimulus.size_pix stimulus.size_pix];
numRects = NaN(4,stimulus.nColours);

for thisItem = 1:stimulus.nColours
    
    numRects(:,thisItem) = CenterRectOnPoint(numRect, screenCentreX+(refJitter(thisItem)*stimulus.refEccentricity_pix),screenCentreY+4.5*stimulus.refEccentricity_pix);
   
end

% Show instruction text

% ----------------------------------------------------------------------- %

%                       Start Experiment Loop                             %

% ----------------------------------------------------------------------- %

for thisBlock = 1:experiment.nBlocks
    
    % Set up what will be cued for this block
    
    % Retrieve block parameters
    
    block.thisCondition = participant.condnOrder(thisBlock);    % Retrieve whether condition is uniform (1) or pattern (2)
    block.thisShape = participant.shapeOrder(thisBlock);        % Retrieve whether shape is square (1) or circle (2)
    
    % Set up block parameters
    
    block.allPairs = NaN(experiment.nTrialsPerBlock,stimulus.nPairs);
    block.allColours = NaN(experiment.nTrialsPerBlock,stimulus.nPairs,2);
    
    block.allProbes = mod(randperm(experiment.nTrialsPerBlock),8)+1;    % Generate a vector of which item will be probed on each trial in this block
    
    block.allResponseKey = NaN(experiment.nTrialsPerBlock,1);     % Creating an array to save responses (can use a blank vector as well)
    block.allResponseColour = NaN(experiment.nTrialsPerBlock,1);
    block.allCorrect = NaN(experiment.nTrialsPerBlock,1);
    
    for thisTrial = 1:experiment.nTrialsPerBlock
    
        % Select the colour pairs to be displayed on each trial of the block
        
        findUniqueColours = 1;

        while findUniqueColours == 1

            findUniquePairs = 1;

            while findUniquePairs == 1

                if block.thisCondition == 1         % Uniform condition

                    block.allPairs(thisTrial,:) = randsample(1:64,4,true,stimulus.unifVector);

                elseif block.thisCondition == 2     % Pattern condition

                    block.allPairs(thisTrial,:) = randsample(1:64,4,true,stimulus.probVector);          % Samples pairs based on assigned probabilities

                end

                if length(unique(block.allPairs(thisTrial,:))) == 4                                 % Ensures four unique numbers sampled

                    break

                end

            end

            block.allColours(thisTrial,:,:) = stimulus.pairVector(block.allPairs(thisTrial,:),:);   % Retrieves the colour pairs that have been sampled

            if length(unique(block.allColours(thisTrial,:))) == 8                                   % Ensures eight unique colours sampled

                break

            end

        end
        
    end
    
    % Show block text
    
    
    
% ----------------------------------------------------------------------- %    
    
%                         TRIAL LOOP STARTS                               %                             

% ----------------------------------------------------------------------- %
   
    for thisTrial = 1:experiment.nTrialsPerBlock
        
        % Start trial with fixation point/blank
        
        if stimulus.fixationOn
                
            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);
        
        end
        
        startTrialTime = Screen('Flip',ptbWindow);
        
        % Build the stimuli display on each trial
                        
            % Retrieve this trial's colour pairs and allocate to the location rects
            
                % Turning the trial's colour pairs into a vector
            
                thisTrialColourPairs = squeeze(block.allColours(thisTrial,:,:))';
                thisTrialColours = reshape(thisTrialColourPairs,1,stimulus.nColours);
            
            for thisItem = 1:stimulus.nColours
                    
                thisColour = thisTrialColours(thisItem);
                
                    if block.thisShape == 1
                                              
                        Screen('FillRect',ptbWindow,[stimulus.colours(thisColour,:)],colourRects(:,thisItem));
                        
                    elseif block.thisShape == 2
                        
                        Screen('FillOval',ptbWindow,[stimulus.colours(thisColour,:)],colourRects(:,thisItem));
                        
                    end
                    
            end
            
        % Draw stimuli display

        if stimulus.fixationOn
                
            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);
        
        end
        
        stimulusDisplayTime = Screen('Flip',ptbWindow,startTrialTime+timing.ITI);
        
        % Flip to blank for delay
        
        if stimulus.fixationOn
                
            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);
        
        end
        
        responseDelayTime = Screen('Flip',ptbWindow,stimulusDisplayTime+timing.memory);
        
        % Build probe display
        
            % Retrieve item to be probed
            
        thisProbeLocation = block.allProbes(thisTrial);
        thisProbeColour = thisTrialColours(thisProbeLocation);
        
            % Build response screen
            
                % Building the memoranda display
            
        for thisItem = 1:stimulus.nColours

            if block.thisShape == 1
            
                if thisItem ~= thisProbeLocation           % Not the probe item

                    Screen('FrameRect',ptbWindow,colour.probeVal,probeRects(:,thisItem),stimulus.thinPenWidth);

                elseif thisItem == thisProbeLocation       % Is the probe item

                    Screen('FrameRect',ptbWindow,colour.probeVal,probeRects(:,thisItem),stimulus.thickPenWidth);

                end

            elseif block.thisShape == 2
                
                if thisItem ~= thisProbeLocation
                    
                    Screen('FrameOval',ptbWindow,colour.probeVal,probeRects(:,thisItem),stimulus.thinPenWidth);
                    
                elseif thisItem == thisProbeLocation
                    
                    Screen('FrameOval',ptbWindow,colour.probeVal,probeRects(:,thisItem),stimulus.thickPenWidth);
                    
                end
                
            end
                
        end

            % Building the reference display
                
        for thisColour = 1:stimulus.nColours
            
            if block.thisShape == 1
                
                Screen('FillRect',ptbWindow,stimulus.colours(thisColour,:),refRects(:,thisColour));

            elseif block.thisShape == 2
                
                Screen('FillOval',ptbWindow,stimulus.colours(thisColour,:),refRects(:,thisColour));
                
            end
            
            DrawFormattedText(ptbWindow,num2str(thisColour),'center','center',colour.textVal,[],[],[],[],[],numRects(:,thisColour)');
            
        end

            % Draw the screen
        
        if stimulus.fixationOn
                
            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);
        
        end
            
        responseTime = Screen('Flip',ptbWindow,responseDelayTime+timing.delay);
        
        % Get participant response
        
        waitResponse = 1;

        while waitResponse

            [keySecs, keyCode] = KbWait(-1,2);
            pressedKey = find(keyCode);

            if isempty(find(equipment.responseKeys == pressedKey)) == 0     % While loop will only break when a response Key is pressed
                waitResponse = 0;
            end

        end
        
        % Save response
        
        block.allResponseKey(thisTrial) = pressedKey;
        block.allResponseColour(thisTrial) = pressedKey - 29;
        
        % Check if correct
        
        if block.allResponseColour(thisTrial) == thisProbeColour
            
            block.allCorrect(thisTrial) = 1;
            
        else
            
            block.allCorrect(thisTrial) = 0;
            
        end
        
    end
    
    % Save a block file
    
    cd(saveDirectory);
    blockFileName = [participant.ID '_VWMStatLearning_' num2str(thisBlock) '.mat'];
    save(blockFileName, 'block', 'experiment', 'equipment', 'colour', 'stimulus', 'timing', 'participant')
    
end

% Save user file

cd(userDirectory);
userFileName = [participant.ID '_VWMStatLearning.mat'];
save(userFileName, 'participant', 'experiment', 'equipment', 'colour', 'stimulus', 'timing');

