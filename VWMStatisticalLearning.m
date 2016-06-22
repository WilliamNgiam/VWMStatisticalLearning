% Statistical Learning in Visual Working Memory

% This experiment is designed to examine the influence of statistical
% learning on visual working memory capacity. This experiment is a
% replication of the retracted Anderson paper

% WXQN started writing this on 1/4/16
% WXQN started piloting on 1/5/16
% WXQN set up code for UChicago on 22/6/16

% -------------------------------------------------------------------------

clear all;
Screen('CloseAll');

% Set up experiment parameters

experiment.nBlocks = 20;            % Number of blocks
experiment.nTrialsPerBlock = 60;    % Number of trials per block

% Set up participant parameters

userDirectory = 'C:\Users\AwhVogelLab\Documents\MATLAB\Will\VWMStatisticalLearning\UserData';
saveDirectory = 'C:\Users\AwhVogelLab\Documents\MATLAB\Will\VWMStatisticalLearning\Data';
experimentDirectory = 'C:\Users\AwhVogelLab\Documents\MATLAB\Will\VWMStatisticalLearning';    

% Set up equipment parameters

equipment.viewDist = 500;           % Viewing distance in mm
equipment.ppm = 3.6;                % Pixels per mm - Measured at UChicago on 22/6/16
equipment.gammaVals = 1.0./[3.0902 2.4049 2.3194];      % Gamma values for CRT in GT519 (recalibrated 25/9/15)

equipment.greyVal = .5;             % Grey background
equipment.blackVal = 0;
equipment.whiteVal = 1;

equipment.responseKeys = 49:56;     % Responses key for 1 to 8 on Windows Keyboard (30:37 on Macbook)

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
stimulus.fixationEccentricity_dva = 1.7;    % Eccentricity of stimulus from fixation point
stimulus.pairEccentricity_dva = 1.0;        % Eccentricity between (center of) colour pairs
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
stimulus.thickPenWidth = 3;         % Thick line for the memoranda rect that is the probe

% Set up temporal parameters

timing.ITI = .75;       % Inter-trial interval
timing.memory = 1;
timing.delay = 1;

% Get participant ID

while true
    participant.ID = upper(input('Enter two-digit participant number: ', 's'));
    if length(participant.ID)==2
        break
    end
end

cd(userDirectory);
participant.userFile = [participant.ID '_VWMStatLearning_Exp1.mat'];
newUser = ~exist(participant.userFile,'file');

% Set up Psychtoolbox Pipeline

AssertOpenGL;

    % Imaging set-up

screenID = max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
Screen('Preference','SkipSyncTests',0);

    % Window set-up
    
[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, equipment.greyVal,[],[],[],[],6);
PsychColorCorrection('SetEncodingGamma', ptbWindow, equipment.gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);

    % Text set-up
    
Screen('TextFont',ptbWindow,'Arial');
Screen('TextSize',ptbWindow,20);
Screen('TextStyle',ptbWindow,1);        % Bold text

global ptb_drawformattedtext_disableClipping;       % Disable clipping of text 
ptb_drawformattedtext_disableClipping = 1;

startExperimentText = ['On each trial, you will be shown eight unique colors for a short time. \n\n' ...
    'You will be prompted to click on which color you think was shown at one of these locations. \n\n' ...
    'Press any key to begin the experiment.'];

startBlockText = ['Press any key to begin the next block.'];

% Enable alpha blending for typical drawing of masked textures
Screen('BlendFunction', ptbWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Enable multisampling to display anti-aliased images
% 
% Screen('OpenWindow',ptbWindow,[],[],[],[],[],6);

% Calculate equipment parameters

equipment.mpd = (equipment.viewDist)*tan(deg2rad(2*stimulus.fixationEccentricity_dva))/stimulus.fixationEccentricity_dva; % Calculate mm per degree of visual angle to the ecccentricity of the stimuli
equipment.ppd = equipment.ppm*equipment.mpd;        % Pixels per degree

% Calculate spatial parameters

stimulus.size_pix = round(stimulus.size_dva*equipment.ppd);                                     % Item size in pixels
stimulus.fixationEccentricity_pix = round(stimulus.fixationEccentricity_dva*equipment.ppd);     % Eccentricity of stimulus in pixels
stimulus.fixationSize_pix = stimulus.fixationSize_dva*equipment.ppd;                            % Fixation dot size in pixels
stimulus.pairEccentricity_pix = round(stimulus.pairEccentricity_dva*equipment.ppd);             % Eccentricity between colour pairs in pixels
stimulus.refEccentricity_pix = round(stimulus.refEccentricity_dva*equipment.ppd);               % Eccentricity between colour references in response screen

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
            
            if whichCondnFirst == 1     % Uniform first
                
                participant.condnOrder = repmat([1,2],1,9);
                participant.condnOrder = [participant.condnOrder 1 1];  % Make the final two blocks uniform conditions
                
            elseif whichCondnFirst == 2 % Pattern first
                
                participant.condnOrder = repmat([2,1],1,9);
                participant.condnOrder = [participant.condnOrder 1 1];  % Make the final two blocks uniform conditions
                
            end     
        
    participant.thisBlock = 1;
    save(participant.userFile, 'participant');
    
else
    
    load(participant.userFile);
    participant.thisBlock = participant.thisBlock+1;
    
end

% Set up location rects

fixRect = [0 0 stimulus.fixationSize_pix stimulus.fixationSize_pix];    % Fixation rect
fixRect = CenterRectOnPoint(fixRect, screenCentreX, screenCentreY);     % Centred in the middle of the screen

colourRect = [0 0 stimulus.size_pix stimulus.size_pix];                 % Colour stimulus rect
awareRect = CenterRectOnPoint(colourRect, screenCentreX, screenCentreY);

% Build location rects of the colours going clockwise starting from north

colourRects = NaN(4,stimulus.nColours);

colourTheta_deg = linspace(0,270,stimulus.nPairs)-90;       %[-90 0 90 180] - Rotation of colour pairs in degrees (North, East, South, West)
colourTheta_rad = deg2rad(colourTheta_deg)';                % Rotation of the colour pairs in radians (North, East, South, West)
[theseColourX, theseColourY] = pol2cart(colourTheta_rad,stimulus.fixationEccentricity_pix*ones(stimulus.nPairs,1));

for thisPair = 1:stimulus.nPairs
    
    for thisColour = 1:2        % Create rect for each colour block in pair, separate them vertically by adding eccentricity to y-value of rect
        
        if thisColour == 1      % Subtracting half the eccentricity
            
            colourRects(:,thisPair*2-1) = CenterRectOnPoint(colourRect, theseColourX(thisPair)+screenCentreX, theseColourY(thisPair)+screenCentreY-(stimulus.pairEccentricity_pix/2));
            
        elseif thisColour == 2  % Adding half the eccentricity
            
            colourRects(:,thisPair*2) = CenterRectOnPoint(colourRect, theseColourX(thisPair)+screenCentreX, theseColourY(thisPair)+screenCentreY+(stimulus.pairEccentricity_pix/2));

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
    
    refRects(:,thisItem) = round(CenterRectOnPoint(refRect, screenCentreX+(refJitter(thisItem)*stimulus.refEccentricity_pix),screenCentreY+3*stimulus.refEccentricity_pix));
   
end

% Build location rects of the reference numbers

numRect = [0 0 stimulus.size_pix stimulus.size_pix];
numRects = NaN(4,stimulus.nColours);

for thisItem = 1:stimulus.nColours
    
    numRects(:,thisItem) = round(CenterRectOnPoint(numRect, screenCentreX+(refJitter(thisItem)*stimulus.refEccentricity_pix),screenCentreY+2.5*stimulus.refEccentricity_pix));
   
end

% Show instruction text

DrawFormattedText(ptbWindow,startExperimentText,'center','center',colour.textVal);
startExperimentTime = Screen('Flip',ptbWindow);
waitResponse = 1;

while waitResponse

    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;

end

Screen('Flip',ptbWindow);

% ----------------------------------------------------------------------- %

%                       Start Experiment Loop                             %

% ----------------------------------------------------------------------- %

for thisBlock = 1:experiment.nBlocks
    
    block.thisBlock = thisBlock;                % Save the block to file.
    block.thisCondnBlock = ceil(thisBlock/2);   % The nth block of condition.
    
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
    block.allResponseTimes = NaN(experiment.nTrialsPerBlock,1);
    
    for thisTrial = 1:experiment.nTrialsPerBlock
    
        % Select the colour pairs to be displayed on each trial of the block
        
        selectedPairs = [];
        
            % Retrieve probability vector and matrix for manipulation
        
        if block.thisCondition == 1         % Uniform condition
            
            trialProbVector = stimulus.unifVector;
            trialProbMatrix = stimulus.unifMatrix;
            
        elseif block.thisCondition == 2     % Pattern condition
            
            trialProbVector = stimulus.probVector;
            trialProbMatrix = stimulus.probMatrix;
            
        end
        
            % Sample one pair at a time
            
        for thisPair = 1:stimulus.nPairs
            
            thisSelectedPair = randsample(1:64,1,true,trialProbVector);
            selectedPairs = [selectedPairs thisSelectedPair];
            
            % Set probability of any pair with the selected pairs' colours
            % to zero
            
            thisPairColours = stimulus.pairVector(thisSelectedPair,:);      % Retrieves the colours selected this pair
            trialProbMatrix(thisPairColours(1),:) = 0;                      % Sets the probabilities in the matrix of each pair ...
            trialProbMatrix(thisPairColours(2),:) = 0;                      % with the selected colour to zero.
            trialProbMatrix(:,thisPairColours(1)) = 0;
            trialProbMatrix(:,thisPairColours(2)) = 0;
            
            % Create new probability vector with newly set zeroes
            
            trialProbVector = [];
            
            for thisRow = 1:stimulus.nColours
                
                trialProbVector = [trialProbVector trialProbMatrix(thisRow,:)];
                
            end
            
        end
        
        block.allPairs(thisTrial,:) = selectedPairs;        % Stores the sampled pairs in to trial matrix for the block
        block.allColours(thisTrial,:,:) = stimulus.pairVector(block.allPairs(thisTrial,:),:);   % Retrieves the colour pairs that have been sampled

    end
    
    % Start block text
    
    DrawFormattedText(ptbWindow,startBlockText,'center','center',colour.textVal);
    startBlockTime = Screen('Flip',ptbWindow);
    waitResponse = 1;

    while waitResponse

        [startTime, keyCode] = KbWait(-1,2);
        waitResponse = 0;

    end
    
    endTrialTime = startTime;
    
% ----------------------------------------------------------------------- %    
    
%                         TRIAL LOOP STARTS                               %                             

% ----------------------------------------------------------------------- %
   
    for thisTrial = 1:experiment.nTrialsPerBlock
        
        % Start trial with fixation point/blank
        HideCursor;
        
        if stimulus.fixationOn
                
            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);
        
        end
        
        startTrialTime = Screen('Flip',ptbWindow,endTrialTime+timing.ITI);
        
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
            
%             DrawFormattedText(ptbWindow,num2str(thisColour),'center','center',colour.textVal,[],[],[],[],[],numRects(:,thisColour)');

        end

            % Draw the screen
        
        if stimulus.fixationOn
                
            Screen('FillOval', ptbWindow, colour.fixVal, fixRect);
        
        end
            
        responseTime = Screen('Flip',ptbWindow,responseDelayTime+timing.delay);
        
        % Get participant response
% For mouse click responses        

        ShowCursor(0);
        SetMouse(screenCentreX,screenCentreY,ptbWindow);
        CheckResponse = zeros(1,stimulus.nColours);
        
        while ~any(CheckResponse)
            
            [~,xClickResponse,yClickResponse] = GetClicks(ptbWindow,0);     % Retrieves x- and y-coordinates of mouse click
            clickSecs = GetSecs;
        
            for thisColour = 1:stimulus.nColours;
                
                CheckResponse(thisColour) = IsInRect(xClickResponse,yClickResponse,refRects(:,thisColour));     % Tests if mouse click is inside aperture of each successive item
                
            end
            
        end
        
        responseColour = find(CheckResponse);

% For keyboard responses   
%        
%         waitResponse = 1;
% 
%         while waitResponse
% 
%             [keySecs, keyCode] = KbWait(-1,2);
%             pressedKey = find(keyCode);
% 
%             if length(pressedKey) ~= 1      % In case, participant accidentally presses two buttons simultaneously
%                 continue
%             end
%                           
%             if isempty(find(equipment.responseKeys == pressedKey)) == 0     % While loop will only break when a response Key is pressed
%                 waitResponse = 0;
%             end
% 
%         end
        
        % Save response
% For mouse click responses
        block.allResponseColour(thisTrial) = responseColour;
        block.allResponseTimes(thisTrial) = clickSecs - responseTime;
        
% For keyboard responses
%        block.allResponseKey(thisTrial) = pressedKey;
%        block.allResponseColour(thisTrial) = pressedKey - 29;
%        block.allResponseTimes(thisTrial) = keySecs - responseTime;
        
        % Check if correct
        
        if block.allResponseColour(thisTrial) == thisProbeColour
            
            block.allCorrect(thisTrial) = 1;
            
        else
            
            block.allCorrect(thisTrial) = 0;
            
        end
        
        endTrialTime = Screen('Flip',ptbWindow);
        
    end
    
    % Save a block file
    
    cd(saveDirectory);
    blockFileName = [participant.ID '_VWMStatLearning_' num2str(thisBlock) '.mat'];
    save(blockFileName, 'block', 'experiment', 'equipment', 'colour', 'stimulus', 'timing', 'participant','startExperimentTime')
    
    % Completed block text

    if mod(block.thisBlock,4) == 0
        
        if block.thisBlock ~= experiment.nBlocks
            
            takeBreakText = ['You have completed ' num2str(block.thisBlock) ' out of ' num2str(experiment.nBlocks) ' blocks.\n\n' ...
                'Please take a break for a few minutes.'];

            DrawFormattedText(ptbWindow,takeBreakText,'center','center',colour.textVal);
            Screen('Flip',ptbWindow);
            WaitSecs(120);

            completedBlockText = ['You have completed ' num2str(block.thisBlock) ' out of ' num2str(experiment.nBlocks) ' blocks.\n\n' ...
                'Press any key to continue.'];
            DrawFormattedText(ptbWindow,completedBlockText,'center','center',colour.textVal);
            Screen('Flip',ptbWindow);
            waitResponse = 1;

            while waitResponse

                [time, keyCode] = KbWait(-1,2);
                waitResponse = 0;

            end
            
        elseif block.thisBlock == experiment.nBlocks
            
            finishAllText = ['You have completed all blocks in the experiment.\n\n' ...
                'Press any key to continue on to the next part.'];
            
        end
        
    else
    
        completedBlockText = ['You have completed ' num2str(block.thisBlock) ' out of ' num2str(experiment.nBlocks) ' blocks.\n\n' ...
            'Press any key to continue.'];

        DrawFormattedText(ptbWindow,completedBlockText,'center','center',colour.textVal);
        Screen('Flip',ptbWindow);
        waitResponse = 1;

        while waitResponse

            [time, keyCode] = KbWait(-1,2);
            waitResponse = 0;

        end
    
    end
    
end

% Test explicit awareness of blocks

% % Instruction text

awarenessText = ['You will be presented with a colour in the middle of the screen.\n\n' ...
    'Click on which colour you think appeared most commonly with that colour.\n\n' ...
    'Press any key to continue.'];

DrawFormattedText(ptbWindow,awarenessText,'center','center',colour.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;

while waitResponse

    [time, keyCode] = KbWait(-1,2);
    waitResponse = 0;

end    

awareness.responses = [];
awareness.responseColours = [];

for thisTestColour = 1:stimulus.nColours
    
    % Retrieve which shape had pattern configurations
    
    if whichCondnFirst == 2
        
        testShape = whichShapeFirst;
        
    elseif whichCondnFirst == 1
        
        testShape = 3-whichShapeFirst;
        
    end
     
    % Draw question text above
    
    questionText = ['Which colour was most likely to appear with this colour shown?'];
    
    DrawFormattedText(ptbWindow,questionText,'center',screenCentreY - 1*stimulus.refEccentricity_pix,colour.textVal);
    
    % Present colour in the middle
    
    if testShape == 1
                                              
        Screen('FillRect',ptbWindow,[stimulus.colours(thisTestColour,:)],awareRect);

    elseif testShape == 2

        Screen('FillOval',ptbWindow,[stimulus.colours(thisTestColour,:)],awareRect);

    end

    % Present reference values
    
    for thisColour = 1:stimulus.nColours
                 
        if testShape == 1

            Screen('FillRect',ptbWindow,stimulus.colours(thisColour,:),refRects(:,thisColour));

        elseif testShape == 2

            Screen('FillOval',ptbWindow,stimulus.colours(thisColour,:),refRects(:,thisColour));

        end

        DrawFormattedText(ptbWindow,num2str(thisColour),'center','center',colour.textVal,[],[],[],[],[],numRects(:,thisColour)');

    end
    
    Screen('Flip',ptbWindow);
    
    % Record response
% For mouse click responses 
    ShowCursor;
    SetMouse(screenCentreX,screenCentreY,ptbWindow);
    CheckResponse = zeros(1,stimulus.nColours);

    while ~any(CheckResponse)

        [~,xClickResponse,yClickResponse] = GetClicks(ptbWindow,0);     % Retrieves x- and y-coordinates of mouse click
        clickSecs = GetSecs;

        for thisColour = 1:stimulus.nColours;

            CheckResponse(thisColour) = IsInRect(xClickResponse,yClickResponse,refRects(:,thisColour));     % Tests if mouse click is inside aperture of each successive item

        end

    end

    responseColour = find(CheckResponse);

% For keyboard responses    
%     waitResponse = 1;
%     
%     while waitResponse
% 
%         [keySecs, keyCode] = KbWait(-1,2);
%         pressedKey = find(keyCode);
% 
%         if length(pressedKey) ~= 1
%             continue
%         end
%         
%         if isempty(find(equipment.responseKeys == pressedKey)) == 0     % While loop will only break when a response Key is pressed
%             waitResponse = 0;
%         end
% 
%     end
% 
     % Save response
% For mouse click responses

    awareness.responseColours = [awareness.responses responseColour];

% For keyboard responses    
%     awareness.responses = [awareness.responses pressedKey];
%     awareness.responseColours = [awareness.responseColours pressedKey-29];
    
end

% Save user file

cd(userDirectory);
userFileName = [participant.ID '_VWMStatLearning_Exp1.mat'];
save(userFileName, 'participant', 'experiment', 'equipment', 'colour', 'stimulus', 'timing', 'awareness');

% Completed Experiment Text

completedBlockText = ['You have completed the experiment.\n' ...
    'Press any key to continue.'];

DrawFormattedText(ptbWindow,completedBlockText,'center','center',colour.textVal);
Screen('Flip',ptbWindow);
waitResponse = 1;

Screen('CloseAll');
close all;

