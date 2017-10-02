function runExp2(noCST, noPST, noDT, noPDT, RANDOMIZE)
%% DUAL ATTENTIONAL (VISUAL DISCRIMINATION) TASK %%
% Central Task = Display of letters in circle (4 Conditions: All L's, All T's, 4 L's
% and 1 T or 4 T's and 1 L) before a mask (T = Central SOA, cSOA) of all F's
% Perihperal Task = Presentation of two halves of a vertically bisected disk (either R/G or G/R) at one
% point in the periphery before a random mondrian mask (T = Peripheral SOA, pSOA) is presented.

% noCST: Number of blocks in central single task condition
% noPST: Number of blocks in peripheral single task condition
% noDT: Number of blocks in dual task condition
% noPDT: Number of blocks in dual task condition with partial report
% RANDOMIZE: Boolean, randomize blocks YES/NO
%
% NB. This experiment is designed for 60Hz presentation!

dbstop if error

nBlocks = noCST + noPST + noDT + 2*noPDT;

UseQUEST = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% LOAD EXPERIMENTAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Exp_Parameters2;

tic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% LOAD PREVIOUS QUESTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path = ['../data/raw/Exp2/' Gral.subjNo '_' Gral.subjID '/' Gral.subjNo '_' Gral.subjID '_']; %#ok<*NODEF>
prevSession = num2str(str2num(Gral.session)-1); %#ok<*ST2NM>
prevRun = num2str(str2num(Gral.run)-1);

if str2num(Gral.session) > 1 && str2num(Gral.run) == 1
    clear q p;
    if exist([path prevSession '_4.mat'],'file')
        load([path prevSession '_4.mat'], 'q', 'p');
    else
        load([path prevSession '_3.mat'], 'q', 'p');
    end
elseif str2num(Gral.session) > 1 && str2num(Gral.run) > 1
    clear q p;
    load([path Gral.session '_' prevRun '.mat'], 'q', 'p');
elseif str2num(Gral.session) == 1 && str2num(Gral.run) > 1
    clear p q;
    load([path Gral.session '_' prevRun '.mat'], 'q', 'p');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

blocks = zeros(1, nBlocks);
blocks(1:noCST) = 1;
blocks((noCST + 1):(noCST + noPST)) = 2;
blocks((noCST + noPST + 1):(noCST + noPST + noDT)) = 3;
blocks((noCST + noPST + noDT + 1):(noCST + noPST + noDT + 2*noPDT)) = 4;

if RANDOMIZE
    random = randperm(nBlocks);
    blocks = blocks(random);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b = 1:nBlocks
    
    % Determine condition
    cond = blocks(b);
    
    % Create new QUEST if trialCount > 3*nTrials
    if q.trialCount >= 3*nTrials && cond == 1
        cSOA_estim = QuestMean(q);
        SDcSOA_guess = 3;
        beta = q.beta;
        delta = q.delta;
        gamma = q.gamma;
        q=QuestCreate(cSOA_estim,SDcSOA_guess,pThreshold,beta,delta,gamma,1,50);
        q.normalizePdf=1;
        
    end
    
    if p.trialCount >= 3*nTrials && cond == 2
        pSOA_estim = QuestMean(p);
        SDpSOA_guess = 3;
        beta = p.beta;
        delta = p.delta;
        gamma = p.gamma;
        p=QuestCreate(pSOA_estim,SDpSOA_guess,pThreshold,beta,delta,gamma,1,50);
        p.normalizePdf=1;
    end
    
    
    % Show instructions
    show_instructions(cond, Cfg);
    
    % Remove old data
    
    empty = cell(1,nTrials);
    [TR(:).c_keyid] = empty{:};
    [TR(:).c_confidence] = empty{:};
    [TR(:).mouseResponsesMain] = empty{:};
    [TR(:).c_response] = empty{:};
    
    [TR(:).p_keyid] = empty{:};
    [TR(:).p_confidence] = empty{:};
    [TR(:).mouseResponsesPer] = empty{:};
    [TR(:).p_response] = empty{:};
    
    
    
    % Randomize Trials
    randi = randperm(nTrials);
    TR = TR(randi);
    
    for tr = 1:nTrials
        
        % Retrieve previus Quest estimates
        TR(tr).cSOA = round(QuestMean(q));
        TR(tr).true_cSOA = QuestMean(q);
        TR(tr).pSOA = round(QuestMean(p));
        TR(tr).true_pSOA = QuestMean(p);
        
        % Run Trials
        TR = show_trial(TR, Cfg, tr, cond);
        
        % Update Quests
        if UseQUEST
            
            if cond == 1
                q=QuestUpdate(q,TR(tr).cSOA,TR(tr).c_response);
            elseif cond == 2
                p=QuestUpdate(p,TR(tr).pSOA,TR(tr).p_response);
            end
            
        end
        
        
    end
    
    
    % Store current block and condition in data
    Data(b).TR = TR; %#ok<*AGROW>
    Data(b).condition = cond;
    Data(b).estim_cSOA = [];
    Data(b).estim_pSOA = [];
    Data(b).c_performance = [];
    Data(b).p_performance = [];
    
    % Store performance and SOAs
    switch cond
        case 1
            Data(b).c_performance = sum([TR(:).c_response])/nTrials;
            Data(b).estim_cSOA = QuestMean(q);
        case 2
            Data(b).p_performance = sum([TR(:).p_response])/nTrials;
            Data(b).estim_pSOA = QuestMean(p);
        case 3
            Data(b).c_performance = sum([TR(:).c_response])/nTrials;
            Data(b).p_performance = sum([TR(:).p_response])/nTrials;
            Data(b).estim_cSOA = QuestMean(q);
            Data(b).estim_pSOA = QuestMean(p);
        case 4
            Data(b).c_performance = sum([TR(:).c_response])/length([TR(:).c_response]);
            Data(b).estim_cSOA = QuestMean(q);
            Data(b).p_performance = sum([TR(:).p_response])/length([TR(:).p_response]);
            Data(b).estim_pSOA = QuestMean(p);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gral.exptotalDuration = toc;

fileName_trial = [path Gral.session '_' Gral.run, '.mat'];
save(fileName_trial, 'Data', 'p', 'q', 'Cfg', 'Gral');

% Last flip to show the final painted square on the screen
Screen('FillRect', Cfg.windowPtr, [100 100 100]); %0 - black or gray colour
for m = 1 : TR(tr).screenInterval % In Frames! (eg 60hz = 15 = 250ms)
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
end

Screen('TextSize',Cfg.windowPtr,20);
Screen('DrawText', Cfg.windowPtr, 'End of Session.', Cfg.xCentre-80, Cfg.yCentre, [255 255 255]);
for m = 1 : TR(tr).screenInterval % In Frames! (eg 60hz = 15 = 250ms)
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
end
WaitSecs(1);
KbWait;


%% close all windows / devices opened
ShowCursor;
sca

end

%% Function definitions


function TR = show_trial(TR, Cfg, tr, cond)


HideCursor;



%% Draw Fixation Cross on Centre of Screen
%Draw the lines
Screen('DrawLines', Cfg.windowPtr, Cfg.crossLines, Cfg.crossWidth, Cfg.crossColour, [Cfg.xCentre, Cfg.yCentre]);

% Duration of fixation cross presentation: 300 +/- 100 ms
jitter = randi([0 2],1);
TR(tr).crosstrialDur = TR(tr).crosstrialDur + jitter * 6;

for m = 1 : (TR(tr).crosstrialDur) % In Frames! (eg 60hz = 15 = 250ms)
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
end

Screen('FillRect', Cfg.windowPtr, [100 100 100]);
for m = 1 : TR(tr).intertrialInt % In Frames! (eg 60hz = 15 = 250ms)
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% PRESENT STUMULI DEPENDING ON cSOA and pSOA %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine number of frames to be presented
if TR(tr).pSOA >= TR(tr).cSOA
    nFrames = TR(tr).pSOA + 8; %6 ;%4 %3 %más frames para que dure más el experimento.
else
    nFrames = TR(tr).cSOA + 8; %5 ; %3 %2
end

% Set parameters for presentation 
%nFrames es el tiempo total del experimento.

fmask_start1 = 1; %2
fmask_dur = 2; %(nFrames - (TR(tr).pSOA + 1)/2)/2; %Antes había un entre dos, ahora ya no. Dividir todo entre dos

%----
l_start = 4; %1
l_dur = (TR(tr).cSOA)  %+ 1); %agregué todas las "f's"
% l_dur = TR(tr).cSOA + fmask_dur + f_start + f_dur; %agregué todas las "f's"

f_start = 1 + fmask_dur + 1; % 2 + fmask_dur %Inicio del experimento es igual a inicio y duración de la máscara.
f_dur = (TR(tr).pSOA) %+ 2); %+ 1 %duración del experimento.
% ---

lmask_start = TR(tr).cSOA + 1;
lmask_dur = (nFrames - TR(tr).cSOA);

fmask_start2 = (TR(tr).cSOA + 1)/2 %f_start + f_dur + 1; %Inicio más duración del experimento.
fmask_dur2 = nFrames - fmask_start2; %agregamos esto.
% l_dur = TR(tr).cSOA;

% MAKE MASK TEXTURES FOR THIS TRIAL
% Specify circle locations

circ_left = TR(tr).rectCirclePeriph(1);
circ_right = TR(tr).rectCirclePeriph(3);
circ_up = TR(tr).rectCirclePeriph(2);
circ_down = TR(tr).rectCirclePeriph(4);

% Start looping through frames
co = 1;

for f = 1:nFrames
    
    % Clear screen
    Screen('FillRect', Cfg.windowPtr, [100 100 100]); %0, now Screen('FillRect', Cfg.windowPtr, [100 100 100])
    
%     if f >= f_start && f < (f_start + f_dur)
if f >= l_start && f < (l_start + l_dur)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Draw Ts and Ls (DEP ON targetTrialType SPECIFIED ABOVE) %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Screen('DrawTextures', Cfg.windowPtr, TR(tr).xTexture, [], [TR(tr).centreLetterRect; TR(tr).rectsMain]', TR(tr).ang(1:5));
    
    %Sets the rewriting of one of the letters (either T or L) to the opposite)
    if TR(tr).targetTrialType == 2
        Screen('DrawTexture', Cfg.windowPtr, TR(tr).T_texture, [], TR(tr).rectsMain(TR(tr).reWriteLetterPos,:), TR(tr).ang(TR(tr).reWriteLetterPos+1));
    elseif TR(tr).targetTrialType == 3
        Screen('DrawTexture', Cfg.windowPtr, TR(tr).L_texture, [], TR(tr).rectsMain(TR(tr).reWriteLetterPos,:), TR(tr).ang(TR(tr).reWriteLetterPos+1));
    end
    
end
    
    %Imreads outside of the loops! (for peripheral masks)
    
    fileName = ['Copy_of_faces/' num2str(TR(tr).gender_mask) num2str(TR(tr).picNo2_mask), '.jpg'];
    dataStruct = imread(fileName);
    dims = size(fileName);
    
    faceData = Screen('MakeTexture', Cfg.windowPtr, dataStruct); % faceData = Screen('MakeTexture', Cfg.windowPtr, dataStruct);
    
    fileName2 = ['Copy_of_faces/' num2str(randi(300)) num2str(TR(tr).picNo2_mask), '.jpg'];
    dataStruct2 = imread(fileName2);
    faceData2 = Screen('MakeTexture', Cfg.windowPtr, dataStruct2);
    
    if f >= fmask_start1 && f < (fmask_start1 + fmask_dur)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% SANDWICHI RIAN PERIPHERAL DISK MASK PRESENTATION %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Peripheral Mask ( displayed in periphery)
        % Set up mask for peripheral Task WOOOOOOO
        
        Screen('DrawTexture', Cfg.windowPtr, faceData, [], [TR(tr).rectCirclePeriph]);
        
        Screen('DrawTexture', Cfg.windowPtr, faceData2, [], [Cfg.width - circ_right,circ_up,Cfg.width - circ_left,circ_down]);
        
        
    end
    
    if f >= f_start && f < (f_start + f_dur)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% INSERT PERIPHERAL DISK PRESENTATION %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if TR(tr).diskType == 1
            TR(tr).rectCirclePeriph
            Screen('FillArc', Cfg.windowPtr, Cfg.green, [Cfg.width - circ_right,circ_up,Cfg.width - circ_left,circ_down], 0, 180);
            Screen('FillArc', Cfg.windowPtr, Cfg.red, [TR(tr).rectCirclePeriph], 180, 180);
        elseif TR(tr).diskType == 2
            Screen('FillArc', Cfg.windowPtr, Cfg.red, [Cfg.width - circ_right,circ_up,Cfg.width - circ_left,circ_down], 0, 180);
            Screen('FillArc', Cfg.windowPtr, Cfg.green, [TR(tr).rectCirclePeriph], 180, 180);
        end
        
    end
    
    
    %Set parameters for the other Mondrian peripheral disk mask
    %presentation :)- % Set up mask for peripheral Task
    
    if co == 1
        n_mask1 = randi(300); 
        n_mask2 = randi(300); 
        n_mask3 = randi(300); 
        n_mask4 = randi(300); 
        
   co = 0;
   
    end
   
    fileName = ['Copy_of_faces/' num2str(n_mask1) num2str(TR(tr).picNo2_mask), '.jpg'];
        dataStruct = imread(fileName);
        faceData = Screen('MakeTexture', Cfg.windowPtr, dataStruct);
        
        fileName2 = ['Copy_of_faces/' num2str(n_mask2) num2str(TR(tr).picNo2_mask), '.jpg'];
        dataStruct2 = imread(fileName2);
        faceData2 = Screen('MakeTexture', Cfg.windowPtr, dataStruct2);
        
        fileName3 = ['Copy_of_faces/' num2str(n_mask3) num2str(TR(tr).picNo2_mask), '.jpg'];
        dataStruct3 = imread(fileName);
        faceData3 = Screen('MakeTexture', Cfg.windowPtr, dataStruct);
        
        fileName4 = ['Copy_of_faces/' num2str(n_mask4) num2str(TR(tr).picNo2_mask), '.jpg'];
        dataStruct4 = imread(fileName2);
        faceData4 = Screen('MakeTexture', Cfg.windowPtr, dataStruct2);
        
    if f >= fmask_start2 && f < (fmask_start2 + fmask_dur2)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% PERIPHERAL DISK MASK PRESENTATION (Mondrian) %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Peripheral Mask (displayed in periphery)
        % fileName = ['Copy_of_faces/' num2str(TR(tr).gender_mask) num2str(TR(tr).picNo2_mask), '.jpg'];
        
        Screen('DrawTexture', Cfg.windowPtr, faceData, [], [TR(tr).rectCirclePeriph]);
        
        Screen('DrawTexture', Cfg.windowPtr, faceData2, [], [Cfg.width - circ_right,circ_up,Cfg.width - circ_left,circ_down]);
        
    end
    
    if f >= lmask_start && f < (lmask_start + lmask_dur)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%% DRAW LETTER MASKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Set mask text images at same angles
        
        Screen('DrawTextures', Cfg.windowPtr, TR(tr).F_texture, [], [TR(tr).centreLetterRect; TR(tr).rectsMain]', TR(tr).ang(1:5));
    end
    
    % Present stimuli on screen
    Screen('Flip', Cfg.windowPtr);  
end

%% intermediate screen

Screen('FillRect', Cfg.windowPtr, [100 100 100]); % 0, now Screen('FillRect', Cfg.windowPtr, [100 100 100])
for m = 1 : TR(tr).screenInterval % In Frames! (eg 60hz = 15 = 250ms)
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% COLLECT RESPONSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ShowCursor;



if cond == 1 || cond == 3 || (cond == 4 && TR(tr).cond_PR == 1)
    responseType = 1; % Need this for DrawResponseScreen script 
    
    DrawResponseScreen2;
    Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
    WaitSecs(.3);
    
    clicks = 0; % flags for the responses: 2afc left, 2afc right, 4 conf left, 4 conf right
    
    % Wait until subject has given a response
    while clicks == 0
        
        [x, y] = getMouseResponse();
        
        % Check whether the click went inside a box area
        for m = 1 : size(polyL, 1)
            idxs_left(m) = inpolygon(x,y,squeeze(polyL(m,1,:)),squeeze(polyL(m,2,:)));
            
            idxs_right(m) = inpolygon(x,y,squeeze(polyR(m,1,:)),squeeze(polyR(m,2,:)));
        end
        
        idx_pos_left = find(idxs_left == 1);
        idx_pos_right = find(idxs_right == 1);
        
        % Left boxes click
        if length(idx_pos_left) == 1 %~isempty(idx_pos_left)
            keyid = 1;
            keyid2 = idx_pos_left;
            
            clicks = 1;
            
            % Paint selected box blue
            Screen('FillPoly', Cfg.windowPtr, [0 0 255], squeeze(polyL(idx_pos_left,:,:))',1);
            for wait = 1:10
                Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
            end
            
        end
        
        if length(idx_pos_right) == 1 %~isempty(idx_pos_right)
            keyid = 2;
            keyid2 = idx_pos_right;
            
            clicks= 1; 
            
            % Paint selected box blue
            Screen('FillPoly', Cfg.windowPtr, [0 0 255], squeeze(polyR(idx_pos_right,:,:))',1);
            for wait = 1:10
                Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
            end
            
        end
    end
    
    % Check response
    if keyid == 1
        response = 'same';
    elseif keyid == 2
        response = 'different';
    end
    
    if TR(tr).targetTrialType < 2
        trialType = 'same';
    else
        trialType = 'different';
    end
    
    TR(tr).c_keyid = keyid;
    TR(tr).c_response = strcmp(response, trialType);
    TR(tr).c_confidence = keyid2;
    
    
    
    TR(tr).mouseResponsesMain = [x y];
end

% Interstimulus interval
if cond == 3
    Screen('FillRect', Cfg.windowPtr, [100 100 100]); %0, now Screen('FillRect', Cfg.windowPtr, [100 100 100])
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
    
end

if cond == 2 || cond == 3 || (cond == 4 && TR(tr).cond_PR == 2)
    responseType = 2; %#ok<*NASGU> % Need this for DrawResponseScreen script
    
    DrawResponseScreen2;
    
    Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
    
    
    clicks = 0; % flags for the responses: 2afc left, 2afc right, 5 conf left, 5 conf right
    
    while clicks == 0
        
        [x, y] = getMouseResponse();
        
        % Check whether the click went inside a box area
        for m = 1 : size(polyL, 1)
            idxs_left(m) = inpolygon(x,y,squeeze(polyL(m,1,:)),squeeze(polyL(m,2,:)));
            
            idxs_right(m) = inpolygon(x,y,squeeze(polyR(m,1,:)),squeeze(polyR(m,2,:)));
        end
        
        idx_pos_left = find(idxs_left == 1);
        idx_pos_right = find(idxs_right == 1);
        
        % Left boxes click
        if length(idx_pos_left) == 1
            keyid = 1;
            keyid2 = idx_pos_left;
            
            clicks = 1;
            
            % Paint selected box blue
            Screen('FillPoly', Cfg.windowPtr, [0 0 255], squeeze(polyL(idx_pos_left,:,:))',1);
            for wait = 1:10
                Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
            end
            
        end
        
        if length(idx_pos_right) == 1
            keyid = 2;
            keyid2 = idx_pos_right;
            
            clicks= 1;
            
            % Paint selected box blue
            Screen('FillPoly', Cfg.windowPtr, [0 0 255], squeeze(polyR(idx_pos_right,:,:))',1);
            for wait = 1:10
                Screen('Flip', Cfg.windowPtr,  [], Cfg.aux_buffer);
            end
            
        end
    end
    
    TR(tr).p_keyid = keyid;
    TR(tr).p_confidence = keyid2;
    TR(tr).mouseResponsesPer = [x y];
    
    % Check response
    
    TR(tr).p_response = (keyid == TR(tr).diskType);
    
end


HideCursor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inter trial interval

Screen('FillRect', Cfg.windowPtr, [100 100 100]); %0, now Screen('FillRect', Cfg.windowPtr, [100 100 100])
for m = 1 : TR(tr).intertrialInt % In Frames! (eg 60hz = 120 = 2 secs)
    Screen('Flip', Cfg.windowPtr, [], Cfg.aux_buffer);
end


%% wait for click to proceed

if tr ~= TR(tr).nTrials
    Screen('TextSize',Cfg.windowPtr,20);
    Screen('TextColor', Cfg.windowPtr ,[255 255 255]);
    DrawFormattedText(Cfg.windowPtr, '<<Click to proceed to the next trial>>','center', 'center');
    Screen('Flip', Cfg.windowPtr, [], []);
    % Wait for mouse click
    [~,~,buttons] = GetMouse;
    while any(buttons) % if already down, wait for release
        [~,~,buttons] = GetMouse;
    end
    while ~any(buttons) % wait for press
        [~,~,buttons] = GetMouse;
    end
    while any(buttons) % wait for release
        [~,~,buttons] = GetMouse;
    end
end

end


function show_instructions(cond, Cfg)

Screen(Cfg.windowPtr, 'TextSize', 20);

if cond == 1
    
    instr1 = 'Please focus on the centre of the screen and decide whether the letters first presented to you are the same or different.';
    
    DrawFormattedText(Cfg.windowPtr, 'Central Task', 'center', 400, [255 255 255], 140, [], [], 2);
    DrawFormattedText(Cfg.windowPtr, instr1, 'center', 500, [255 255 255], 140, [], [], 2);
    DrawFormattedText(Cfg.windowPtr, '<<Click to begin with the first trial>>','center', 1000);
    Screen('Flip', Cfg.windowPtr, [], []);
    WaitSecs(.3);
    % wait for click to proceed
    [~,~,buttons] = GetMouse;
    while any(buttons) % if already down, wait for release
        [~,~,buttons] = GetMouse;
    end
    while ~any(buttons) % wait for press
        [~,~,buttons] = GetMouse;
    end
    while any(buttons) % wait for release
        [~,~,buttons] = GetMouse;
    end
    
    
elseif cond == 2
    
    instr2 = 'Please focus on the centre of the screen and decide whether the disk first presented to you in the periphery is red/green or green/red.';
    
    DrawFormattedText(Cfg.windowPtr, 'Peripheral Task', 'center', 400, [255 255 255], 80, [], [], 2);
    DrawFormattedText(Cfg.windowPtr, instr2, 'center', 500, [255 255 255], 80, [], [], 2);
    DrawFormattedText(Cfg.windowPtr, '<<Click to begin with the first trial>>','center', 750);
    Screen('Flip', Cfg.windowPtr, [], []);
    WaitSecs(.3);
    % wait for click to proceed
    [~,~,buttons] = GetMouse;
    while any(buttons) % if already down, wait for release
        [~,~,buttons] = GetMouse;
    end
    while ~any(buttons) % wait for press
        [~,~,buttons] = GetMouse;
    end
    while any(buttons) % wait for release
        [~,~,buttons] = GetMouse;
    end
    
    
elseif cond == 3
    
    instr3 = 'Please focus on the centre of the screen and decide whether the letters first presented to you in the centre are the same or different as well as whether the disk first presented to you in the periphery is red/green or green/red. After each trial you will be asked to respond both to the letter and the disk task.';
    
    DrawFormattedText(Cfg.windowPtr, 'Dual Task', 'center', 400, [255 255 255], 80, [], [], 2);
    DrawFormattedText(Cfg.windowPtr, instr3, 'center', 500, [255 255 255], 80, [], [], 2);
    DrawFormattedText(Cfg.windowPtr, '<<Click to begin with the first trial>>','center', 750);
    Screen('Flip', Cfg.windowPtr, [], []);
    WaitSecs(.3);
    % wait for click to proceed
    [~,~,buttons] = GetMouse;
    while any(buttons) % if already down, wait for release
        [~,~,buttons] = GetMouse;
    end
    while ~any(buttons) % wait for press
        [~,~,buttons] = GetMouse;
    end
    while any(buttons) % wait for release
        [~,~,buttons] = GetMouse;
    end
    
elseif cond == 4
    
    instr3 = 'Please focus on the centre of the screen and decide whether the letters first presented to you in the centre are the same or different as well as whether the disk first presented to you in the periphery is red/green or green/red. After each trial you will be asked to respond either to the letter task or the disk task.';
    
    DrawFormattedText(Cfg.windowPtr, 'Dual Task with Partial Report', 'center', 400, [255 255 255], 80, [], [], 2);
    DrawFormattedText(Cfg.windowPtr, instr3, 'center', 500, [255 255 255], 80, [], [], 2);
    DrawFormattedText(Cfg.windowPtr, '<<Click to begin with the first trial>>','center', 750);
    Screen('Flip', Cfg.windowPtr, [], []);
    WaitSecs(.3);
    % wait for click to proceed
    [~,~,buttons] = GetMouse;
    while any(buttons) % if already down, wait for release
        [~,~,buttons] = GetMouse;
    end
    while ~any(buttons) % wait for press
        [~,~,buttons] = GetMouse;
    end
    while any(buttons) % wait for release
        [~,~,buttons] = GetMouse;
    end
    
    
end

end