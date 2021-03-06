%% RESPONSE SCREEN

SetMouse(Cfg.x,Cfg.y);
ShowCursor;

%Draw the response boxes
Screen('FillRect',  Cfg.windowPtr, [100 100 100]); %Screen('FillRect',  Cfg.windowPtr, Cfg.color.black);


Screen('TextSize',Cfg.windowPtr,20);


Screen('FillRect', Cfg.windowPtr, Cfg.color.white, Cfg.bigrect_{1} );
Screen('FillRect', Cfg.windowPtr, [100 100 100], Cfg.smallrect_{1} ); %Cfg.color.black (hasta la l�nea 26, lo cambi� por [100 100 100])
Screen('FillRect', Cfg.windowPtr, [100 100 100], Cfg.cleavage_{1});
Screen('DrawLine', Cfg.windowPtr,[100 100 100],Cfg.smallrect_{1}(1),Cfg.smallrect_{1}(4),Cfg.bigrect_{1}(1),Cfg.bigrect_{1}(4),2);
Screen('DrawLine', Cfg.windowPtr,[100 100 100],Cfg.smallrect_{1}(1),Cfg.smallrect_{1}(2),Cfg.bigrect_{1}(1),Cfg.bigrect_{1}(2),2);
Screen('DrawLine', Cfg.windowPtr,[100 100 100],Cfg.smallrect_{1}(1),Cfg.smallrect_{1}(2),Cfg.bigrect_{1}(1),Cfg.bigrect_{1}(2),2);
Screen('DrawLine', Cfg.windowPtr,[100 100 100],Cfg.smallrect_{1}(3),Cfg.smallrect_{1}(2),Cfg.bigrect_{1}(3),Cfg.bigrect_{1}(2),2);
Screen('DrawLine', Cfg.windowPtr,[100 100 100],Cfg.smallrect_{1}(3),Cfg.smallrect_{1}(4),Cfg.bigrect_{1}(3),Cfg.bigrect_{1}(4),2);

Screen('DrawLine', Cfg.windowPtr,[100 100 100],Cfg.x,Cfg.smallrect_{1}(4),Cfg.x,Cfg.bigrect_{1}(4),2);
Screen('DrawLine', Cfg.windowPtr,[100 100 100],Cfg.x,Cfg.smallrect_{1}(2),Cfg.x,Cfg.bigrect_{1}(2),2);

Screen('DrawLine', Cfg.windowPtr,[100 100 100],Cfg.smallrect_{1}(1),Cfg.y,Cfg.bigrect_{1}(1),Cfg.y,1.3);
Screen('DrawLine', Cfg.windowPtr,[100 100 100],Cfg.smallrect_{1}(3),Cfg.y,Cfg.bigrect_{1}(3),Cfg.y,1.3);

polyL(1,:,:)=[Cfg.smallrect_{1}(1) Cfg.cleavage_{1}(1) Cfg.cleavage_{1}(1) Cfg.bigrect_{1}(1);Cfg.smallrect_{1}(4) Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4) Cfg.bigrect_{1}(4)];
polyL(2,:,:)=[Cfg.smallrect_{1}(1) Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1) Cfg.bigrect_{1}(1);Cfg.smallrect_{1}(4) Cfg.y Cfg.y Cfg.bigrect_{1}(4)];
polyL(3,:,:)=[Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1) Cfg.bigrect_{1}(1) Cfg.smallrect_{1}(1);Cfg.y Cfg.y Cfg.bigrect_{1}(2) Cfg.smallrect_{1}(2)];
polyL(4,:,:)=[Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1) Cfg.cleavage_{1}(1) Cfg.cleavage_{1}(1);Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2) Cfg.bigrect_{1}(2) Cfg.smallrect_{1}(2)];

polyR(1,:,:)=[Cfg.cleavage_{1}(3) Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3) Cfg.cleavage_{1}(3);Cfg.smallrect_{1}(4) Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4) Cfg.bigrect_{1}(4)];
polyR(2,:,:)=[Cfg.smallrect_{1}(3) Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3) Cfg.bigrect_{1}(3);Cfg.smallrect_{1}(4) Cfg.y Cfg.y Cfg.bigrect_{1}(4)];
polyR(3,:,:)=[Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3) Cfg.bigrect_{1}(3) Cfg.smallrect_{1}(3);Cfg.y Cfg.y Cfg.bigrect_{1}(2) Cfg.smallrect_{1}(2)];
polyR(4,:,:)=[Cfg.smallrect_{1}(3) Cfg.cleavage_{1}(3) Cfg.cleavage_{1}(3) Cfg.bigrect_{1}(3);Cfg.smallrect_{1}(2) Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2) Cfg.bigrect_{1}(2)];


% Draw 'sureness'

Screen('TextSize',Cfg.windowPtr,36);

DrawFormattedText(Cfg.windowPtr, 'sure', 'center',  Cfg.bigrect_{1}(2) - 60 , [255 255 255]);
DrawFormattedText(Cfg.windowPtr, 'not sure', 'center',  Cfg.bigrect_{1}(4) + 20 , [255 255 255]);


% Boxes for text
text_box=[0 0 1 1];
offset=(Cfg.x-mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)]))/2.;

global ptb_drawformattedtext_disableClipping; % needed for DrawFormattedText with winRect ... otherwise text is not drawn
ptb_drawformattedtext_disableClipping=1;

rect_1L=text_box+[mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)])+offset mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)]) mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)])+offset mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])];
rect_2L=text_box+[mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)]) mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])-offset mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)]) mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])-offset];
rect_3L=text_box+[mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)]) mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])+offset mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)]) mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])+offset];
rect_4L=text_box+[mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)])+offset mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)]) mean([Cfg.smallrect_{1}(1) Cfg.bigrect_{1}(1)])+offset mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])];

rect_1R=text_box+[mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)])-offset mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)]) mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)])-offset mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])];
rect_2R=text_box+[mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)]) mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])-offset mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)]) mean([Cfg.smallrect_{1}(4) Cfg.bigrect_{1}(4)])-offset];
rect_3R=text_box+[mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)]) mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])+offset mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)]) mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])+offset];
rect_4R=text_box+[mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)])-offset mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)]) mean([Cfg.smallrect_{1}(3) Cfg.bigrect_{1}(3)])-offset mean([Cfg.smallrect_{1}(2) Cfg.bigrect_{1}(2)])];

Screen('TextSize',Cfg.windowPtr,30);

DrawFormattedText(Cfg.windowPtr, '1', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_1L);
DrawFormattedText(Cfg.windowPtr, '2', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_2L);
DrawFormattedText(Cfg.windowPtr, '3', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_3L);
DrawFormattedText(Cfg.windowPtr, '4', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_4L);

DrawFormattedText(Cfg.windowPtr, '1', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_1R);
DrawFormattedText(Cfg.windowPtr, '2', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_2R);
DrawFormattedText(Cfg.windowPtr, '3', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_3R);
DrawFormattedText(Cfg.windowPtr, '4', 'center', 'center',  [0 0 0],[],[],[],[],[],rect_4R);


% Draw text for the left boxes
move_x = 10; 
move_y= 25;

Screen('TextSize',Cfg.windowPtr,40);


if responseType == 1
    Screen('DrawText', Cfg.windowPtr, 'S', Cfg.x-7*move_x, Cfg.y-move_y, [255 255 255]);
    Screen('DrawText', Cfg.windowPtr, 'D',  Cfg.x+4*move_x, Cfg.y-move_y, [255 255 255]);    
elseif responseType == 2
    
        % set position for disks to be presented
        loc_left1 = [Cfg.x-4.5*move_x; Cfg.y]; % ###
        loc_left2= [Cfg.x-5.5*move_x; Cfg.y]; % ###
        loc_right1 = [Cfg.x+5.5*move_x; Cfg.y];
        loc_right2 = [Cfg.x+4.5*move_x; Cfg.y];
        
        disk_loc_left1 = vertcat(loc_left1 - .5*Cfg.pixelsPerDegree,loc_left1 + .5*Cfg.pixelsPerDegree);
        disk_loc_left2 = vertcat(loc_left2 - .5*Cfg.pixelsPerDegree,loc_left2 + .5*Cfg.pixelsPerDegree);
        disk_loc_right1 = vertcat(loc_right1 - .5*Cfg.pixelsPerDegree,loc_right1 + .5*Cfg.pixelsPerDegree);
        disk_loc_right2 = vertcat(loc_right2 - .5*Cfg.pixelsPerDegree,loc_right2 + .5*Cfg.pixelsPerDegree);
 

        Screen('FillArc', Cfg.windowPtr, Cfg.green, disk_loc_left1, 0, 180); % ###
        Screen('FillArc', Cfg.windowPtr, Cfg.red, disk_loc_left2, 180, 180);

        Screen('FillArc', Cfg.windowPtr, Cfg.red, disk_loc_right1, 0, 180);
        Screen('FillArc', Cfg.windowPtr, Cfg.green, disk_loc_right2, 180, 180);

      
end

