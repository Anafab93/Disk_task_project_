% blackback = imread('1.jpg');
% image(blackback)
% axis equal
% axis tight
%
% imshow('1.jpg')
%
% dims = size(blackback) % "y", "x" y "z" (en ese orden)
%
% grayback = blackback;
% dims = size(blackback);
%
% for ii = 1:dims(1)
%     for jj = 1:dims(2)
%         if blackback(ii,jj,3) == 0 && blackback(ii,jj,2) == 0 && blackback(ii,jj,1) == 0
%             grayback(ii,jj,1) = 100;
%             grayback(ii,jj,2) = 100;
%             grayback(ii,jj,3) = 100;
%         end
%     end
% end
%
% imwrite(grayback,'graybackback.jpg','jpg')
%
% imshow('graybackback.jpg')
%% Todos - Mondrian mask

happypath = 'C:\Ana\Julians_task\scripts\faces\';
happypathout = 'C:\Ana\Julians_task\scripts\Copy_of_faces\';

%
% blackback = imread('1.jpg');
% image(blackback)

% imshow('1.jpg')

% "y", "x" y "z" (en ese orden)
% grayback = blackback;
% 
% for ii = 1:300
%     grayback = imread([happypath,num2str(ii),'.jpg']);
%     dims = size(grayback);
%     axis equal
%     axis tight
%     
%     for uu = 1:dims(1)
%         for jj = 1:dims(2)
%             if grayback(uu,jj,3) <= 35 && grayback(uu,jj,2) <= 35 && grayback(uu,jj,1) <= 35
%                 grayback(uu,jj,1) = 100;
%                 grayback(uu,jj,2) = 100;
%                 grayback(uu,jj,3) = 100;
%             end
%         end
%     end
%     imwrite(grayback,[happypathout,num2str(ii),'.jpg'],'jpg')
% end

%% Letter

happypath = 'C:\Ana\Julians_task\scripts\faces\';
happypathout = 'C:\Ana\Julians_task\scripts\Copy_of_faces\';

%
% blackback = imread('1.jpg');
% image(blackback)

% imshow('1.jpg')

% "y", "x" y "z" (en ese orden)
% grayback = blackback;

for ii = 1:3
    switch ii
        case 1
            nombre = 'L';
        case 2
            nombre = 'F';
        case 3
            nombre = 'T';
    end
    
    grayback = imread([happypath,nombre,'.jpg']);
    dims = size(grayback);
    axis equal
    axis tight
    
    for uu = 1:dims(1)
        for jj = 1:dims(2)
            if grayback(uu,jj,3) <= 35 && grayback(uu,jj,2) <= 35 && grayback(uu,jj,1) <= 35
                grayback(uu,jj,1) = 100;
                grayback(uu,jj,2) = 100;
                grayback(uu,jj,3) = 100;
            end
        end
    end
    imwrite(grayback,[happypathout,nombre,'.jpg'],'jpg')
end

%% %
% pathToImages = 'C:\Ana\Julians_task\disk_masks';
% d = dir(fullfile(pathToImages, imageNames));
% nImages = length(d);
%
% for iImage = 1:nImages
%     [im, map] = imread(fullfile(pathToImages, d(iImage).name),'bmp');
%     if dims(im) < 3
%         %convert to color
%         img = ind2rgb(im,map);
%     end
%
% imwrite(im,fullfile(pathToImages,filename));
% end

% %%
%
% blackback = imread('1.jpg');
% % image(blackback)
% axis equal
% axis tight
%
% % imshow('1.jpg')
%
% % dims = size(blackback) % "y", "x" y "z" (en ese orden)
%
%
% [rows, columns, numberOfColorChannels] = size(blackback)
% % Extract the individual red, green, and blue color channels.
%
% redChannel = blackback(:, :, 1);
% greenChannel = blackback(:, :, 2);
% blueChannel = blackback(:, :, 3);
%
% if numberOfColorChannels == 3 || isempty(map)
%     %RGB or standard gray scale
%     blackPixels = blackback == 0;
%     blackback(blackPixels) = newGrayLevel;
%     %Write out
%     imwrite(blackback,'zaza1.jpg');
% else
%     grayback = map;
%     for gl = 1:size(map,1)
%         if map(gl,1) == 0 && map(gl, 2) == 0 && map(gl,3) == 0
%             %Change from white to our gray level
%             grayback(gl,1) = newGrayLevel;
%             grayback(gl,2) = newGrayLevel;
%             grayback(gl,3) = newGrayLevel;
%         end
%     end
%     %Apply new colormap
%     imwrite(blackback, grayback, 'graybackmaster.jpg');
% end
%%
% blackback = imread('1.jpg');
% image(blackback)
% axis equal
% axis tight
%
% % imshow('1.jpg')
%
% dims = size(blackback) % "y", "x" y "z" (en ese orden)
%
% newImage = ~all(blackback == 0, 3);