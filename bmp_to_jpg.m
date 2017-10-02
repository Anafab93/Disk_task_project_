
% % Demo to create JPG copies of BMP files in a different folder.
% % By Image Analyst
% % inputFolder = fileparts(which('cameraman.tif')) % Determine where demo folder is (works with all versions).
% inputFolder = fullfile(pwd, 'C:\Ana\Julians_task\disk_masks'); 
% filePattern = fullfile(inputFolder, '*.bmp')
% % Get list of all BMP files in input folder
% bmpFiles = dir(filePattern)
% % Create the output folder:
 outputFolder = fullfile(pwd, 'C:\Ana\Julians_task\new_disk_masks_jpg'); 
% mkdir(outputFolder);
% 
% % if ~exist(outputFolder, 'dir')
% % 	mkdir(outputFolder);
% % end
% figure;

happypath = 'C:\Ana\Julians_task\disk_masks\';
happypathout = 'C:\Ana\Julians_task\scripts\faces\';

for ii = 1:300
    m = imread([happypath,num2str(ii),'.bmp']);
    imwrite(m,[happypathout,num2str(ii),'.jpg'],'jpg')
%     saveas(gcf,[happypathout,num2str(ii),'.jpg'])
end

% % Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % Loop over all bmp files, making a jpg version 
% % of them in the output folder.
% 
% for k = 1 : length(bmpFiles)
% 	% Read in .bmp file
% 	baseFileName = bmpFiles(k).name;
% 	fullFileNameInput = fullfile(inputFolder, baseFileName);
% 	rgbImage = imread(fullFileNameInput);
% 	subplot(1, 2, 1);
% 	imshow(rgbImage);
% 	title('Original image', 'FontSize', 30);
% 	drawnow;
% 	
% 	% Prepare output file name
% 	fullFileNameOutput = fullfile(outputFolder, baseFileName);
% 	% Converts to JPEG and gives it the .jpg extension
% 	fullFileNameOutput = strrep(lower(fullFileNameOutput), '.bmp', 'jpg')
% 	imwrite(rgbImage,fullFileNameOutput);
% end
% 
% winopen(outputFolder);
