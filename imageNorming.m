clear

% reading in images into cell array, with each cell being an image 
imgMatrix{1} = imread('10 deg flat flat 1:5 2*.jpg');
imgMatrix{2} = imread('10 deg flat rough 1:4*.jpg');
imgMatrix{3} = imread('10 deg rough flat 1:5*.jpg');
imgMatrix{4} = imread('10 deg rough rough 1:4*.jpg');
imgMatrix{5} = imread('20 deg flat flat 1:5 redo*.jpg');
imgMatrix{6} = imread('20 deg flat rough 1:5 redo*.jpg');
imgMatrix{7} = imread('20 deg rough flat 1:4 redo*.jpg');
imgMatrix{8} = imread('20 deg rough rough 0''4" redo*.jpg');
imgMatrix{9} = imread('30 deg flat flat 1:4*.jpg');
imgMatrix{10} = imread('30 deg flat rough 0''6"*.jpg');
imgMatrix{11} = imread('30 deg rough flat 0''4"*.jpg');
imgMatrix{12} = imread('30 deg rough rough 0''8"*.jpg');

% angles of images 
angleMatrix = [10, 10, 10, 10, 20, 20, 20, 20, 30, 30, 30, 30];

% integration times 
intTimeMatrix = [1/5, 1/4, 1/5, 1/4, 1/5, 1/5, 1/4, 0.4, 1/4, 0.6, 0.4, 0.8];
figure('Units','normalized','Position',[0 0 0.8 0.5333])

tiledlayout(3, 4)
    for ii=1:1:size(imgMatrix,2)
        nexttile(ii);
        imagesc(imgMatrix{ii})
        axis off        
    end
    set(gcf, 'Color', 'white')


% converting image data to double format, normalizing by 255 (to be in the
% range of [0 1] for using imagesc with doubles and recalculating
% intensities of each image for smallest integration time in matrix 'intTimeMatrix'

for ii=1:1:size(imgMatrix,2)
     doubleMatrix{ii} = double(imgMatrix{ii})/255; % turn into double, so we don't lose information in the multiplication/division
     doubleMatrixSameTimeMin{ii} = doubleMatrix{ii}*min(intTimeMatrix)/intTimeMatrix(ii);
end



% multiplier matrix to scale intensities by a multiple of smallest
% integration time 
multiplierMatrix = [1, 2, 3, 5, 1, 2, 3, 5, 1, 2, 3, 5];

% plotting data number in images shows multiplier for intensity to show
% data with multiple of minimum integration time given by entry in
% 'multiplierMatrix'
figure('Units','normalized','Position',[0 0 0.8 0.5333])
tiledlayout(3, 4)
    for ii=1:1:size(doubleMatrixSameTimeMin,2)
        nexttile(ii);
        imagesc(doubleMatrixSameTimeMin{ii}.*multiplierMatrix(ii))
        dim = [0.228 + 0.22*mod(ii-1,4) 0.915 - 0.295*floor((ii-1)/4) 0 0];
        str = num2str(multiplierMatrix(ii)) + "x";
        annotation("textbox", dim,'String',str,'FitBoxToText','on', 'Fontsize', 16, 'Color', [1 1 1], 'EdgeColor', 'none');
        axis off        
    end
    annotation("textbox", [0.05 0.835 0 0],'String','10º','FitBoxToText','on', 'Fontsize', 16, 'Color', [0 0 0], 'EdgeColor', 'none');
    annotation("textbox", [0.05 0.54 0 0],'String','20º','FitBoxToText','on', 'Fontsize', 16, 'Color', [0 0 0], 'EdgeColor', 'none');
    annotation("textbox", [0.05 0.245 0 0],'String','30º','FitBoxToText','on', 'Fontsize', 16, 'Color', [0 0 0], 'EdgeColor', 'none');
    
    annotation("textbox", [0.145 0.98 0 0],'String','flat-flat','FitBoxToText','on', 'Fontsize', 16, 'Color', [0 0 0], 'EdgeColor', 'none');
    annotation("textbox", [0.36 0.98 0 0],'String','flat-rough','FitBoxToText','on', 'Fontsize', 16, 'Color', [0 0 0], 'EdgeColor', 'none');
    annotation("textbox", [0.58 0.98 0 0],'String','rough-flat','FitBoxToText','on', 'Fontsize', 16, 'Color', [0 0 0], 'EdgeColor', 'none');
    annotation("textbox", [0.79 0.98 0 0],'String','rough-rough','FitBoxToText','on', 'Fontsize', 16, 'Color', [0 0 0], 'EdgeColor', 'none');
    
    set(gcf, 'Color', 'white')




