function [bodypart]= importDeepLabCutfile_VV_freebeh(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [bodypart1, bodypart2, bodypart3, bodypart4, bodypart5, bodypart6]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [bodypart1, bodypart2, bodypart3, bodypart4, bodypart5, bodypart6]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [bodypart1, bodypart2, bodypart3, bodypart4, bodypart5, bodypart6]
%   =
%   importfile('yourdeeplabcutoutput_DeepCut_resnet50_ReachDec10shuffle1_1030000.csv',4, 1169);
%

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 4;
    endRow = inf;
end

%% Format string for each line of text

% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
bodypart(1,:) = dataArray{:, 1};
bodypart(2,:) = dataArray{:, 2};
bodypart(3,:) = dataArray{:, 3};
bodypart(4,:) = dataArray{:, 4};
bodypart(5,:) = dataArray{:, 5};
bodypart(6,:) = dataArray{:, 6};
bodypart(7,:) = dataArray{:, 7};
bodypart(8,:) = dataArray{:, 8};
bodypart(9,:) = dataArray{:, 9};
bodypart(10,:) = dataArray{:, 10};
bodypart(11,:) = dataArray{:, 11};
bodypart(12,:) = dataArray{:, 12};
bodypart(13,:) = dataArray{:, 13};
bodypart(14,:) = dataArray{:, 14};
bodypart(15,:) = dataArray{:, 15};
bodypart(16,:) = dataArray{:, 16};
bodypart(17,:) = dataArray{:, 17};
bodypart(18,:) = dataArray{:, 18};
bodypart(19,:) = dataArray{:, 19};
bodypart(20,:) = dataArray{:, 20};
bodypart(21,:) = dataArray{:, 21};
bodypart(22,:) = dataArray{:, 22};
bodypart(23,:) = dataArray{:, 23};
bodypart(24,:) = dataArray{:, 24};
bodypart(25,:) = dataArray{:, 25};
bodypart(26,:) = dataArray{:, 26};
bodypart(27,:) = dataArray{:, 27};
bodypart(28,:) = dataArray{:, 28};
bodypart(29,:) = dataArray{:, 29};
bodypart(30,:) = dataArray{:, 30};
bodypart(31,:) = dataArray{:, 31};
bodypart(32,:) = dataArray{:, 32};
bodypart(33,:) = dataArray{:, 33};
bodypart(34,:) = dataArray{:, 34};
bodypart(35,:) = dataArray{:, 35};
bodypart(36,:) = dataArray{:, 36};
bodypart(37,:) = dataArray{:, 37};
bodypart(38,:) = dataArray{:, 38};
bodypart(39,:) = dataArray{:, 39};
bodypart(40,:) = dataArray{:, 40};
bodypart(41,:) = dataArray{:, 41};
bodypart(42,:) = dataArray{:, 42};
bodypart(43,:) = dataArray{:, 43};
bodypart(44,:) = dataArray{:, 44};
bodypart(45,:) = dataArray{:, 45};
bodypart(46,:) = dataArray{:, 46};
bodypart(47,:) = dataArray{:, 47};
bodypart(48,:) = dataArray{:, 48};
bodypart(49,:) = dataArray{:, 49};
bodypart(50,:) = dataArray{:, 50};
bodypart(51,:) = dataArray{:, 51};
bodypart(52,:) = dataArray{:, 52};
bodypart(53,:) = dataArray{:, 53};
bodypart(54,:) = dataArray{:, 54};
bodypart(55,:) = dataArray{:, 55};
bodypart(56,:) = dataArray{:, 56};
bodypart(57,:) = dataArray{:, 57};
bodypart(58,:) = dataArray{:, 58};
bodypart(59,:) = dataArray{:, 59};
bodypart(60,:) = dataArray{:, 60};
bodypart(61,:) = dataArray{:, 61};

