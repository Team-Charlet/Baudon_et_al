% This code is used to plot heatmap, Hypnogram(indicate behavior state),
% and normalized dFF across time.
% In foot schock test, Hypnogram: 4=(CS+ with FS) 2=(CS- without FS) 1=(interval)
% Input： CalciumTraces_Clean_AllSessions is a array，containing matrixs of
% normalized dFF (10 HZ) of each animal （frames x neuron numbers）
% In animal Matrixs,  each colum indicate the time course of dFF changes of
% one neuron. 


% Downsampling the calcium traces from 10 HZ to 1s
% Initialize the cell array for downsampled data
CalciumTraces_Clean_AllSessions_downsampled = cell(size(CalciumTraces_Clean_AllSessions));

% Loop through each matrix (each animal)
for i = 1:length(CalciumTraces_Clean_AllSessions)
    % Get the current matrix
    currentMatrix = CalciumTraces_Clean_AllSessions{i};

    % Get the number of cells (columns) and time points (rows)
    [numTimePoints, numCells] = size(currentMatrix);

    % Calculate the number of downsampled time points
    numDownsampledTimePoints = floor(numTimePoints / 10);

    % Initialize the downsampled matrix
    downsampledMatrix = zeros(numDownsampledTimePoints, numCells);

    % Loop through each cell (column)
    for j = 1:numCells
        % Downsample the column by averaging every 10 values
        for k = 1:numDownsampledTimePoints
            downsampledMatrix(k, j) = mean(currentMatrix((k-1)*10+1:k*10, j));
        end
    end

    % Store the downsampled matrix
    CalciumTraces_Clean_AllSessions_downsampled{i} = downsampledMatrix;
end

% downsampledData now contains the downsampled datasets

%% Combine and order the dataset
%Combine all the animals together in one dataset, concatenated by columns
% Each column represents a cell from any animal, and rows represent time points
% Initialize an empty matrix for concatenating data
combinedMatrix = [];

% Iterate through each animal's data and concatenate horizontally
for i = 1:length(CalciumTraces_Clean_AllSessions_downsampled)
    combinedMatrix = [combinedMatrix, CalciumTraces_Clean_AllSessions_downsampled{i}];
end

CalciumTraces_Clean_AllSessions_downsampled = combinedMatrix;

%order the data set by average activity

% Define the time range for analysis (in seconds)
startTimeSec = 590; % Modify as needed, start time in seconds
endTimeSec = 595; % Modify as needed, end time in seconds

% Convert time range to row indices (assuming 1 second per row)
startIndex = startTimeSec + 1;
endIndex = min(endTimeSec, size(CalciumTraces_Clean_AllSessions_downsampled, 1));

% Extract the data for the desired time range
dataInRange = CalciumTraces_Clean_AllSessions_downsampled(startIndex:endIndex, :);

% Calculate the average activity for each cell in the time range
averageActivityInRange = mean(dataInRange, 1);

% Sort the cells by their average activity in the range
% Sort the cells by their average activity
[~, sortedIndices] = sort(averageActivityInRange, 'ascend'); % 'descend' for highest to lowest

% Reorder the data matrix based on the sorted indices
sortedData = CalciumTraces_Clean_AllSessions_downsampled(:, sortedIndices);

CalciumTraces_Clean_AllSessions_downsampled_SortedByAverage = sortedData;

%% Plot the heatmap for all cell

%  Data is a combined matrix with cells in columns and time in rows

% Define the time range for the heatmap (in seconds)
startTimeSec = 0; % Modify as needed, e.g., start time in seconds
endTimeSec = 900; % Modify as needed, e.g., end time in seconds

% Convert time range to row indices (assuming 1 second per row)
startIndex = startTimeSec + 1;
endIndex = min(endTimeSec, size(CalciumTraces_Clean_AllSessions_downsampled_SortedByAverage, 1));

% Extract the data for the desired time range
dataToPlot = CalciumTraces_Clean_AllSessions_downsampled_SortedByAverage(startIndex:endIndex, :);
hypnoToPlot = Hypnogram(startIndex:endIndex);


% Plotting the heatmap
figure;
imagesc(dataToPlot');
colormap jet; % You can change the colormap if needed
colorbar;
xlabel('Time (seconds)');
ylabel('Cell ID');
caxis([0, 30]);

% Adjust x-axis to show time in seconds
xticks = get(gca, 'XTick');
set(gca, 'XTickLabel', round((xticks + startIndex - 1), 2)); % Adjust labels to reflect actual time

title('Heatmap of Cell Activity Over Time');



%plot heatmap and hypnogram
figure;
% Subplot for the heatmap
ax1 = subplot(2, 1, 1); % Two rows, one column, first plot
imagesc(dataToPlot');
colormap jet;
cb = colorbar;
ylabel('Cell ID');
title('Heatmap of Cell Activity');
caxis([0, 30]);

% Adjust x-axis to show time in seconds
xticks = get(ax1, 'XTick');
set(ax1, 'XTickLabel', round((xticks + startIndex - 1), 2));

% Subplot for the Hypnogram
ax2 = subplot(2, 1, 2); % Two rows, one column, second plot
plot(hypnoToPlot);
xlabel('Time (seconds)');
ylabel('Hypnogram Value');
title('Hypnogram');

% Manually adjust the width of the first subplot to match the second
pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
pos1(3) = pos2(3); % Set the width of ax1 equal to ax2
set(ax1, 'Position', pos1);

% Ensure the x-axis of both plots aligns properly
linkaxes([ax1, ax2], 'x');
% Adjust x-axis limits to the time range
xlim(ax1, [1 length(hypnoToPlot)]);
xlim(ax2, [1 length(hypnoToPlot)]);

%% DFF visualization
space = 30; % space between dataArrayEnv traces
colbar = [1,5];   % min and max color data
fs = 10;

fprintf('PLOT FIGURES - not saved\n')

columns = size(CalciumTraces_Clean_AllSessions_downsampled_SortedByAverage,2);
tmp = CalciumTraces_Clean_AllSessions_downsampled_SortedByAverage;

% plot
figure
for n = 1:columns
    plot(tmp(:,n)+(n-1)*space,'k')
    hold on
end
xlim([0,size(CalciumTraces_Clean_AllSessions_downsampled_SortedByAverage,1)]); ylim([-20,columns*space+20])
yt = 0:space*2:columns*space;
if isempty(yt)
    warning('yt vector is empty. Check the values of "columns" and "space".');
else
    yticks(yt)
end

tmp2 = yt/5; tmp2(1) = 1; yticklabels(tmp2)
xticks(0:60*fs:size(CalciumTraces_Clean_AllSessions_downsampled_SortedByAverage,1))
xticklabels(0:size(CalciumTraces_Clean_AllSessions_downsampled_SortedByAverage,1)/fs/60)
xlabel('Time [min]'); ylabel('Neuron ID')
box off