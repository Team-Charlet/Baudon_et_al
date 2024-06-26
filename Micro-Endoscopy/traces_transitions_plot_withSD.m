
n_mice = numel(Mouse_Names);
n_TotalSessions = numel(Hypnogram_AllSessions);
Transition_Duration = 20; % [s]
Transition_Duration = Transition_Duration.*Opts.General.FrameRate;
SD_Session = 10;
% Compute for each mouse separately the average Calcium Trace around
% transitions.
%Adapt from sleep analysis script: 
% Hypnogram 1=Awake = Interval in both FC and stress Exp, 
% Hypnogram 2=NREM  = CS- =lift by tail; 
% Hypnogram 4=REM   = CS+ = shake cage 

for i_mouse = 1:n_mice
    % Get Current Mouse
    Mouse_Name_Current_Mouse = Mouse_Names{i_mouse};
    i_current_session = 0;
    for i_session = 1:n_TotalSessions
        if strcmpi(Hypnogram_AllSessions(i_session).MouseName, Mouse_Name_Current_Mouse)
            i_current_session = i_current_session +1;
            Hypnogram_Current_Mouse(i_current_session) = Hypnogram_AllSessions(i_session);
            CalciumTraces_Current_Mouse{i_current_session} = CalciumTraces_Clean_AllSessions{i_session};
        end
    end
    
    % Consider each Session Separately
    n_sessions = numel(Hypnogram_Current_Mouse);
    for i_session = 1:n_sessions
        Hypnogram_Current_Session = Hypnogram_Current_Mouse(i_session);
        CalciumTraces_Current_Session = CalciumTraces_Current_Mouse{i_session};
        
        % Compute average trace
        CalciumTraces_Mean_Current_Session = nanmean(CalciumTraces_Current_Session, 2);
        RecordingDuration = numel(CalciumTraces_Mean_Current_Session);
        
        % Consider only stable states.
        n_states = numel(Hypnogram_Current_Session.StateDuration);
        StableState_Index = [];
        for i_state = 2:n_states
            if Hypnogram_Current_Session.StateDuration(i_state) >= (Opts.General.MinStableStateDuration)*(Opts.General.FrameRate)
                StableState_Index = [StableState_Index, i_state];
            end
        end
        tmp_Transitions_Array = [Hypnogram_Current_Session.StateChanges];
        tmp_StatesTag_Array = [Hypnogram_Current_Session.StateType];
        StableState_Start_Array = tmp_Transitions_Array(StableState_Index);
        StableState_Tag_Array = tmp_StatesTag_Array(StableState_Index);
        n_transitions = numel(StableState_Start_Array);
        StableStatePrevious_Index_Array = NaN(n_transitions, 1);
        StableStatePrevious_Tag_Array = NaN(n_transitions, 1);
        ToRemove = [];
        for i_state = 1:n_transitions
            if StableState_Start_Array(i_state) <= Transition_Duration + 1 || StableState_Start_Array(i_state) >= RecordingDuration - (Transition_Duration - 1)
                ToRemove = [ToRemove, i_state];
                continue
            end
            if StableState_Index(i_state) > 1
                StableStatePrevious_Index_Array(i_state) = StableState_Index(i_state) - 1;
                StableStatePrevious_Tag_Array(i_state) = tmp_StatesTag_Array(StableState_Index(i_state) - 1);
            end
        end
        % Ignore states that come too early or late in the recordings.
        StableState_Index(ToRemove) = [];
        StableState_Start_Array(ToRemove) = [];
        StableState_Tag_Array(ToRemove) = [];
        StableStatePrevious_Index_Array(ToRemove) = [];
        StableStatePrevious_Tag_Array(ToRemove) = [];
        n_transitions = numel(StableState_Start_Array);
        
        % Get the Traces at each transition
        Transition_Segments_tmp = NaN(n_transitions, (1 + (Transition_Duration*2)));
        Transition_Segments.Awake2NREM = Transition_Segments_tmp;
        Transition_Segments.NREM2REM = Transition_Segments_tmp;
        Transition_Segments.NREM2Awake = Transition_Segments_tmp;
        Transition_Segments.REM2Awake = Transition_Segments_tmp;
        
        for i_state_start = 1:n_transitions
            try
                tmp = CalciumTraces_Mean_Current_Session((StableState_Start_Array(i_state_start) - Transition_Duration):(StableState_Start_Array(i_state_start) + Transition_Duration));
                tmp = tmp - nanmin(tmp); % Subtract the minimum for more comparable data.
                tmp = tmp./nanmax(tmp); % Normalize for more comparable data
                if StableState_Tag_Array(i_state_start) == Opts.General.TagAWAKE
                    if StableStatePrevious_Tag_Array(i_state_start) == 2
                        Transition_Segments.NREM2Awake(i_state_start, 1:(1 + (Transition_Duration*2))) = tmp;
                    elseif StableStatePrevious_Tag_Array(i_state_start) == 4
                        Transition_Segments.REM2Awake(i_state_start, 1:(1 + (Transition_Duration*2))) = tmp;
                    end
                elseif StableState_Tag_Array(i_state_start) == Opts.General.TagNoNREM
                        Transition_Segments.Awake2NREM(i_state_start, 1:(1 + (Transition_Duration*2))) = tmp;
                elseif StableState_Tag_Array(i_state_start) == Opts.General.TagREM
                    Transition_Segments.NREM2REM(i_state_start, 1:(1 + (Transition_Duration*2))) = tmp;
                end
            catch
                keyboard
            end
        end
        
        % Mean over all transitions
        Transition_Segments_Mean_perSession(i_session).Awake2NREM = nanmean(Transition_Segments.Awake2NREM, 1);
        Transition_Segments_Mean_perSession(i_session).NREM2REM = nanmean(Transition_Segments.NREM2REM, 1);
        Transition_Segments_Mean_perSession(i_session).NREM2Awake = nanmean(Transition_Segments.NREM2Awake, 1);
        Transition_Segments_Mean_perSession(i_session).REM2Awake = nanmean(Transition_Segments.REM2Awake, 1);       
        
    end
    Transition_Segments_Mean_perMouse{i_mouse} = Transition_Segments_Mean_perSession;
    
    
    clear Hypnogram_Current_Mouse
    clear CalciumTraces_Current_Mouse
end


% Mean everything per mouse
AllTracesMatrix_Awake2NREM = NaN(n_sessions, (2*Transition_Duration)+1, n_mice);
AllTracesMatrix_NREM2REM = NaN(n_sessions, (2*Transition_Duration)+1, n_mice);
AllTracesMatrix_NREM2Awake = NaN(n_sessions, (2*Transition_Duration)+1, n_mice);
AllTracesMatrix_REM2Awake = NaN(n_sessions, (2*Transition_Duration)+1, n_mice);

for i_session = 1:n_sessions
    for i_mouse = 1:n_mice
        current_mouse_Transitions = Transition_Segments_Mean_perMouse{i_mouse};
        tmp_transitions_Awake2NREM = cell(n_sessions, 1);
        tmp_transitions_NREM2REM = cell(n_sessions, 1);
        tmp_transitions_NREM2Awake = cell(n_sessions, 1);
        tmp_transitions_REM2Awake = cell(n_sessions, 1);
        
        for i = 1:n_sessions
            tmp_transitions_Awake2NREM{i} = current_mouse_Transitions(i).Awake2NREM;
            tmp_transitions_NREM2REM{i} = current_mouse_Transitions(i).NREM2REM;
            tmp_transitions_NREM2Awake{i} = current_mouse_Transitions(i).NREM2Awake;
            tmp_transitions_REM2Awake{i} = current_mouse_Transitions(i).REM2Awake;
        end
        tmp_transitions_Awake2NREM = cell2mat(tmp_transitions_Awake2NREM);
        tmp_transitions_NREM2REM = cell2mat(tmp_transitions_NREM2REM);
        tmp_transitions_NREM2Awake = cell2mat(tmp_transitions_NREM2Awake);
        tmp_transitions_REM2Awake = cell2mat(tmp_transitions_REM2Awake);
        AllTracesMatrix_Awake2NREM(:,:, i_mouse) = tmp_transitions_Awake2NREM;
        AllTracesMatrix_NREM2REM(:,:, i_mouse) = tmp_transitions_NREM2REM;
        AllTracesMatrix_NREM2Awake(:,:, i_mouse) = tmp_transitions_NREM2Awake;
        AllTracesMatrix_REM2Awake(:,:, i_mouse) = tmp_transitions_REM2Awake;
    end
end
Traces_MeansPerMice_Awake2NREM = nanmean(AllTracesMatrix_Awake2NREM, 3);
Traces_MeansPerMice_NREM2REM = nanmean(AllTracesMatrix_NREM2REM, 3);
Traces_MeansPerMice_NREM2Awake = nanmean(AllTracesMatrix_NREM2Awake, 3);
Traces_MeansPerMice_REM2Awake = nanmean(AllTracesMatrix_REM2Awake, 3);

% Calculate SD and SEM for each transition type
SD_PerMice_Awake2NREM = nanstd(AllTracesMatrix_Awake2NREM, [], 3);
SEM_PerMice_Awake2NREM = SD_PerMice_Awake2NREM / sqrt(n_mice);

SD_PerMice_NREM2REM = nanstd(AllTracesMatrix_NREM2REM, [], 3);
SEM_PerMice_NREM2REM = SD_PerMice_NREM2REM / sqrt(n_mice);

SD_PerMice_NREM2Awake = nanstd(AllTracesMatrix_NREM2Awake, [], 3);
SEM_PerMice_NREM2Awake = SD_PerMice_NREM2Awake / sqrt(n_mice);

SD_PerMice_REM2Awake = nanstd(AllTracesMatrix_REM2Awake, [], 3);
SEM_PerMice_REM2Awake = SD_PerMice_REM2Awake / sqrt(n_mice);

% Put minimum to zero.
for i_session = 1:n_sessions
    Traces_MeansPerMice_Awake2NREM(i_session, :) = Traces_MeansPerMice_Awake2NREM(i_session, :) - nanmin(Traces_MeansPerMice_Awake2NREM(i_session, :));
    Traces_MeansPerMice_NREM2REM(i_session, :) = Traces_MeansPerMice_NREM2REM(i_session, :) - nanmin(Traces_MeansPerMice_NREM2REM(i_session, :));
    Traces_MeansPerMice_NREM2Awake(i_session, :) = Traces_MeansPerMice_NREM2Awake(i_session, :) - nanmin(Traces_MeansPerMice_NREM2Awake(i_session, :));
    Traces_MeansPerMice_REM2Awake(i_session, :) = Traces_MeansPerMice_REM2Awake(i_session, :) - nanmin(Traces_MeansPerMice_REM2Awake(i_session, :));
end

grandavg_Awake2NREM = smooth(nanmean(Traces_MeansPerMice_Awake2NREM, 1));
grandavg_NREM2REM = smooth(nanmean(Traces_MeansPerMice_NREM2REM, 1));
grandavg_NREM2Awake = smooth(nanmean(Traces_MeansPerMice_NREM2Awake, 1));
grandavg_REM2Awake = smooth(nanmean(Traces_MeansPerMice_REM2Awake, 1));

grandavg_Awake2NREM = grandavg_Awake2NREM(:)';
grandavg_NREM2REM = grandavg_NREM2REM(:)';
grandavg_NREM2Awake = grandavg_NREM2Awake(:)';
grandavg_REM2Awake = grandavg_REM2Awake(:)';


% Apply zero-basing and smoothing to SEM
SEM_PerMice_Awake2NREM_Processed = smooth(nanmean(SD_PerMice_Awake2NREM, 1) / sqrt(n_mice)) - nanmin(smooth(nanmean(SD_PerMice_Awake2NREM, 1) / sqrt(n_mice)));
SEM_PerMice_NREM2REM_Processed = smooth(nanmean(SD_PerMice_NREM2REM, 1) / sqrt(n_mice)) - nanmin(smooth(nanmean(SD_PerMice_NREM2REM, 1) / sqrt(n_mice)));
SEM_PerMice_NREM2Awake_Processed = smooth(nanmean(SD_PerMice_NREM2Awake, 1) / sqrt(n_mice)) - nanmin(smooth(nanmean(SD_PerMice_NREM2Awake, 1) / sqrt(n_mice)));
SEM_PerMice_REM2Awake_Processed = smooth(nanmean(SD_PerMice_REM2Awake, 1) / sqrt(n_mice)) - nanmin(smooth(nanmean(SD_PerMice_REM2Awake, 1) / sqrt(n_mice)));

SEM_PerMice_Awake2NREM_Processed = SEM_PerMice_Awake2NREM_Processed(:)';
SEM_PerMice_NREM2REM_Processed = SEM_PerMice_NREM2REM_Processed(:)';
SEM_PerMice_NREM2Awake_Processed = SEM_PerMice_NREM2Awake_Processed(:)';
SEM_PerMice_REM2Awake_Processed = SEM_PerMice_REM2Awake_Processed(:)';

%% Plot with SEM as Shaded Area

% Ensure you have these standard deviation arrays calculated similarly to your mean arrays
% std_Awake2NREM, std_NREM2REM, std_NREM2Awake, std_REM2Awake

% Options and figure setup remain the same
% ...
% Options
Time_Array = (-Transition_Duration:Transition_Duration)./Opts.General.FrameRate; % [s], distance from transition


n_subplot_rows = 1;
n_subplot_columns = 4;
% xticks_array = -20:10:20;
yticks_array = 0:0.2:1;
axis_limits = [-20, 20, 0, 1];
axis_FontSize = 16;
TitleFontSize = 16;
SupTitleFontSize = 26;
GrandMean_LineWidth = 2.5;
figure('units','normalized','outerposition',[0 0 1 1]);

% Int to CS- (Awake2NREM)
subplot(n_subplot_rows, n_subplot_columns, 1);
hold on;

% Plot SEM as shaded area
fill([Time_Array fliplr(Time_Array)], ...
     [grandavg_Awake2NREM + SEM_PerMice_Awake2NREM_Processed fliplr(grandavg_Awake2NREM - SEM_PerMice_Awake2NREM_Processed)], ...
     [0.9 0.9 0.9], 'LineStyle', 'none'); % Light grey color for SEM

% Plot the grand average line
plot(Time_Array, grandavg_Awake2NREM, 'r', 'LineWidth', GrandMean_LineWidth);

% Plotting individual sessions
for i_session = 1:n_sessions
    plot(Time_Array, Traces_MeansPerMice_Awake2NREM(i_session, :), 'k');
end

ax = gca;
axis square
box on
grid on
% xticks(xticks_array)
% yticks(yticks_array)
axis(axis_limits)
ax.FontSize = axis_FontSize;
title('Int to CS-', 'FontSize', TitleFontSize)

xlabel('Distance from State Transition [s]')
ylabel('Average \DeltaF/F Normalized')


% Int to CS+ (NREM2REM)
subplot(n_subplot_rows, n_subplot_columns, 2);

hold on;

% Plot SEM as shaded area
fill([Time_Array fliplr(Time_Array)], ...
     [grandavg_NREM2REM + SEM_PerMice_NREM2REM_Processed fliplr(grandavg_NREM2REM - SEM_PerMice_NREM2REM_Processed)], ...
     [0.9 0.9 0.9], 'LineStyle', 'none'); % Light grey color for SEM

% Plot the grand average line
plot(Time_Array, grandavg_NREM2REM, 'r', 'LineWidth', GrandMean_LineWidth);

for i_session = 1:n_sessions
    if i_session == SD_Session
        plot(Time_Array, Traces_MeansPerMice_NREM2REM(i_session, :), 'k')
    else
        plot(Time_Array, Traces_MeansPerMice_NREM2REM(i_session, :), 'k')
    end
end

ax = gca;
axis square
box on
grid on
% xticks(xticks_array)
% yticks(yticks_array)
axis(axis_limits)
ax.FontSize = axis_FontSize;
title('Int to CS+', 'FontSize', TitleFontSize)

% CS-toInt (NREM2Awake)

subplot(n_subplot_rows, n_subplot_columns, 3);
hold on;

% Plot SEM as shaded area
fill([Time_Array fliplr(Time_Array)], ...
     [grandavg_NREM2Awake + SEM_PerMice_NREM2Awake_Processed fliplr(grandavg_NREM2Awake - SEM_PerMice_NREM2Awake_Processed)], ...
     [0.9 0.9 0.9], 'LineStyle', 'none'); % Light grey color for SEM


for i_session = 1:n_sessions
    if i_session == SD_Session
        plot(Time_Array, Traces_MeansPerMice_NREM2Awake(i_session, :), 'k')
    else
        plot(Time_Array, Traces_MeansPerMice_NREM2Awake(i_session, :), 'k')
    end
end
plot(Time_Array, grandavg_NREM2Awake, 'r', 'LineWidth', GrandMean_LineWidth)

ax = gca;
axis square
box on
grid on
% xticks(xticks_array)
% yticks(yticks_array)
axis(axis_limits)
ax.FontSize = axis_FontSize;
title('CS- to Int', 'FontSize', TitleFontSize)

% CS+ to Int (REM2Awake)
subplot(n_subplot_rows, n_subplot_columns, 4);
hold on;

% Plot SEM as shaded area
fill([Time_Array fliplr(Time_Array)], ...
     [grandavg_REM2Awake + SEM_PerMice_REM2Awake_Processed fliplr(grandavg_REM2Awake - SEM_PerMice_REM2Awake_Processed)], ...
     [0.9 0.9 0.9], 'LineStyle', 'none'); % Light grey color for SEM

% Plot the grand average line
plot(Time_Array, grandavg_REM2Awake, 'r', 'LineWidth', GrandMean_LineWidth);


for i_session = 1:n_sessions
    if i_session == SD_Session
        plot(Time_Array, Traces_MeansPerMice_REM2Awake(i_session, :), 'k')
    else
        plot(Time_Array, Traces_MeansPerMice_REM2Awake(i_session, :), 'k')
    end
end

ax = gca;
axis square
box on
grid on
% xticks(xticks_array)
% yticks(yticks_array)
axis(axis_limits)
ax.FontSize = axis_FontSize;
title('CS+ to Int', 'FontSize', TitleFontSize)

suptitle_text = sprintf('%s\nTraces at Transitions for each Recording Session\nAverage over mice.', Opts.CellType);
try
    if verLessThan('matlab','9.5')
        h_suptitle = suptitle(suptitle_text);
        h_suptitle.FontSize = SupTitleFontSize;
        h_suptitle.FontWeight = 'bold';
    else
        h_suptitle = sgtitle(suptitle_text, 'FontSize', SupTitleFontSize, 'FontWeight', 'bold');
    end
catch
    warning ('Could not add suptitle.')
end

