%% ===== Section 1 - Windowing your data =====
clear; clc;

filePath = "C:\Users\callo\OneDrive\Desktop\University\2025-26\Wearbles\GrandChallengeData.mat";

S = load(filePath, "clean_labels");
labels = S.clean_labels;
clear S;

participants = fieldnames(labels);

% Storage
durations = [];
isFallAfter = [];

%% ---- Loop all participants and trials ----
for p = 1:numel(participants)

    pID = participants{p};
    trials = fieldnames(labels.(pID));

    for t = 1:numel(trials)

        trialName = trials{t};
        T = labels.(pID).(trialName);   % Mx3: [label, time_ms, group]

        if isempty(T) || size(T,2) < 3
            continue
        end

        %%Old code (Signal Processing Assignment)
        % grp = T(:,3);
        %time_s = T(:,2)/1000;

        % Find near-fall rows
        %nf_rows = find(grp == 1);

        %%New Code (Final Assignment)
        eventLabel = T(:,1);
        time_s = T(:,2)/1000;

        % Find near-fall start rows (label == 2)
        nf_start_rows = find(eventLabel == 2);
    
        if isempty(nf_start_rows)
            continue
        end

        % Split into contiguous blocks
        %breaks = [1; find(diff(nf_rows)>1)+1; numel(nf_rows)+1];

        %for b = 1:numel(breaks)-1

        %     rowsChunk = nf_rows(breaks(b):breaks(b+1)-1);
        % 
        %     % Use first and last time as period
        %     tStart = time_s(rowsChunk(1));
        %     tEnd   = time_s(rowsChunk(end));
        % 
        %     dur = tEnd - tStart;
        % 
        %     % Check if a fall occurs AFTER this near-fall
        %     fall_rows = find(grp == 4);
        % 
        %     fallAfter = false;
        %     if ~isempty(fall_rows)
        %         fallAfter = any(time_s(fall_rows) > tEnd);
        %     end
        % 
        %     durations(end+1,1) = dur;
        %     isFallAfter(end+1,1) = fallAfter;
        % 
        % end
        %%New Code (Final Assignment)
        for i = 1:numel(nf_start_rows)

            r_start = nf_start_rows(i);
            t_start = time_s(r_start);
        
            % Find all rows that come after this near-fall start
            future_rows = find(time_s > t_start);
        
            if isempty(future_rows)
                continue
            end
        
            % Find the first label that is 3 (recovery) or 4 (fall)
            future_labels = eventLabel(future_rows);
            term_idx = find(future_labels == 3 | future_labels == 4, 1, 'first');
        
            if isempty(term_idx)
                continue
            end
        
            r_end   = future_rows(term_idx);
            t_end   = time_s(r_end);
            termLbl = eventLabel(r_end);
        
            dur = t_end - t_start;
        
            % If it ended with a 4 (fall), mark it as fell
            fellAfter = (termLbl == 4);
        
            durations(end+1,1)   = dur;
            isFallAfter(end+1,1) = fellAfter;
        
        end
    end
end

%% =============================
% 1) Count number of total near-fall periods
%This counts the number of tracked near-fall periods from the cycled code
nTotal = numel(durations);

fprintf("\nTotal near-fall periods across all participants and trials: %d\n", nTotal);

%% =============================
% 2) Average & Standard Deviation

meanDur = mean(durations);
stdDur  = std(durations);

fprintf("Average duration: %.3f s\n", meanDur);
fprintf("Standard deviation: %.3f s\n", stdDur);

%% =============================
% 3) Median & IQR

medianDur = median(durations);

%Calculate the 25th and 75th percentiles
Q1 = quantile(durations, 0.25);
Q3 = quantile(durations,0.75);
IQR = Q3 - Q1;

%I can use the IQR function from matlab to check IQR manually
IQR_matlab = iqr(durations);

fprintf("25th percentile (Q1): %.3f s and 75th percentile (Q3): %.3f s\n", Q1, Q3);
fprintf("Median duration: %.3f s\n", medianDur);
fprintf("Interquartile range (IQR) is: %.3f s\n", IQR);

%% =============================
% 4) Compare near-falls prior to fall and near-falls with recovery

dur_fall = durations(isFallAfter==1);
dur_rec  = durations(isFallAfter==0);

mean_fall = mean(dur_fall);
mean_rec  = mean(dur_rec);

fprintf("\nMean duration (Near-Fall BEFORE Fall): %.3f s\n", mean_fall);
fprintf("Mean duration (Recovered Near-Fall): %.3f s\n", mean_rec);

% Check the two with a t-test 
[~,pValue] = ttest2(dur_fall, dur_rec);

fprintf("Two-sample t-test p-value: %.5f\n", pValue);


%% =============================
% 5) Participant Results Table (Recovered vs Fell) for assignment
pList = strings(0,1);
recCounts = [];
fallCounts = [];

for p = 1:numel(participants)

    pID = participants{p};
    trials = fieldnames(labels.(pID));
    nRec_p = 0;
    nFall_p = 0;

    for t = 1:numel(trials)

        trialName = trials{t};
        T = labels.(pID).(trialName);

        if isempty(T) || size(T,2) < 3
            continue
        end

        eventLabel = T(:,1);
        time_s = T(:,2)/1000;

        nf_start_rows = find(eventLabel == 2);

        if isempty(nf_start_rows)
            continue
        end

        for i = 1:numel(nf_start_rows)

            r_start = nf_start_rows(i);
            t_start = time_s(r_start);

            future_rows = find(time_s > t_start);
            if isempty(future_rows), continue; end

            future_labels = eventLabel(future_rows);
            term_idx = find(future_labels == 3 | future_labels == 4, 1, 'first');
            if isempty(term_idx), continue; end

            termLbl = future_labels(term_idx);

            if termLbl == 4
                nFall_p = nFall_p + 1;
            else
                nRec_p = nRec_p + 1;
            end

        end
    end

    pList(end+1,1) = string(pID);
    recCounts(end+1,1) = nRec_p;
    fallCounts(end+1,1) = nFall_p;

end

% Create table
resultsTable = table(pList, recCounts, fallCounts, recCounts+fallCounts, ...
    'VariableNames', {'Participant', ...
                      'NearFall_Recovered', ...
                      'NearFall_Fell', ...
                      'Total'});

% Add overall row
overallRow = table("OVERALL", sum(recCounts), sum(fallCounts), ...
    sum(recCounts)+sum(fallCounts), ...
    'VariableNames', {'Participant', ...
                      'NearFall_Recovered', ...
                      'NearFall_Fell', ...
                      'Total'});

resultsTable = [resultsTable; overallRow];

fprintf("\n=============================\n");
disp(resultsTable);
fprintf("=============================\n");
%% Near-fall duration distribution

% Example: your durations vector should already exist
% durations = [...];

% Summary statistics
mean_val = mean(durations);
median_val = median(durations);

Q1 = prctile(durations,25);
Q3 = prctile(durations,75);

% Plot histogram
figure
histogram(durations,25,'FaceColor',[0.3 0.3 0.3],'EdgeColor','none')
hold on

% Plot vertical lines
xline(mean_val,'r','LineWidth',1,'DisplayName','Mean')

xline(median_val,'b','LineWidth',1,'DisplayName','Median')

xline(Q1,'k--','LineWidth',1,'DisplayName','Q1')

xline(Q3,'k--','LineWidth',1,'DisplayName','Q3')

xlabel('Near-fall Duration (s)')
ylabel('Number of Events')
title('Distribution of Near-Fall Durations')

legend('Near-fall durations','Mean','Median','Q1','Q3')

grid on
hold off

%% Recovery near-falls

mean_rec = mean(dur_rec);
median_rec = median(dur_rec);
Q1_rec = prctile(dur_rec,25);
Q3_rec = prctile(dur_rec,75);

figure
histogram(dur_rec,20,'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')
hold on

xline(mean_rec,'r','LineWidth',2,'DisplayName','Mean')
xline(median_rec,'b','LineWidth',2,'DisplayName','Median')
xline(Q1_rec,'k--','LineWidth',2,'DisplayName','Q1')
xline(Q3_rec,'k--','LineWidth',2,'DisplayName','Q3')

xlabel('Near-fall Duration (s)')
ylabel('Number of Events')
title('Distribution of Near-fall Durations (Recovery)')

legend('Near-fall durations','Mean','Median','Q1','Q3','Location','best')

grid on
hold off

%% Near-falls that transition to fall

mean_fall = mean(dur_fall);
median_fall = median(dur_fall);
Q1_fall = prctile(dur_fall,25);
Q3_fall = prctile(dur_fall,75);

figure
histogram(dur_fall,20,'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')
hold on

xline(mean_fall,'r','LineWidth',2,'DisplayName','Mean')
xline(median_fall,'b','LineWidth',2,'DisplayName','Median')
xline(Q1_fall,'k--','LineWidth',2,'DisplayName','Q1')
xline(Q3_fall,'k--','LineWidth',2,'DisplayName','Q3')

xlabel('Near-fall Duration (s)')
ylabel('Number of Events')
title('Distribution of Near-fall Durations (Transition to Fall)')

legend('Near-fall durations','Mean','Median','Q1','Q3','Location','best')

grid on
hold off