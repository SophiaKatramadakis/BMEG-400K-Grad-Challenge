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

        grp = T(:,3);
        time_s = T(:,2)/1000;

        % Find near-fall rows
        nf_rows = find(grp == 1);

        if isempty(nf_rows)
            continue
        end

        % Split into contiguous blocks
        breaks = [1; find(diff(nf_rows)>1)+1; numel(nf_rows)+1];

        for b = 1:numel(breaks)-1

            rowsChunk = nf_rows(breaks(b):breaks(b+1)-1);

            % Use first and last time as period
            tStart = time_s(rowsChunk(1));
            tEnd   = time_s(rowsChunk(end));

            dur = tEnd - tStart;

            % Check if a fall occurs AFTER this near-fall
            fall_rows = find(grp == 4);

            fallAfter = false;
            if ~isempty(fall_rows)
                fallAfter = any(time_s(fall_rows) > tEnd);
            end

            durations(end+1,1) = dur;
            isFallAfter(end+1,1) = fallAfter;

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