%% ===== P1960: Events x (AccelX, GyroX) grid | overlay Back + Left Thigh + Right Thigh =====
clear; clc; close all;

%% Information to remember (assignment notation)
% IMU columns:
%   1 = time (s)
%   2–4 = accel X,Y,Z (m/s^2)
%   5–7 = gyro  X,Y,Z (deg/s)

%% 1) Load MAT file
filePath = "C:\Users\callo\OneDrive\Desktop\University\2025-26\Wearbles\GrandChallengeData.mat";
S = load(filePath);

disp("Top-level variables in MAT file:");
disp(fieldnames(S))

%% 2) Auto-find raw sensor struct
names = fieldnames(S);
raw = [];

pID = "P1960";

for k = 1:numel(names)
    nm = names{k};
    if strcmp(nm, "clean_labels"), continue; end
    if ~isstruct(S.(nm)), continue; end
    if ~isfield(S.(nm), pID), continue; end

    trialFields = fieldnames(S.(nm).(pID));
    if isempty(trialFields), continue; end

    t1 = trialFields{1};
    if isstruct(S.(nm).(pID).(t1))
        sf = fieldnames(S.(nm).(pID).(t1));
        if any(strcmp(sf,"Back")) && any(strcmp(sf,"Left_Thigh")) && any(strcmp(sf,"Right_Thigh"))
            raw = S.(nm);
            fprintf("Using RAW sensor struct: S.%s\n", nm);
            break
        end
    end
end

if isempty(raw)
    error("Could not auto-find the raw sensor struct. Check fieldnames(S).");
end

%% 3) Events (rows) - base names
eventList = ["Walk", "Sit", "Slip"];

%% 4) Plot settings
ds = 10;                % downsample for speed (set ds=1 for full resolution)
demeanSignals = true;   % center each channel (helps compare shapes)
tMax = 30;              % plot 30 seconds only

%% 5) Grid figure: rows=events, col1=AccelX overlay, col2=GyroX overlay
figure('Color','w');
tl = tiledlayout(numel(eventList), 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf("%s | First %.0f s | X-axis only: Back vs Left Thigh vs Right Thigh", pID, tMax));

for r = 1:numel(eventList)
    eventBase = eventList(r);

    % Pick the correct trial field (expects T#_Event naming)
    trialField = pickTrialField_Tprefix(raw, pID, eventBase);

    if strlength(trialField) == 0
        nexttile; axis off; text(0,0.5,"No matching trial for: " + eventBase, 'FontSize', 11);
        nexttile; axis off; text(0,0.5,"No matching trial for: " + eventBase, 'FontSize', 11);
        continue;
    end

    fprintf("Event '%s' -> using trial field '%s'\n", eventBase, trialField);

    X = raw.(pID).(trialField);

    Back = X.Back;
    LTh  = X.Left_Thigh;
    RTh  = X.Right_Thigh;

    % Downsample (plotting only)
    B = Back(1:ds:end,:);
    L = LTh(1:ds:end,:);
    R = RTh(1:ds:end,:);

    % ---- Crop to first 30 seconds of THIS trial ----
    tStart = B(1,1);                   % trial starts at this time (seconds)
    tEnd   = tStart + tMax;            % first 30 seconds
    keep   = (B(:,1) >= tStart) & (B(:,1) <= tEnd);

    B = B(keep,:);
    L = L(keep,:);
    R = R(keep,:);

    % Use Back time as reference
    t = B(:,1);

    % Extract X-only signals
    aBx = B(:,2);  aLx = L(:,2);  aRx = R(:,2);  % accel X
    wBx = B(:,5);  wLx = L(:,5);  wRx = R(:,5);  % gyro X

    % Remove mean (centers signals to compare shape)
    if demeanSignals
        aBx = aBx - mean(aBx);  aLx = aLx - mean(aLx);  aRx = aRx - mean(aRx);
        wBx = wBx - mean(wBx);  wLx = wLx - mean(wLx);  wRx = wRx - mean(wRx);
    end

    % Duration label (actual plotted duration; may be <30 if trial is shorter)
    durPlot = t(end) - t(1);
    rowLabel = sprintf("%s | plotted %.1fs", trialField, durPlot);

    %% ---- Column 1: ACCEL X overlay (3 lines) ----
    nexttile;
    plot(t, aBx); hold on;
    plot(t, aLx);
    plot(t, aRx);
    grid on;
    ylabel(rowLabel, 'FontWeight','bold');
    title(sprintf("%s (Accel X)", upper(eventBase)));  % event name above the row (left tile)
    xlabel("Time (s)"); ylabel("Accel X (m/s^2)");
    legend("Back","Left Thigh","Right Thigh", "Location","best");

    %% ---- Column 2: GYRO X overlay (3 lines) ----
    nexttile;
    plot(t, wBx); hold on;
    plot(t, wLx);
    plot(t, wRx);
    grid on;
    title("Gyro X");
    xlabel("Time (s)"); ylabel("Gyro X (deg/s)");
    legend("Back","Left Thigh","Right Thigh", "Location","best");
end

%% ================== Helper Function to Ensure Finding Trial ==================
function trialField = pickTrialField_Tprefix(raw, pID, eventBase)
% Chooses trials that follow the naming convention T#_Event (e.g., T12_Sit).
% For each event, it prefers fields that END WITH "_Event". If multiple exist, picks the last.

    trialFields = string(fieldnames(raw.(pID)));

    % Enforce the suffix pattern "_Event"
    suffix = "_" + eventBase;

    hits = trialFields(endsWith(trialFields, suffix, 'IgnoreCase', true));

    if isempty(hits)
        % Fallback: contains match (in case naming isn't consistent)
        hits = trialFields(contains(trialFields, eventBase, 'IgnoreCase', true));
        if isempty(hits)
            trialField = "";
            return;
        end
    end

    % Pick last match if repeated (common pattern)
    trialField = hits(end);
end
