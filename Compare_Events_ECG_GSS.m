%% ===== P1960: ECG grid + GSS grid (separate figures) | first 30 seconds =====
clear; clc; close all;

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
        if any(strcmp(sf,"ECG")) && any(strcmp(sf,"GSS"))
            raw = S.(nm);
            fprintf("Using RAW sensor struct: S.%s\n", nm);
            break
        end
    end
end

if isempty(raw)
    error("Could not auto-find the raw sensor struct containing ECG/GSS.");
end

%% 3) Events (rows)
eventList = ["Walk", "Sit", "Slip"];

%% 4) Plot settings
ds = 10;      % downsample for speed (set ds=1 for full resolution)
tMax = 30;    % plot 30 seconds only

%% =========================
% Figure 1: ECG only
%% =========================
figure('Color','w');
tl1 = tiledlayout(numel(eventList), 1, 'TileSpacing','compact', 'Padding','compact');
title(tl1, sprintf("%s | ECG | First %.0f s", pID, tMax));

for r = 1:numel(eventList)
    eventBase = eventList(r);
    trialField = pickTrialField_Tprefix(raw, pID, eventBase);

    if strlength(trialField) == 0
        nexttile; axis off; text(0,0.5,"No matching trial for: " + eventBase, 'FontSize', 11);
        continue;
    end

    X = raw.(pID).(trialField);
    ECG = X.ECG;               % Nx2: [time(s), voltage]
    E = ECG(1:ds:end,:);

    % Crop first 30 s of that trial
    tStart = E(1,1) + 10;
    tEnd   = tStart + 30;
    keep   = (E(:,1) >= tStart) & (E(:,1) <= tEnd);
    E = E(keep,:);

    t = E(:,1);
    v = E(:,2);

    durPlot = t(end) - t(1);

    nexttile;
    plot(t, v); grid on;
    title(sprintf("%s (ECG) | %s | plotted %.1fs", upper(eventBase), trialField, durPlot));
    xlabel("Time (s)"); ylabel("ECG (V)");
end

%% =========================
% Figure 2: GSS only
%% =========================
figure('Color','w');
tl2 = tiledlayout(numel(eventList), 1, 'TileSpacing','compact', 'Padding','compact');
title(tl2, sprintf("%s | GSS | First %.0f s", pID, tMax));

for r = 1:numel(eventList)
    eventBase = eventList(r);
    trialField = pickTrialField_Tprefix(raw, pID, eventBase);

    if strlength(trialField) == 0
        nexttile; axis off; text(0,0.5,"No matching trial for: " + eventBase, 'FontSize', 11);
        continue;
    end

    X = raw.(pID).(trialField);
    GSS = X.GSS;               % Nx2: [time(s), voltage]
    G = GSS(1:ds:end,:);

    % Crop first 30 s of that trial
    tStart = G(1,1) + 10;
    tEnd   = tStart + tMax;
    keep   = (G(:,1) >= tStart) & (G(:,1) <= tEnd);
    G = G(keep,:);

    t = G(:,1);
    v = G(:,2);

    durPlot = t(end) - t(1);

    nexttile;
    plot(t, v); grid on;
    title(sprintf("%s (GSS) | %s | plotted %.1fs", upper(eventBase), trialField, durPlot));
    xlabel("Time (s)"); ylabel("GSS (V)");
end

%% ================== Helper Function to Ensure Finding Trial ==================
function trialField = pickTrialField_Tprefix(raw, pID, eventBase)
% Chooses trials that follow the naming convention T#_Event (e.g., T12_Sit).
% Prefers fields that END WITH "_Event". If multiple exist, picks the last.

    trialFields = string(fieldnames(raw.(pID)));
    suffix = "_" + eventBase;

    hits = trialFields(endsWith(trialFields, suffix, 'IgnoreCase', true));
    if isempty(hits)
        hits = trialFields(contains(trialFields, eventBase, 'IgnoreCase', true));
        if isempty(hits)
            trialField = "";
            return;
        end
    end

    trialField = hits(end);
end
