%% ===== GCAssignment1: Visualize two participants (P1960 and P2070) =====
clear; clc;

%% 1) Load MAT file
filePath = "/Users/sophiakatramadakis/Documents/MATLAB/GrandChallengeData.mat";
S = load(filePath);

disp("Top-level variables in MAT file:");
disp(fieldnames(S))

%% 2) Assign labels + auto-find raw sensor struct
if ~isfield(S, "clean_labels")
    error("Expected S.clean_labels based on your screenshot, but it wasn't found. Check fieldnames(S).");
end
labels = S.clean_labels;

names = fieldnames(S);
raw = [];

for k = 1:numel(names)
    nm = names{k};
    if strcmp(nm, "clean_labels"), continue; end
    if ~isstruct(S.(nm)), continue; end
    if ~isfield(S.(nm), "P1960"), continue; end

    trialFields = fieldnames(S.(nm).P1960);
    if isempty(trialFields), continue; end

    t1 = trialFields{1};
    if isstruct(S.(nm).P1960.(t1))
        sf = fieldnames(S.(nm).P1960.(t1));
        if any(strcmp(sf,"Back")) && any(strcmp(sf,"ECG")) && any(strcmp(sf,"GSS"))
            raw = S.(nm);
            fprintf("Using RAW sensor struct: S.%s\n", nm);
            break
        end
    end
end

if isempty(raw)
    error("Could not auto-find the raw sensor struct. Paste fieldnames(S) and I’ll point to it.");
end

%% 3) Participants to view
participantsToPlot = ["P1960","P2070","P2995"];

%% 4) Trials to use (same logic for both participants)
perturbTrial = "Slip";     % we will take the first near-fall + first fall block inside Slip
adlTrial     = "Walk";     % plot first 20 seconds

%% 5) Loop participants
for pID = participantsToPlot
    fprintf("\n=== Plotting participant %s ===\n", pID);

    % --- Build near-fall and fall windows from labels (in Slip) ---
    if ~isfield(labels, pID) || ~isfield(labels.(pID), perturbTrial)
        warning("No labels found for %s %s. Skipping.", pID, perturbTrial);
        continue;
    end

    T = labels.(pID).(perturbTrial);  % Mx3: [label, time_ms, group]

    % Near-fall: first group==1 block (expect 3 rows: 1,2,3)
    nf_rows = find(T(:,3)==1);
    if numel(nf_rows) >= 3
        nf_block = T(nf_rows(1:3),:);
        t_nf0 = nf_block(1,2)/1000;
        t_nf1 = nf_block(3,2)/1000 + 2;
        plotSensorsWindow(raw, labels, pID, perturbTrial, t_nf0, t_nf1, true, "NEAR-FALL window (Slip)");
    else
        warning("No near-fall block (group 1) found in %s %s.", pID, perturbTrial);
    end

    % Fall: first group==4 block (expect 4 rows: 1,2,4,5)
    fall_rows = find(T(:,3)==4);
    if numel(fall_rows) >= 4
        fall_block = T(fall_rows(1:4),:);
        t_f0 = fall_block(1,2)/1000;
        t_f1 = fall_block(4,2)/1000 + 2;
        plotSensorsWindow(raw, labels, pID, perturbTrial, t_f0, t_f1, true, "FALL window (Slip)");
    else
        warning("No fall block (group 4) found in %s %s.", pID, perturbTrial);
    end

    % ADL: first 20 seconds of Walk (no labels)
    plotSensorsWindow(raw, labels, pID, adlTrial, 0, 20, false, "ADL window (Walk) first 20s");
end

%% ================== Functions ==================
function plotSensorsWindow(raw, labels, pID, trialBaseName, t0, t1, showMarkers, subtitleTxt)

    % Handle repeated trials like Slip_2 / Walk_2 by picking the last match
    trialFields = string(fieldnames(raw.(pID)));
    if any(strcmp(trialFields, trialBaseName))
        trialField = trialBaseName;
    else
        hits = trialFields(contains(trialFields, trialBaseName, 'IgnoreCase', true));
        if isempty(hits)
            error("Trial '%s' not found in sensor data for %s.", trialBaseName, pID);
        end
        trialField = hits(end);
    end

    X = raw.(pID).(trialField);

    Back = X.Back;  LTh = X.Left_Thigh;  RTh = X.Right_Thigh;
    ECG  = X.ECG;   GSS = X.GSS;

    % Timing markers
    timing_s = []; timingLabel = [];
    if showMarkers
        base = regexprep(trialBaseName, "_\d+$", "");
        if isfield(labels.(pID), base)
            T = labels.(pID).(base);
            timingLabel = T(:,1);
            timing_s = T(:,2)/1000;
        end
    end

    crop = @(M) M(M(:,1)>=t0 & M(:,1)<=t1, :);

    Bc = crop(Back); Lc = crop(LTh); Rc = crop(RTh);
    Ec = crop(ECG);  Gc = crop(GSS);

    aB = vecnorm(Bc(:,2:4),2,2); wB = vecnorm(Bc(:,5:7),2,2);
    aL = vecnorm(Lc(:,2:4),2,2); wL = vecnorm(Lc(:,5:7),2,2);
    aR = vecnorm(Rc(:,2:4),2,2); wR = vecnorm(Rc(:,5:7),2,2);

    figure('Color','w');
    tl = tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
    title(tl, pID + " | " + trialField + " | " + subtitleTxt + sprintf(" (%.1f–%.1fs)", t0, t1));

    nexttile; plot(Bc(:,1), aB); hold on; plot(Bc(:,1), wB);
    ylabel("Back |a|,|w|"); legend("|a|","|w|","Location","best");
    addLines(timing_s,timingLabel,t0,t1);

    nexttile; plot(Lc(:,1), aL); hold on; plot(Lc(:,1), wL);
    ylabel("L thigh |a|,|w|"); legend("|a|","|w|","Location","best");
    addLines(timing_s,timingLabel,t0,t1);

    nexttile; plot(Rc(:,1), aR); hold on; plot(Rc(:,1), wR);
    ylabel("R thigh |a|,|w|"); legend("|a|","|w|","Location","best");
    addLines(timing_s,timingLabel,t0,t1);

    nexttile; plot(Ec(:,1), Ec(:,2));
    ylabel("ECG (V)");
    addLines(timing_s,timingLabel,t0,t1);

    nexttile; plot(Gc(:,1), Gc(:,2));
    ylabel("GSS (V)"); xlabel("Time (s)");
    addLines(timing_s,timingLabel,t0,t1);
end

function addLines(timing_s, timingLabel, t0, t1)
    if isempty(timing_s), return; end
    inWin = timing_s>=t0 & timing_s<=t1;
    ts = timing_s(inWin);
    ls = timingLabel(inWin);

    for i = 1:numel(ts)
        xline(ts(i),'--',"L"+ls(i),'LabelOrientation','horizontal');
    end
end
