%% Information to remember
% IMU columns:
%   1 = time (s)
%   2–4 = accel X,Y,Z (m/s^2)
%   5–7 = gyro  X,Y,Z (deg/s)

%% ===== GCAssignment1: Gait Analysis of 3 participants (X-axis) =====
clear; clc; close all;

%% 1) Load MAT file
filePath = "C:\Users\callo\OneDrive\Desktop\University\2025-26\Wearbles\GrandChallengeData.mat";
S = load(filePath);

disp("Top-level variables in MAT file:");
disp(fieldnames(S))

%% 2) Participants to analyze
participants = ["P1960"];

%% 3) Auto-find raw sensor struct (use first participant as reference)
names = fieldnames(S);
raw = [];
refPID = participants(1);

for k = 1:numel(names)
    nm = names{k};
    if strcmp(nm, "clean_labels"), continue; end
    if ~isstruct(S.(nm)), continue; end
    if ~isfield(S.(nm), refPID), continue; end

    trialFields = fieldnames(S.(nm).(refPID));
    if isempty(trialFields), continue; end

    t1 = trialFields{1};
    if isstruct(S.(nm).(refPID).(t1))
        sf = fieldnames(S.(nm).(refPID).(t1));
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

%% 4) Trial settings
trialBaseName = "Walk";
t0 = 10;     % seconds
t1 = 70;     % seconds

% Choose which signal defines gait cycles:
cycleFrom = "accel";  % "accel" or "gyro"

%% 5) Loop through participants
for p = 1:numel(participants)
    pID = participants(p);
    fprintf("\n=== %s ===\n", pID);
    plotGaitX_withPeaks(raw, pID, trialBaseName, t0, t1, cycleFrom);
end

%% ================== Function ==================
function plotGaitX_withPeaks(raw, pID, trialBaseName, t0, t1, cycleFrom)

    % --- Find the trial field (handles names like T14_Walk, Walk_2, etc.) ---
    trialFields = string(fieldnames(raw.(pID)));

    if any(strcmpi(trialFields, trialBaseName))
        trialField = trialBaseName;
    else
        hits = trialFields(contains(trialFields, "_" + trialBaseName, 'IgnoreCase', true) | ...
                           contains(trialFields, trialBaseName, 'IgnoreCase', true));
        if isempty(hits)
            error("Trial '%s' not found for %s.", trialBaseName, pID);
        end
        trialField = hits(end);
    end

    X = raw.(pID).(trialField);

    % --- Pull thigh IMUs ---
    LTh = X.Left_Thigh;
    RTh = X.Right_Thigh;

    % --- Crop to time window ---
    crop = @(M) M(M(:,1)>=t0 & M(:,1)<=t1, :);
    Lc = crop(LTh);
    Rc = crop(RTh);

    % --- Extract X-direction signals ---
    tL = Lc(:,1);
    tR = Rc(:,1);

    aXL = Lc(:,2);   % accel X
    aXR = Rc(:,2);   % accel X
    wXL = Lc(:,5);   % gyro X
    wXR = Rc(:,5);   % gyro X

    % Center signals
    aXL = aXL - mean(aXL);  aXR = aXR - mean(aXR);
    wXL = wXL - mean(wXL);  wXR = wXR - mean(wXR);

    % Light smoothing to reduce tiny noise without removing gait peaks
    aXLs = smoothdata(aXL,'movmean',25);
    aXRs = smoothdata(aXR,'movmean',25);
    wXLs = smoothdata(wXL,'movmean',25);
    wXRs = smoothdata(wXR,'movmean',25);

    %% --- Peak detection (LEFT thigh defines cycle markers) ---
    % Actual duration in your chosen window:
    duration_s = tL(end) - tL(1);

    % Estimate sampling rate
    dt = median(diff(tL));
    fs = 1/dt;

    % Choose the signal for peak detection
    if cycleFrom == "gyro"
        sig = wXLs;
        sigName = "Left thigh gyro X";
    else
        sig = aXLs;
        sigName = "Left thigh accel X";
    end

    % Peak detection settings:
    % - MinPeakDistance prevents counting multiple peaks inside one step
    % - MinPeakProminence focuses on the "highest/most meaningful" peaks
    minPeakDist_samp = max(1, round(0.8 * fs));     % ~0.8 s separation
    minProm = 1.0 * std(sig);                       % research suggests 0.4-1.0

    % Use ~ to ignore peak heights when these can be ignored (cuts down
    % time)
    [pks, locs] = findpeaks(sig, ...
        'MinPeakDistance', minPeakDist_samp, ...
        'MinPeakProminence', minProm);

    peakTimes = tL(locs);
    nPeaks = numel(peakTimes);

    % Cadence: #peaks * 2 legs * (60 / duration)
    cadence = nPeaks * 2 * (60 / duration_s);

    % Print results
    fprintf("Trial: %s | Window: %.1f–%.1f s | Duration: %.2f s\n", trialField, t0, t1, duration_s);
    fprintf("Cycle markers from: %s\n", sigName);
    fprintf("Detected peaks (cycles) = %d\n", nPeaks);
    fprintf("Cadence = %d * (60 / %.2f) = %.2f peaks/min\n", nPeaks, duration_s, cadence);

    %% --- Plot: 4x1 grid + peak lines ---
    figure('Color','w');
    tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('%s | %s | Gait in X (%.1f–%.1fs) | Peaks=%d | Cadence=%.1f/min', ...
        pID, trialField, t0, t1, nPeaks, cadence));

    nexttile;
    plot(tL, aXLs); grid on; hold on;
    addPeakLines(peakTimes);
    plot(peakTimes, aXLs(locs), 'o'); % show peak points on left accel
    title('Left Thigh Accel X'); xlabel('Time (s)'); ylabel('Accel X (m/s^2)');

    nexttile;
    plot(tR, aXRs); grid on; hold on;
    addPeakLines(peakTimes);
    title('Right Thigh Accel X'); xlabel('Time (s)'); ylabel('Accel X (m/s^2)');

    nexttile;
    plot(tL, wXLs); grid on; hold on;
    addPeakLines(peakTimes);
    title('Left Thigh Gyro X'); xlabel('Time (s)'); ylabel('Gyro X (deg/s)');

    nexttile;
    plot(tR, wXRs); grid on; hold on;
    addPeakLines(peakTimes);
    title('Right Thigh Gyro X'); xlabel('Time (s)'); ylabel('Gyro X (deg/s)');
end

function addPeakLines(peakTimes)
    for i = 1:numel(peakTimes)
        xline(peakTimes(i), '--');
    end
end
