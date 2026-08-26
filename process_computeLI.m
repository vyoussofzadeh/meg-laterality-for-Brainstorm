function varargout = process_computeLI(varargin )

% PROCESS_COMPUTELI_HCP: Compute the Lateralization Index (LI) using various methods (source magnitude, counting, bootstrapping)
% on MEG source-level data projected onto the HCP atlas.
%
% This process integrates into Brainstorm's processing pipeline and allows selecting:
%   - Time interval method (specific, averaged, window-based)
%   - LI computation method (source magnitude, counting, bootstrapping)
%   - Thresholding approach (global-max, time-max, region-max)
%   - Additional output options such as saving .mat files or generating plots for the window-based approach.
%
% REQUIREMENTS:
%   - Brainstorm environment
%   - MEG source results
%   - HCP MMP1.0 atlas integrated into the subject anatomy
%
% Author: Vahab Youssof Zadeh, 2025
% Last Update: 08/26/26
% Changes: Fixed/optimized window analysis and unified MAT/TSV export across interval modes.

eval(macro_method);

end

function sProcess = GetDescription() 
% Configure the Brainstorm process menu for lateralization-index analysis.

sProcess.Comment     = 'Lateralization index (LI): HCP / DK';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 337;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/Scouts';
sProcess.InputTypes  = {'results'};
sProcess.OutputTypes = {'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

%% 1. Atlas
sProcess.options.atlasSection.Comment = '<B>1. Atlas</B>';
sProcess.options.atlasSection.Type    = 'label';
sProcess.options.atlasSection.Value   = {};

sProcess.options.atlas.Comment = 'Cortical atlas:';
sProcess.options.atlas.Type    = 'combobox';
sProcess.options.atlas.Value   = {1, {'HCP-MMP1.0 (symmetric)', 'Desikan-Killiany (DK)'}};

sProcess.options.atlas_note.Comment = ['<I>HCP:</I> requires the symmetric MMP atlas ' ...
    '(<I>mmp_in_mni_symmetrical_1.nii</I>). <I>DK:</I> uses the subject''s built-in Desikan-Killiany atlas.'];
sProcess.options.atlas_note.Type  = 'label';
sProcess.options.atlas_note.Value = {};

%% 2. Time analysis
sProcess.options.timeSection.Comment = '<B>2. Time analysis</B>';
sProcess.options.timeSection.Type    = 'label';
sProcess.options.timeSection.Value   = {};

sProcess.options.window.Comment = 'Analysis mode:';
sProcess.options.window.Type    = 'combobox';
sProcess.options.window.Value   = {1, {'Single time interval', 'Average across interval', 'Sliding windows'}};

sProcess.options.poststim_custom.Comment = 'Analysis interval:';
sProcess.options.poststim_custom.Type    = 'poststim';
sProcess.options.poststim_custom.Value   = [];

sProcess.options.windowSettings.Comment = '<I>Sliding-window settings</I>';
sProcess.options.windowSettings.Type    = 'label';
sProcess.options.windowSettings.Value   = {};

sProcess.options.twindow1.Comment = 'Window length:';
sProcess.options.twindow1.Type    = 'value';
sProcess.options.twindow1.Value   = {300, 'ms', 0};

sProcess.options.toverlap.Comment = 'Window overlap:';
sProcess.options.toverlap.Type    = 'value';
sProcess.options.toverlap.Value   = {50, '%', 0, 99, 1};

sProcess.options.window_note.Comment = ['<I>Example:</I> A 300 ms window with 50% overlap advances by 150 ms. ' ...
    'Window settings are used only for Sliding windows.'];
sProcess.options.window_note.Type  = 'label';
sProcess.options.window_note.Value = {};

%% 3. LI method
sProcess.options.methodSection.Comment = '<B>3. LI calculation method</B>';
sProcess.options.methodSection.Type    = 'label';
sProcess.options.methodSection.Value   = {};

sProcess.options.methodSource.Comment = 'Source magnitude (suprathreshold amplitude)';
sProcess.options.methodSource.Type    = 'checkbox';
sProcess.options.methodSource.Value   = 0;

sProcess.options.methodCounting.Comment = 'Vertex count (suprathreshold observations)';
sProcess.options.methodCounting.Type    = 'checkbox';
sProcess.options.methodCounting.Value   = 0;

sProcess.options.methodBootstrap.Comment = 'Bootstrap LI with 95% confidence interval';
sProcess.options.methodBootstrap.Type    = 'checkbox';
sProcess.options.methodBootstrap.Value   = 0;

sProcess.options.method_note.Comment = '<I>You may select one or multiple LI methods.</I>';
sProcess.options.method_note.Type  = 'label';
sProcess.options.method_note.Value = {};

%% 4. Source values
sProcess.options.effectSection.Comment = '<B>4. Source values</B>';
sProcess.options.effectSection.Type    = 'label';
sProcess.options.effectSection.Value   = {};

sProcess.options.effect.Comment = 'Values to analyze:';
sProcess.options.effect.Type    = 'combobox';
sProcess.options.effect.Value   = {1, {'Positive', 'Negative (sign reversed)', 'Absolute magnitude'}};

%% 5. Thresholding for source/count methods
sProcess.options.thresholdSection.Comment = '<B>5. Thresholding</B> <I>(source magnitude and vertex count)</I>';
sProcess.options.thresholdSection.Type    = 'label';
sProcess.options.thresholdSection.Value   = {};

sProcess.options.threshtype1.Comment = 'Reference maximum:';
sProcess.options.threshtype1.Type    = 'combobox';
sProcess.options.threshtype1.Value   = {1, {'Whole-recording maximum', 'Selected-interval maximum', 'ROI maximum'}};

sProcess.options.ratio4threshold.Comment = 'Threshold level:';
sProcess.options.ratio4threshold.Type    = 'value';
sProcess.options.ratio4threshold.Value   = {50, '%', 0, 100, 1, 1};

sProcess.options.threshold_note.Comment = ['<I>Sliding windows:</I> source/count methods use a consistent per-ROI maximum across windows. ' ...
    'Bootstrap LI instead uses the threshold divisions below.'];
sProcess.options.threshold_note.Type  = 'label';
sProcess.options.threshold_note.Value = {};

%% 6. Bootstrap
sProcess.options.bootstrapSection.Comment = '<B>6. Bootstrap settings</B>';
sProcess.options.bootstrapSection.Type    = 'label';
sProcess.options.bootstrapSection.Value   = {};

sProcess.options.divs.Comment = 'Threshold divisions:';
sProcess.options.divs.Type    = 'value';
sProcess.options.divs.Value   = {10, '', 1, 100, 1};

sProcess.options.n_resampling.Comment = 'Bootstrap samples per threshold:';
sProcess.options.n_resampling.Type    = 'value';
sProcess.options.n_resampling.Value   = {20, '', 1, 1000, 1};

sProcess.options.RESAMPLE_RATIO.Comment = 'Sample fraction:';
sProcess.options.RESAMPLE_RATIO.Type    = 'value';
sProcess.options.RESAMPLE_RATIO.Value   = {75, '%', 1, 100, 1, 1};

sProcess.options.bootstrap_note.Comment = ['<I>Used only for Bootstrap LI.</I> Runtime increases with the number of windows, ' ...
    'threshold divisions, and bootstrap samples.'];
sProcess.options.bootstrap_note.Type  = 'label';
sProcess.options.bootstrap_note.Value = {};

%% 7. Output
sProcess.options.outputSection.Comment = '<B>7. Output</B>';
sProcess.options.outputSection.Type    = 'label';
sProcess.options.outputSection.Value   = {};

sProcess.options.savedir.Comment = 'Results folder:';
sProcess.options.savedir.Type    = 'text';
sProcess.options.savedir.Value   = '';

sProcess.options.sname.Comment = 'Output file name:';
sProcess.options.sname.Type    = 'text';
sProcess.options.sname.Value   = 'LI_results';

% Keep the legacy field name "saveMat" so existing Brainstorm pipelines
% continue to load, but use it as the general save-results switch.
sProcess.options.saveMat.Comment = 'Save LI results';
sProcess.options.saveMat.Type    = 'checkbox';
sProcess.options.saveMat.Value   = 0;

sProcess.options.saveFormat.Comment = 'File format:';
sProcess.options.saveFormat.Type    = 'combobox';
sProcess.options.saveFormat.Value   = {1, {'MAT file (.mat)', 'Tab-delimited text (.tsv)'}};

sProcess.options.plotResults.Comment = 'Plot LI across sliding windows';
sProcess.options.plotResults.Type    = 'checkbox';
sProcess.options.plotResults.Value   = 1;

sProcess.options.output_note.Comment = ['<I>Saving applies to Single interval, Average across interval, and Sliding windows. ' ...
    'Temporal plotting applies only to Sliding windows.</I>'];
sProcess.options.output_note.Type  = 'label';
sProcess.options.output_note.Value = {};

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
Comment = sProcess.Comment;
end

%%
function OutputFiles = Run(sProcess, sInput)

% Initialize output
OutputFiles = {};

% Load the results file
sResultP = in_bst_results(sInput.FileName, 1);
if isempty(sResultP)
    error('Could not load the selected results file.');
end

% Stop any progress bars
bst_progress('stop');

% Extract user-defined options
plotResults = sProcess.options.plotResults.Value;
saveResults = sProcess.options.saveMat.Value; % Legacy option field, now applies to all modes
if isfield(sProcess.options, 'saveFormat')
    saveFormat = sProcess.options.saveFormat.Value{1};
else
    saveFormat = 1; % Backward-compatible default: MAT
end
savedir = sProcess.options.savedir.Value;
sname   = sProcess.options.sname.Value;

% Extract atlas selection
atlasChoice = sProcess.options.atlas.Value{1};
% atlasChoice == 1: HCP;
% atlasChoice == 2: DK

% Extract time interval method
Tinterval = selectTimeInterval(sProcess.options.window.Value{1});

% Determine the time range based on user selection
if Tinterval == 1 || Tinterval == 3
    timerange = sProcess.options.poststim_custom.Value{1};
else
    timerange = 1; % For averaged or other intervals, default/legacy behavior
end

% Extract effect type, threshold parameters
effect          = selectEffectType(sProcess.options.effect.Value{1});
Ratio4Threshold = sProcess.options.ratio4threshold.Value{1}/100;
Threshtype      = selectThresholdType(sProcess.options.threshtype1.Value{1});

% Process ImageGridAmp based on selected effect
ImageGridAmp = processImageGridAmp(sResultP.ImageGridAmp, effect);

% Determine max values for various time windows
[AllMax, GlobalMax, t1, t2] = determineMaxValues(Tinterval, ImageGridAmp, sResultP, timerange);

% If user chose Averaged Time Interval and sufficient data points:
if sProcess.options.window.Value{1} == 2 && length(sResultP.Time) > 10
    timerange = sProcess.options.poststim_custom.Value{1};
    [~, t1] = min(abs(sResultP.Time - timerange(1)));
    [~, t2] = min(abs(sResultP.Time - timerange(2)));
    ImageGridAmp = mean(ImageGridAmp(:, t1:t2), 2);
end

% If window-based, define time windows
if sProcess.options.window.Value{1} == 3
    % Brainstorm converts the menu value with unit "ms" to seconds before Run.
    winLengthSec = sProcess.options.twindow1.Value{1};
    overlapPct = sProcess.options.toverlap.Value{1};
    
    % Reject windows shorter than one data sample: they repeatedly select the
    % same time point and usually indicate a unit-conversion mistake.
    sampleSteps = abs(diff(sResultP.Time(:)));
    sampleSteps = sampleSteps(isfinite(sampleSteps) & sampleSteps > 0);
    if ~isempty(sampleSteps) && winLengthSec < median(sampleSteps)
        error(['Window length (%.3f ms) is shorter than one data sample (%.3f ms). ' ...
            'Increase the window length.'], 1000*winLengthSec, 1000*median(sampleSteps));
    end
    
    cfg = [];
    cfg.strt = timerange(1);
    cfg.spt = timerange(2);
    cfg.overlap = overlapPct/100;
    cfg.linterval = winLengthSec;
    wi = generateTimeWindows(cfg);
    
    nWindows = size(wi, 1);
    if nWindows == 0
        error('The selected interval is shorter than the requested window length.');
    end
    if sProcess.options.methodBootstrap.Value == 1 && nWindows > 1000
        error(['The current settings generate %d bootstrap windows, which is likely unintended. ' ...
            'Increase the window length or reduce the overlap.'], nWindows);
    end
else
    wi = [];
end

if sProcess.options.window.Value{1} == 3
    fprintf('\nSelected sliding-window analysis:\n');
    fprintf('  Interval: [%.4f, %.4f] s\n', timerange(1), timerange(2));
    fprintf('  Window length: %.1f ms | Overlap: %.1f%% | Windows: %d\n', ...
        1000*winLengthSec, overlapPct, nWindows);
    fprintf('  First window: [%.4f, %.4f] s', wi(1,1), wi(1,2));
    if nWindows > 1
        fprintf(' | Last window: [%.4f, %.4f] s', wi(end,1), wi(end,2));
    end
    fprintf('\n');
elseif sProcess.options.window.Value{1} == 1
    fprintf('\nSelected interval: [%.4f, %.4f] s\n', timerange(1), timerange(2));
end

% Convert atlas and define ROIs based on user choice
switch atlasChoice
    case 1
        % HCP atlas selected
        [sScout, ~] = convertHCPScout(sResultP);
        [RoiLabels, RoiIndices] = defineROIs_HCP(sScout);
    case 2
        % Desikan-Killiany atlas selected
        [sScout, ~] = convertDesikanKillianyScout(sResultP);
        [RoiLabels, RoiIndices] = defineROIs_DK(); % DK ROIs defined in defineROIs()
    otherwise
        error('Invalid atlas selection.');
end

% Prepare configuration for LI computation
cfg_LI = [];
cfg_LI.Tinterval       = Tinterval;
cfg_LI.ImageGridAmp    = ImageGridAmp;
cfg_LI.timerange       = timerange;
cfg_LI.RoiLabels       = RoiLabels;
cfg_LI.RoiIndices      = RoiIndices;
cfg_LI.sScout          = sScout;
cfg_LI.AllMax          = AllMax;
cfg_LI.GlobalMax       = GlobalMax;
cfg_LI.Threshtype      = Threshtype;
cfg_LI.Ratio4Threshold = Ratio4Threshold;
cfg_LI.t1              = t1;
cfg_LI.t2              = t2;
cfg_LI.plotResults     = plotResults;
cfg_LI.saveResults     = saveResults;
cfg_LI.saveFormat      = saveFormat;
cfg_LI.saveMat         = saveResults; % Legacy alias
cfg_LI.savedir         = savedir;
cfg_LI.sname           = sname;
cfg_LI.Time            = sResultP.Time;
cfg_LI.Comment         = sResultP.Comment;
cfg_LI.windows         = wi; % Only meaningful if window-based

% Determine which methods are selected
methodsSelected = [sProcess.options.methodSource.Value, sProcess.options.methodCounting.Value, sProcess.options.methodBootstrap.Value];
if ~any(methodsSelected)
    warning('No LI computation method selected. Please choose at least one method.');
    return;
end

% Source Magnitude Method
if sProcess.options.methodSource.Value == 1
    disp('Computing LI using Source Magnitude Method...')
    cfg_LI.method = 1;
    computeLI(cfg_LI);
    pause(0.2);
end

% Counting Method
if sProcess.options.methodCounting.Value == 1
    disp('Computing LI using Counting Method...')
    cfg_LI.method = 2;
    computeLI(cfg_LI);
    pause(0.2);
end

% Bootstrapping Method
if sProcess.options.methodBootstrap.Value == 1
    disp('Computing LI with Bootstrapping...')
    cfg_LI.method = 3;
    cfg_LI.divs           = sProcess.options.divs.Value{1};
    cfg_LI.n_resampling   = sProcess.options.n_resampling.Value{1};
    cfg_LI.RESAMPLE_RATIO = sProcess.options.RESAMPLE_RATIO.Value{1} / 100;
    computeLI(cfg_LI);
    pause(0.2);
end

if Tinterval == 1 || Tinterval == 3
    disp(['LI assessed for interval: ', num2str(timerange(1)), '-', num2str(timerange(2)), ' sec using the selected atlas.']);
end
disp('If needed, edit process_computeLI.m after ensuring Brainstorm is running.')
disp('LI analysis is completed!')

% Set OutputFiles if you create any output files or results
OutputFiles = {}; % Adjust this if you generate and want to return any results

end

% === HELPER FUNCTIONS ===
function Tinterval = selectTimeInterval(Tinterval)
% Prompt user to select time interval

% Ensure valid selection
if isempty(Tinterval) || ~any(Tinterval == [1, 2, 3])
    error('Invalid time interval selection. Choose 1, 2, or 3.');
end
end

function effect = selectEffectType(effect)
% Prompt user to select effect type

% Ensure valid selection
if isempty(effect) || ~any(effect == [1, 2, 3])
    error('Invalid effect type selection. Choose 1, 2, or 3.');
end
end

function Threshtype = selectThresholdType(Threshtype)
% Prompt user to select threshold type

% Ensure valid selection
if isempty(Threshtype) || ~any(Threshtype == [1, 2, 3])
    error('Invalid threshold type selection. Choose 1, 2, or 3.');
end
end

function ImageGridAmp = processImageGridAmp(ImageGridAmp, effect)
% Apply the desired effect on the ImageGridAmp

switch effect
    case 1
        % Positive values: No change needed, as ImageGridAmp remains the same
    case 2
        ImageGridAmp = -ImageGridAmp;
    case 3
        ImageGridAmp = abs(ImageGridAmp);
    otherwise
        error('Invalid effect type. Choose 1 (Positive), 2 (Negative), or 3 (Absolute).');
end
end

function [AllMax, GlobalMax, t1, t2] = determineMaxValues(Tinterval, ImageGridAmp, sResultP, timerange)
% Compute the maximum values AllMax and GlobalMax

switch Tinterval
    case 2
        GlobalMax = max(ImageGridAmp(:));  % Max value over all time points
        AllMax = max(ImageGridAmp(:));     % Max value over the time window of interest
        t1 = []; t2 = [];
    case {1,3}
        t1 = find(sResultP.Time >= timerange(1), 1);
        t2 = find(sResultP.Time >= timerange(2), 1);
        AllMax = max(max(ImageGridAmp(:, t1:t2)));   % Max value over the time window of interest
        GlobalMax = max(max(ImageGridAmp));           % Max value over all time points
end
end

function [sScout, ProtocolInfo] = convertHCPScout(sResultP)

ProtocolInfo = bst_get('ProtocolInfo');
SurfaceFile = load(fullfile(ProtocolInfo.SUBJECTS, sResultP.SurfaceFile));

Scouts = [];
sScout = [];
for i = 1:length(SurfaceFile.Atlas)
    if contains(SurfaceFile.Atlas(i).Name, {'mmp_in_mni_symmetrical_1.nii'})
        Scouts = SurfaceFile.Atlas(i).Scouts;
    end
end
sScout.Scouts = Scouts;

expectedRegions = {'L_V1_ROI'	'R_V1_ROI'	'L_FEF_ROI'	'R_FEF_ROI'	'L_OP4_ROI'	'R_OP4_ROI'	'L_OP1_ROI'	'R_OP1_ROI'	'L_OP2-3_ROI'	'R_OP2-3_ROI'	'L_52_ROI'	'R_52_ROI'	'L_RI_ROI'	'R_RI_ROI'	'L_PFcm_ROI'	'R_PFcm_ROI'	'L_PoI2_ROI'	'R_PoI2_ROI'	'L_TA2_ROI'	'R_TA2_ROI'	'L_FOP4_ROI'	'R_FOP4_ROI'	'L_MI_ROI'	'R_MI_ROI'	'L_PEF_ROI'	'R_PEF_ROI'	'L_Pir_ROI'	'R_Pir_ROI'	'L_AVI_ROI'	'R_AVI_ROI'	'L_AAIC_ROI'	'R_AAIC_ROI'	'L_FOP1_ROI'	'R_FOP1_ROI'	'L_FOP3_ROI'	'R_FOP3_ROI'	'L_FOP2_ROI'	'R_FOP2_ROI'	'L_PFt_ROI'	'R_PFt_ROI'	'L_AIP_ROI'	'R_AIP_ROI'	'L_EC_ROI'	'R_EC_ROI'	'L_PreS_ROI'	'R_PreS_ROI'	'L_55b_ROI'	'R_55b_ROI'	'L_H_ROI'	'R_H_ROI'	'L_ProS_ROI'	'R_ProS_ROI'	'L_PeEc_ROI'	'R_PeEc_ROI'	'L_STGa_ROI'	'R_STGa_ROI'	'L_PBelt_ROI'	'R_PBelt_ROI'	'L_A5_ROI'	'R_A5_ROI'	'L_PHA1_ROI'	'R_PHA1_ROI'	'L_PHA3_ROI'	'R_PHA3_ROI'	'L_STSda_ROI'	'R_STSda_ROI'	'L_STSdp_ROI'	'R_STSdp_ROI'	'L_V3A_ROI'	'R_V3A_ROI'	'L_STSvp_ROI'	'R_STSvp_ROI'	'L_TGd_ROI'	'R_TGd_ROI'	'L_TE1a_ROI'	'R_TE1a_ROI'	'L_TE1p_ROI'	'R_TE1p_ROI'	'L_TE2a_ROI'	'R_TE2a_ROI'	'L_TF_ROI'	'R_TF_ROI'	'L_TE2p_ROI'	'R_TE2p_ROI'	'L_PHT_ROI'	'R_PHT_ROI'	'L_PH_ROI'	'R_PH_ROI'	'L_TPOJ1_ROI'	'R_TPOJ1_ROI'	'L_RSC_ROI'	'R_RSC_ROI'	'L_TPOJ2_ROI'	'R_TPOJ2_ROI'	'L_TPOJ3_ROI'	'R_TPOJ3_ROI'	'L_DVT_ROI'	'R_DVT_ROI'	'L_PGp_ROI'	'R_PGp_ROI'	'L_IP2_ROI'	'R_IP2_ROI'	'L_IP1_ROI'	'R_IP1_ROI'	'L_IP0_ROI'	'R_IP0_ROI'	'L_PFop_ROI'	'R_PFop_ROI'	'L_PF_ROI'	'R_PF_ROI'	'L_PFm_ROI'	'R_PFm_ROI'	'L_POS2_ROI'	'R_POS2_ROI'	'L_PGi_ROI'	'R_PGi_ROI'	'L_PGs_ROI'	'R_PGs_ROI'	'L_V6A_ROI'	'R_V6A_ROI'	'L_VMV1_ROI'	'R_VMV1_ROI'	'L_VMV3_ROI'	'R_VMV3_ROI'	'L_PHA2_ROI'	'R_PHA2_ROI'	'L_V4t_ROI'	'R_V4t_ROI'	'L_FST_ROI'	'R_FST_ROI'	'L_V3CD_ROI'	'R_V3CD_ROI'	'L_LO3_ROI'	'R_LO3_ROI'	'L_V7_ROI'	'R_V7_ROI'	'L_VMV2_ROI'	'R_VMV2_ROI'	'L_31pd_ROI'	'R_31pd_ROI'	'L_31a_ROI'	'R_31a_ROI'	'L_VVC_ROI'	'R_VVC_ROI'	'L_25_ROI'	'R_25_ROI'	'L_s32_ROI'	'R_s32_ROI'	'L_pOFC_ROI'	'R_pOFC_ROI'	'L_PoI1_ROI'	'R_PoI1_ROI'	'L_Ig_ROI'	'R_Ig_ROI'	'L_FOP5_ROI'	'R_FOP5_ROI'	'L_IPS1_ROI'	'R_IPS1_ROI'	'L_p10p_ROI'	'R_p10p_ROI'	'L_p47r_ROI'	'R_p47r_ROI'	'L_TGv_ROI'	'R_TGv_ROI'	'L_MBelt_ROI'	'R_MBelt_ROI'	'L_LBelt_ROI'	'R_LBelt_ROI'	'L_A4_ROI'	'R_A4_ROI'	'L_STSva_ROI'	'R_STSva_ROI'	'L_TE1m_ROI'	'R_TE1m_ROI'	'L_PI_ROI'	'R_PI_ROI'	'L_a32pr_ROI'	'R_a32pr_ROI'	'L_FFC_ROI'	'R_FFC_ROI'	'L_p24_ROI'	'R_p24_ROI'	'L_V3B_ROI'	'R_V3B_ROI'	'L_MST_ROI'	'R_MST_ROI'	'L_LO1_ROI'	'R_LO1_ROI'	'L_LO2_ROI'	'R_LO2_ROI'	'L_PIT_ROI'	'R_PIT_ROI'	'L_MT_ROI'	'R_MT_ROI'	'L_A1_ROI'	'R_A1_ROI'	'L_PSL_ROI'	'R_PSL_ROI'	'L_SFL_ROI'	'R_SFL_ROI'	'L_PCV_ROI'	'R_PCV_ROI'	'L_STV_ROI'	'R_STV_ROI'	'L_7Pm_ROI'	'R_7Pm_ROI'	'L_V6_ROI'	'R_V6_ROI'	'L_7m_ROI'	'R_7m_ROI'	'L_POS1_ROI'	'R_POS1_ROI'	'L_23d_ROI'	'R_23d_ROI'	'L_v23ab_ROI'	'R_v23ab_ROI'	'L_d23ab_ROI'	'R_d23ab_ROI'	'L_31pv_ROI'	'R_31pv_ROI'	'L_5m_ROI'	'R_5m_ROI'	'L_5mv_ROI'	'R_5mv_ROI'	'L_23c_ROI'	'R_23c_ROI'	'L_5L_ROI'	'R_5L_ROI'	'L_V2_ROI'	'R_V2_ROI'	'L_24dd_ROI'	'R_24dd_ROI'	'L_24dv_ROI'	'R_24dv_ROI'	'L_7AL_ROI'	'R_7AL_ROI'	'L_SCEF_ROI'	'R_SCEF_ROI'	'L_6ma_ROI'	'R_6ma_ROI'	'L_7Am_ROI'	'R_7Am_ROI'	'L_7PL_ROI'	'R_7PL_ROI'	'L_7PC_ROI'	'R_7PC_ROI'	'L_LIPv_ROI'	'R_LIPv_ROI'	'L_VIP_ROI'	'R_VIP_ROI'	'L_V3_ROI'	'R_V3_ROI'	'L_MIP_ROI'	'R_MIP_ROI'	'L_1_ROI'	'R_1_ROI'	'L_2_ROI'	'R_2_ROI'	'L_3a_ROI'	'R_3a_ROI'	'L_6d_ROI'	'R_6d_ROI'	'L_6mp_ROI'	'R_6mp_ROI'	'L_6v_ROI'	'R_6v_ROI'	'L_p24pr_ROI'	'R_p24pr_ROI'	'L_33pr_ROI'	'R_33pr_ROI'	'L_a24pr_ROI'	'R_a24pr_ROI'	'L_V4_ROI'	'R_V4_ROI'	'L_p32pr_ROI'	'R_p32pr_ROI'	'L_a24_ROI'	'R_a24_ROI'	'L_d32_ROI'	'R_d32_ROI'	'L_8BM_ROI'	'R_8BM_ROI'	'L_p32_ROI'	'R_p32_ROI'	'L_10r_ROI'	'R_10r_ROI'	'L_47m_ROI'	'R_47m_ROI'	'L_8Av_ROI'	'R_8Av_ROI'	'L_8Ad_ROI'	'R_8Ad_ROI'	'L_9m_ROI'	'R_9m_ROI'	'L_V8_ROI'	'R_V8_ROI'	'L_8BL_ROI'	'R_8BL_ROI'	'L_9p_ROI'	'R_9p_ROI'	'L_10d_ROI'	'R_10d_ROI'	'L_8C_ROI'	'R_8C_ROI'	'L_44_ROI'	'R_44_ROI'	'L_45_ROI'	'R_45_ROI'	'L_47l_ROI'	'R_47l_ROI'	'L_a47r_ROI'	'R_a47r_ROI'	'L_6r_ROI'	'R_6r_ROI'	'L_IFJa_ROI'	'R_IFJa_ROI'	'L_4_ROI'	'R_4_ROI'	'L_IFJp_ROI'	'R_IFJp_ROI'	'L_IFSp_ROI'	'R_IFSp_ROI'	'L_IFSa_ROI'	'R_IFSa_ROI'	'L_p9-46v_ROI'	'R_p9-46v_ROI'	'L_46_ROI'	'R_46_ROI'	'L_a9-46v_ROI'	'R_a9-46v_ROI'	'L_9-46d_ROI'	'R_9-46d_ROI'	'L_9a_ROI'	'R_9a_ROI'	'L_10v_ROI'	'R_10v_ROI'	'L_a10p_ROI'	'R_a10p_ROI'	'L_3b_ROI'	'R_3b_ROI'	'L_10pp_ROI'	'R_10pp_ROI'	'L_11l_ROI'	'R_11l_ROI'	'L_13l_ROI'	'R_13l_ROI'	'L_OFC_ROI'	'R_OFC_ROI'	'L_47s_ROI'	'R_47s_ROI'	'L_LIPd_ROI'	'R_LIPd_ROI'	'L_6a_ROI'	'R_6a_ROI'	'L_i6-8_ROI'	'R_i6-8_ROI'	'L_s6-8_ROI'	'R_s6-8_ROI'	'L_43_ROI'	'R_43_ROI'};

% Handle case when number of anatomical regions are not identical to atlas regions
actualRegions = {sScout.Scouts.Label};
missingRegions = setdiff(expectedRegions, actualRegions);

% Assuming sScout.Scouts is not empty and has at least one scout
if ~isempty(sScout.Scouts)
    % Identify all fields from the first scout as a template
    fieldNames = fieldnames(sScout.Scouts(1));
    % Prepare an empty scout template with all fields
    emptyScout = cell2struct(cell(length(fieldNames), 1), fieldNames, 1);
    
    % Default empty values for known fields
    emptyScout.Label = ''; % Update as necessary
    emptyScout.Vertices = [];
    emptyScout.Seed = 0; % Or any appropriate 'empty' value
    
    % Now, handle missing regions with this updated emptyScout
    for i = 1:length(missingRegions)
        emptyScout.Label = missingRegions{i};
        % Insert the empty scout at the correct position
        idx = find(strcmp(expectedRegions, missingRegions{i}));
        sScout.Scouts = [sScout.Scouts(1:idx-1), emptyScout, sScout.Scouts(idx:end)];
    end
else
    warning('sScout.Scouts is empty, cannot determine structure fields.');
end
end

function [RoiLabels, RoiIndices] = defineROIs_HCP(~)
% Define regions of interest (ROIs)

allRegions = {'L_V1_ROI'	'R_V1_ROI'	'L_FEF_ROI'	'R_FEF_ROI'	'L_OP4_ROI'	'R_OP4_ROI'	'L_OP1_ROI'	'R_OP1_ROI'	'L_OP2-3_ROI'	'R_OP2-3_ROI'	'L_52_ROI'	'R_52_ROI'	'L_RI_ROI'	'R_RI_ROI'	'L_PFcm_ROI'	'R_PFcm_ROI'	'L_PoI2_ROI'	'R_PoI2_ROI'	'L_TA2_ROI'	'R_TA2_ROI'	'L_FOP4_ROI'	'R_FOP4_ROI'	'L_MI_ROI'	'R_MI_ROI'	'L_PEF_ROI'	'R_PEF_ROI'	'L_Pir_ROI'	'R_Pir_ROI'	'L_AVI_ROI'	'R_AVI_ROI'	'L_AAIC_ROI'	'R_AAIC_ROI'	'L_FOP1_ROI'	'R_FOP1_ROI'	'L_FOP3_ROI'	'R_FOP3_ROI'	'L_FOP2_ROI'	'R_FOP2_ROI'	'L_PFt_ROI'	'R_PFt_ROI'	'L_AIP_ROI'	'R_AIP_ROI'	'L_EC_ROI'	'R_EC_ROI'	'L_PreS_ROI'	'R_PreS_ROI'	'L_55b_ROI'	'R_55b_ROI'	'L_H_ROI'	'R_H_ROI'	'L_ProS_ROI'	'R_ProS_ROI'	'L_PeEc_ROI'	'R_PeEc_ROI'	'L_STGa_ROI'	'R_STGa_ROI'	'L_PBelt_ROI'	'R_PBelt_ROI'	'L_A5_ROI'	'R_A5_ROI'	'L_PHA1_ROI'	'R_PHA1_ROI'	'L_PHA3_ROI'	'R_PHA3_ROI'	'L_STSda_ROI'	'R_STSda_ROI'	'L_STSdp_ROI'	'R_STSdp_ROI'	'L_V3A_ROI'	'R_V3A_ROI'	'L_STSvp_ROI'	'R_STSvp_ROI'	'L_TGd_ROI'	'R_TGd_ROI'	'L_TE1a_ROI'	'R_TE1a_ROI'	'L_TE1p_ROI'	'R_TE1p_ROI'	'L_TE2a_ROI'	'R_TE2a_ROI'	'L_TF_ROI'	'R_TF_ROI'	'L_TE2p_ROI'	'R_TE2p_ROI'	'L_PHT_ROI'	'R_PHT_ROI'	'L_PH_ROI'	'R_PH_ROI'	'L_TPOJ1_ROI'	'R_TPOJ1_ROI'	'L_RSC_ROI'	'R_RSC_ROI'	'L_TPOJ2_ROI'	'R_TPOJ2_ROI'	'L_TPOJ3_ROI'	'R_TPOJ3_ROI'	'L_DVT_ROI'	'R_DVT_ROI'	'L_PGp_ROI'	'R_PGp_ROI'	'L_IP2_ROI'	'R_IP2_ROI'	'L_IP1_ROI'	'R_IP1_ROI'	'L_IP0_ROI'	'R_IP0_ROI'	'L_PFop_ROI'	'R_PFop_ROI'	'L_PF_ROI'	'R_PF_ROI'	'L_PFm_ROI'	'R_PFm_ROI'	'L_POS2_ROI'	'R_POS2_ROI'	'L_PGi_ROI'	'R_PGi_ROI'	'L_PGs_ROI'	'R_PGs_ROI'	'L_V6A_ROI'	'R_V6A_ROI'	'L_VMV1_ROI'	'R_VMV1_ROI'	'L_VMV3_ROI'	'R_VMV3_ROI'	'L_PHA2_ROI'	'R_PHA2_ROI'	'L_V4t_ROI'	'R_V4t_ROI'	'L_FST_ROI'	'R_FST_ROI'	'L_V3CD_ROI'	'R_V3CD_ROI'	'L_LO3_ROI'	'R_LO3_ROI'	'L_V7_ROI'	'R_V7_ROI'	'L_VMV2_ROI'	'R_VMV2_ROI'	'L_31pd_ROI'	'R_31pd_ROI'	'L_31a_ROI'	'R_31a_ROI'	'L_VVC_ROI'	'R_VVC_ROI'	'L_25_ROI'	'R_25_ROI'	'L_s32_ROI'	'R_s32_ROI'	'L_pOFC_ROI'	'R_pOFC_ROI'	'L_PoI1_ROI'	'R_PoI1_ROI'	'L_Ig_ROI'	'R_Ig_ROI'	'L_FOP5_ROI'	'R_FOP5_ROI'	'L_IPS1_ROI'	'R_IPS1_ROI'	'L_p10p_ROI'	'R_p10p_ROI'	'L_p47r_ROI'	'R_p47r_ROI'	'L_TGv_ROI'	'R_TGv_ROI'	'L_MBelt_ROI'	'R_MBelt_ROI'	'L_LBelt_ROI'	'R_LBelt_ROI'	'L_A4_ROI'	'R_A4_ROI'	'L_STSva_ROI'	'R_STSva_ROI'	'L_TE1m_ROI'	'R_TE1m_ROI'	'L_PI_ROI'	'R_PI_ROI'	'L_a32pr_ROI'	'R_a32pr_ROI'	'L_FFC_ROI'	'R_FFC_ROI'	'L_p24_ROI'	'R_p24_ROI'	'L_V3B_ROI'	'R_V3B_ROI'	'L_MST_ROI'	'R_MST_ROI'	'L_LO1_ROI'	'R_LO1_ROI'	'L_LO2_ROI'	'R_LO2_ROI'	'L_PIT_ROI'	'R_PIT_ROI'	'L_MT_ROI'	'R_MT_ROI'	'L_A1_ROI'	'R_A1_ROI'	'L_PSL_ROI'	'R_PSL_ROI'	'L_SFL_ROI'	'R_SFL_ROI'	'L_PCV_ROI'	'R_PCV_ROI'	'L_STV_ROI'	'R_STV_ROI'	'L_7Pm_ROI'	'R_7Pm_ROI'	'L_V6_ROI'	'R_V6_ROI'	'L_7m_ROI'	'R_7m_ROI'	'L_POS1_ROI'	'R_POS1_ROI'	'L_23d_ROI'	'R_23d_ROI'	'L_v23ab_ROI'	'R_v23ab_ROI'	'L_d23ab_ROI'	'R_d23ab_ROI'	'L_31pv_ROI'	'R_31pv_ROI'	'L_5m_ROI'	'R_5m_ROI'	'L_5mv_ROI'	'R_5mv_ROI'	'L_23c_ROI'	'R_23c_ROI'	'L_5L_ROI'	'R_5L_ROI'	'L_V2_ROI'	'R_V2_ROI'	'L_24dd_ROI'	'R_24dd_ROI'	'L_24dv_ROI'	'R_24dv_ROI'	'L_7AL_ROI'	'R_7AL_ROI'	'L_SCEF_ROI'	'R_SCEF_ROI'	'L_6ma_ROI'	'R_6ma_ROI'	'L_7Am_ROI'	'R_7Am_ROI'	'L_7PL_ROI'	'R_7PL_ROI'	'L_7PC_ROI'	'R_7PC_ROI'	'L_LIPv_ROI'	'R_LIPv_ROI'	'L_VIP_ROI'	'R_VIP_ROI'	'L_V3_ROI'	'R_V3_ROI'	'L_MIP_ROI'	'R_MIP_ROI'	'L_1_ROI'	'R_1_ROI'	'L_2_ROI'	'R_2_ROI'	'L_3a_ROI'	'R_3a_ROI'	'L_6d_ROI'	'R_6d_ROI'	'L_6mp_ROI'	'R_6mp_ROI'	'L_6v_ROI'	'R_6v_ROI'	'L_p24pr_ROI'	'R_p24pr_ROI'	'L_33pr_ROI'	'R_33pr_ROI'	'L_a24pr_ROI'	'R_a24pr_ROI'	'L_V4_ROI'	'R_V4_ROI'	'L_p32pr_ROI'	'R_p32pr_ROI'	'L_a24_ROI'	'R_a24_ROI'	'L_d32_ROI'	'R_d32_ROI'	'L_8BM_ROI'	'R_8BM_ROI'	'L_p32_ROI'	'R_p32_ROI'	'L_10r_ROI'	'R_10r_ROI'	'L_47m_ROI'	'R_47m_ROI'	'L_8Av_ROI'	'R_8Av_ROI'	'L_8Ad_ROI'	'R_8Ad_ROI'	'L_9m_ROI'	'R_9m_ROI'	'L_V8_ROI'	'R_V8_ROI'	'L_8BL_ROI'	'R_8BL_ROI'	'L_9p_ROI'	'R_9p_ROI'	'L_10d_ROI'	'R_10d_ROI'	'L_8C_ROI'	'R_8C_ROI'	'L_44_ROI'	'R_44_ROI'	'L_45_ROI'	'R_45_ROI'	'L_47l_ROI'	'R_47l_ROI'	'L_a47r_ROI'	'R_a47r_ROI'	'L_6r_ROI'	'R_6r_ROI'	'L_IFJa_ROI'	'R_IFJa_ROI'	'L_4_ROI'	'R_4_ROI'	'L_IFJp_ROI'	'R_IFJp_ROI'	'L_IFSp_ROI'	'R_IFSp_ROI'	'L_IFSa_ROI'	'R_IFSa_ROI'	'L_p9-46v_ROI'	'R_p9-46v_ROI'	'L_46_ROI'	'R_46_ROI'	'L_a9-46v_ROI'	'R_a9-46v_ROI'	'L_9-46d_ROI'	'R_9-46d_ROI'	'L_9a_ROI'	'R_9a_ROI'	'L_10v_ROI'	'R_10v_ROI'	'L_a10p_ROI'	'R_a10p_ROI'	'L_3b_ROI'	'R_3b_ROI'	'L_10pp_ROI'	'R_10pp_ROI'	'L_11l_ROI'	'R_11l_ROI'	'L_13l_ROI'	'R_13l_ROI'	'L_OFC_ROI'	'R_OFC_ROI'	'L_47s_ROI'	'R_47s_ROI'	'L_LIPd_ROI'	'R_LIPd_ROI'	'L_6a_ROI'	'R_6a_ROI'	'L_i6-8_ROI'	'R_i6-8_ROI'	'L_s6-8_ROI'	'R_s6-8_ROI'	'L_43_ROI'	'R_43_ROI'};

AngROIs = {'V7_ROI','IPS1_ROI', 'TPOJ3_ROI', 'PGp_ROI','IP1_ROI','PGi_ROI','PGs_ROI', 'V6A_ROI'};

% Define the new variable as a cell array of ROI names
FrontROI = {
    '10d_ROI', '10r_ROI', '10v_ROI', '11l_ROI', '13l_ROI', '23d_ROI', '33pr_ROI', '44_ROI', '45_ROI', ...
    '46_ROI', '47l_ROI', '47m_ROI', '47s_ROI', '55b_ROI', '8Ad_ROI', '8Av_ROI', '8BL_ROI', '8BM_ROI', ...
    '8C_ROI', '9-46d_ROI', '9a_ROI', '9m_ROI', '9p_ROI', 'AVI_ROI', 'FOP5_ROI', 'IFJa_ROI', 'IFJp_ROI', ...
    'IFSa_ROI', 'IFSp_ROI', 'SFL_ROI', 'a10p_ROI', 'a32pr_ROI', 'a47r_ROI', 'a9-46v_ROI', 'd32_ROI', ...
    'i6-8_ROI', 'p10p_ROI', 'p47r_ROI', 'p9-46v_ROI', 's32_ROI', 's6-8_ROI'
    };

TempROI = {
    'L_FFC_ROI', 'L_EC_ROI', 'L_PreS_ROI', 'L_H_ROI', 'L_PeEc_ROI', 'L_STGa_ROI', 'L_A5_ROI', ...
    'L_PHA1_ROI', 'L_PHA3_ROI', 'L_STSda_ROI', 'L_STSdp_ROI', 'L_STSvp_ROI', 'L_TGd_ROI', ...
    'L_TE1a_ROI', 'L_TE1p_ROI', 'L_TE2a_ROI', 'L_TF_ROI', 'L_TE2p_ROI', 'L_PHT_ROI', 'L_PH_ROI', ...
    'L_PHA2_ROI', 'L_VVC_ROI', 'L_TGv_ROI', 'L_STSva_ROI', 'L_TE1m_ROI'
    };

% Optionally, create a new variable with 'L_' removed
TempROI = cellfun(@(x) x(3:end), TempROI, 'UniformOutput', false);

LatROIs = {
    'L_PEF_ROI', 'L_V7_ROI', 'L_IPS1_ROI', 'L_7PL_ROI', 'L_MIP_ROI', 'L_47m_ROI', 'L_8Av_ROI', 'L_8C_ROI', ...
    'L_44_ROI', 'L_45_ROI', 'L_47l_ROI', 'L_a47r_ROI', 'L_IFJa_ROI', 'L_IFJp_ROI', 'L_IFSp_ROI', 'L_IFSa_ROI', ...
    'L_p9-46v_ROI', 'L_13l_ROI', 'L_47s_ROI', 'L_i6-8_ROI', 'L_AVI_ROI', 'L_AAIC_ROI', 'L_STGa_ROI', 'L_A5_ROI', ...
    'L_STSda_ROI', 'L_STSdp_ROI', 'L_STSvp_ROI', 'L_TGd_ROI', 'L_TE1a_ROI', 'L_TE1p_ROI', 'L_TE2a_ROI', 'L_TE2p_ROI', ...
    'L_PHT_ROI', 'L_PH_ROI', 'L_TPOJ3_ROI', 'L_PGp_ROI', 'L_IP1_ROI', 'L_IP0_ROI', 'L_PGi_ROI', 'L_PGs_ROI', ...
    'L_FOP5_ROI', 'L_p47r_ROI', 'L_TGv_ROI', 'L_STSva_ROI', 'L_TE1m_ROI'
    };

% Remove 'L_' prefix from each ROI
LatROIs = cellfun(@(x) x(3:end), LatROIs, 'UniformOutput', false);

% Combine all ROIs and their labels
RoiLabels = {'Angular', 'Frontal', 'Temporal', 'Lateral'};
RoiGroups = {AngROIs, FrontROI, TempROI, LatROIs};

RoiIndices_L = cell(length(RoiGroups), 1);
for i = 1:length(RoiGroups)
    % Prepend 'L_' to each ROI name in the current group if it does not already start with 'L_'
    modifiedROIs = cellfun(@(roi) ['L_' roi], RoiGroups{i}, 'UniformOutput', false);
    
    % Remove any redundant 'L_' if it already exists
    modifiedROIs = strrep(modifiedROIs, 'L_L_', 'L_');
    
    % Find indices of ROIs in allRegions
    [isPresent, index] = ismember(modifiedROIs, allRegions);
    RoiIndices_L{i} = index;
    
    % Collect names of ROIs that were not detected
    NotDetectedROIs{i} = modifiedROIs(~isPresent);
end

RoiIndices_R = cell(length(RoiGroups), 1);
for i = 1:length(RoiGroups)
    % Prepend 'L_' to each ROI name in the current group if it does not already start with 'L_'
    modifiedROIs = cellfun(@(roi) ['R_' roi], RoiGroups{i}, 'UniformOutput', false);
    
    % Remove any redundant 'L_' if it already exists
    modifiedROIs = strrep(modifiedROIs, 'R_R_', 'R_');
    
    %     RoiIndices_R{i} = find(ismember(allRegions, modifiedROIs));
    % Find indices of ROIs in allRegions
    [isPresent, index] = ismember(modifiedROIs, allRegions);
    RoiIndices_R{i} = index;
    
    % Collect names of ROIs that were not detected
    NotDetectedROIs{i} = modifiedROIs(~isPresent);
end

% Assuming RoiIndices_L and RoiIndices_R are already defined as described
% Initialize the merged cell array
RoiIndices = cell(size(RoiIndices_L));

% Loop through each cell to merge the corresponding L and R indices
for i = 1:length(RoiIndices_L)
    % Combine the indices from left and right ROIs
    % This assumes that the order and number of groups in L and R are the same
    RoiIndices{i} = unique([RoiIndices_L{i}, RoiIndices_R{i}]);
end

end

function [Summ_LI, Summ_CI, LI_label_out, L_vertices_total, R_vertices_total, CI_strings, CI_widths]  = computeLI_bootstrap(cfg_LI)
% Computes the Laterality Index (LI) using bootstrapping and exports the results, including vertex counts.

% Perform bootstrapping for each ROI
RoiLabels = cfg_LI.RoiLabels;
TotROI = length(cfg_LI.RoiIndices);
Summ_LI = zeros(1, TotROI); % Initialize the vector for final LIs
L_vertices_total = zeros(1, TotROI); % Initialize vectors for left vertices count
R_vertices_total = zeros(1, TotROI); % Initialize vectors for right vertices count
LI_label_out = cell(1, TotROI);

% Pre-allocate output arrays
Summ_LI = zeros(1, TotROI);
Summ_CI = zeros(TotROI, 2);  % CI is assumed to have two values: lower and upper.
L_vertices_total = zeros(1, TotROI);
R_vertices_total = zeros(1, TotROI);
LI_label_out = cell(1, TotROI);

for ii = 1:TotROI
    
    % Set up configuration for bootstrapping
    cfg_main = [];
    cfg_main.atlas = cfg_LI.sScout;
    cfg_main.RoiIndices = cfg_LI.RoiIndices{ii};
    cfg_main.divs = cfg_LI.divs;
    cfg_main.n_resampling = cfg_LI.n_resampling;
    cfg_main.RESAMPLE_RATIO = cfg_LI.RESAMPLE_RATIO;
    cfg_main.t1 = cfg_LI.t1;
    cfg_main.t2 = cfg_LI.t2;
    cfg_main.ImageGridAmp = cfg_LI.ImageGridAmp;
    cfg_main.Tinterval = cfg_LI.Tinterval;
    
    % Call bootstrapping function for the current ROI
    [weighted_li, ~, L_vertices_above_thresh, R_vertices_above_thresh, CI] = do_LI_bootstrap(cfg_main);
    
    % Store results
    Summ_LI(ii) = weighted_li;
    Summ_CI(ii, :) = CI;  % Store both lower and upper CI values
    L_vertices_total(ii) = sum(L_vertices_above_thresh);
    R_vertices_total(ii) = sum(R_vertices_above_thresh);
    LI_label_out{ii} = RoiLabels{ii};
end

% Pre-compute a formatted CI string and CI width for each ROI
TotROI = length(LI_label_out);
CI_strings = cell(TotROI,1);
CI_widths = zeros(TotROI,1);
for i = 1:TotROI
    lowerCI = Summ_CI(i,1);
    upperCI = Summ_CI(i,2);
    CI_strings{i} = sprintf('[%.2f, %.2f]', lowerCI, upperCI);
    CI_widths(i) = upperCI - lowerCI;
end


if cfg_LI.report == 1
    
    % Save or display results
    savedir = cfg_LI.savedir;   % Directory to save results
    sname = cfg_LI.sname;       % Use sname from cfg_LI for the filename
    filename = fullfile(savedir, sname);
    
    % Open file for writing
    fid = fopen(filename, 'w');
    % Updated header: now includes CI and CI width
    fprintf(fid, 'ROI\tLI\tCI\tCI_Width\tL_Vertices\tR_Vertices\n');
    for i = 1:TotROI
        fprintf(fid, '%s\t%f\t%s\t%f\t%d\t%d\n', ...
            LI_label_out{i}, Summ_LI(i), CI_strings{i}, CI_widths(i), L_vertices_total(i), R_vertices_total(i));
    end
    fclose(fid);
    
    disp('Results saved to: ');
    disp(filename);
    
    % Before creating the table, combine L and R vertices into one column
    LR_vertices_str = cell(TotROI, 1);
    for i = 1:TotROI
        LR_vertices_str{i} = sprintf('%d-%d', L_vertices_total(i), R_vertices_total(i));
    end
    
    % Convert each to a column vector (if needed)
    Summ_LI = Summ_LI(:);
    RoiLabels = RoiLabels(:);
    CI_strings = CI_strings(:);
    CI_widths = CI_widths(:);
    LR_vertices_str = LR_vertices_str(:);
    
    % Now create the table with the combined LR column
    T = table(RoiLabels, Summ_LI, CI_strings, CI_widths, LR_vertices_str, ...
        'VariableNames', {'ROI', 'LI', 'CI_95', 'CI_Width', 'Vertices_LR'});
    
    disp(T);
    
end
end

function computeLI(cfg_LI)
% Compute the Laterality Index (LI) and associated tasks

Tinterval = cfg_LI.Tinterval;
timerange = cfg_LI.timerange;
RoiLabels = cfg_LI.RoiLabels;
RoiIndices  = cfg_LI.RoiIndices;
sScout  = cfg_LI.sScout;
Ratio4Threshold = cfg_LI.Ratio4Threshold;
savedir = cfg_LI.savedir;
windows = cfg_LI.windows;

% Adjustments for dimensions of verticies when concatenating
for i = 1:length(sScout.Scouts)
    % Check if Vertices is nx1, transpose only in this case
    if size(sScout.Scouts(i).Vertices, 1) > 1
        sScout.Scouts(i).Vertices = sScout.Scouts(i).Vertices';
    end
end

%%
if isempty(windows) && cfg_LI.method ~= 3
    
    [Summ_LI, LI_label_out, L_count, R_count, threshold] = computeLI_Svalue_Count(cfg_LI);
    Summ_LI = round(Summ_LI, 2); L_count = round(L_count, 2); R_count = round(R_count, 2);
    
    cfg_LI.Summ_LI = Summ_LI;
    cfg_LI.L_count = L_count;
    cfg_LI.R_count = R_count;
    cfg_LI.threshold = threshold;
    cfg_LI.LI_label_out = LI_label_out;
    
    if cfg_LI.saveResults == 1
        export_LI_results(cfg_LI);
    end
    report_LI(cfg_LI)
    
elseif isempty(windows) && cfg_LI.method == 3
    
    cfg_LI.report = 0;
    [Summ_LI, Summ_CI, LI_label_out, L_count, R_count, CI_strings, CI_widths]  = computeLI_bootstrap(cfg_LI);
    Summ_LI = round(Summ_LI, 2); L_count = round(L_count, 2); R_count = round(R_count, 2);

    cfg_LI.Summ_LI = Summ_LI;
    cfg_LI.Summ_CI = Summ_CI;
    cfg_LI.L_count = L_count;
    cfg_LI.R_count = R_count;
    cfg_LI.R_count = R_count;
    cfg_LI.CI_strings = CI_strings;
    cfg_LI.CI_widths = CI_widths;
    cfg_LI.LI_label_out = LI_label_out;
    
    if cfg_LI.saveResults == 1
        export_LI_results(cfg_LI);
    end
    report_LI(cfg_LI)
    
elseif ~isempty(windows) && cfg_LI.Tinterval == 3
    
    nWindows = size(windows, 1);
    nRois = numel(RoiIndices);
    final_LI = zeros(nWindows, nRois);
    
    if cfg_LI.method == 3
        % The bootstrap path computes its thresholds internally and does not
        % use globmax_rois, so avoid the expensive global-max preprocessing.
        final_CI = zeros(nWindows, nRois, 2);
    else
        final_CI = [];
        cfg_LI.globmax_rois = compute_globmax_rois(cfg_LI);
    end
    
    for j = 1:nWindows
        [~, cfg_LI.t1] = min(abs(cfg_LI.Time - windows(j,1)));
        [~, cfg_LI.t2] = min(abs(cfg_LI.Time - windows(j,2)));
        if cfg_LI.method == 1 || cfg_LI.method == 2
            [Summ_LI, ~, ~, ~, ~] = computeLI_Svalue_Count(cfg_LI);
        elseif cfg_LI.method == 3
            cfg_LI.report = 0;
            progressEvery = max(1, ceil(nWindows / 20));
            if j == 1 || j == nWindows || mod(j-1, progressEvery) == 0
                fprintf('Bootstrap window %d/%d: [%.4f, %.4f] s\n', ...
                    j, nWindows, windows(j,1), windows(j,2));
            end
            [Summ_LI, Summ_CI, ~] = computeLI_bootstrap(cfg_LI);
            final_CI(j,:,:) = reshape(Summ_CI, [1, nRois, 2]);
        end
        final_LI(j,:) = Summ_LI;
    end
    
    cfg_LI.final_LI = final_LI;
    if cfg_LI.method == 3, cfg_LI.final_CI = final_CI; end
    
    if cfg_LI.plotResults == 1
        plot_LI(cfg_LI);
    end
    LI_summary = report_tLI(cfg_LI);
    
    cfg_LI.LI_summary = LI_summary;
    
    if cfg_LI.saveResults == 1
        export_LI_results(cfg_LI);
    end
    
end
end

function [weighted_li, num_threshvals, L_vertices_above_thresh, R_vertices_above_thresh, CI] = do_LI_bootstrap(cfg_main)
divs = cfg_main.divs;
n_resampling = cfg_main.n_resampling;
RESAMPLE_RATIO = cfg_main.RESAMPLE_RATIO;
RoiIndices = cfg_main.RoiIndices;
MIN_NUM_THRESH_VOXELS = divs / RESAMPLE_RATIO; % Adjust as needed

% For the left scouts (RoiIndices(1), (3), (5), ...):
LHscoutCell = arrayfun(@(x) cfg_main.atlas.Scouts(x).Vertices(:), ...
                       RoiIndices(1:2:end), ...
                       'UniformOutput', false);
LHscout = vertcat(LHscoutCell{:});

% For the right scouts (RoiIndices(2), (4), (6), ...):
RHscoutCell = arrayfun(@(x) cfg_main.atlas.Scouts(x).Vertices(:), ...
                       RoiIndices(2:2:end), ...
                       'UniformOutput', false);
RHscout = vertcat(RHscoutCell{:});

% Extract only the ROI data. The previous implementation copied the
% full-brain window once per ROI and then discarded nearly all of it.
sourceData = cfg_main.ImageGridAmp;
if cfg_main.Tinterval == 2 || size(sourceData,2) == 1
    LHvals = mean(sourceData(LHscout,:), 2);
    RHvals = mean(sourceData(RHscout,:), 2);
else
    time_idx = cfg_main.t1:cfg_main.t2;
    if numel(time_idx) > 500
        time_idx = time_idx(1:10:end);
    end
    LHvals = sourceData(LHscout, time_idx);
    RHvals = sourceData(RHscout, time_idx);
end

ROIMax = max([LHvals(:); RHvals(:)]);
threshvals = linspace(0, ROIMax, divs);

weighted_li = 0;
num_threshvals = 0;

cumulative_L_vertices_above_thresh = 0;
cumulative_R_vertices_above_thresh = 0;

% Preallocate all bootstrap LI distributions across thresholds.
All_TB_LIs_weighted = zeros(n_resampling, numel(threshvals));

for thresh_idx = 1:numel(threshvals)
    thresh = threshvals(thresh_idx);
    l_above_thresh = LHvals >= thresh;
    r_above_thresh = RHvals >= thresh;
    
    l_count = sum(l_above_thresh(:));
    r_count = sum(r_above_thresh(:));
    if l_count < MIN_NUM_THRESH_VOXELS && r_count < MIN_NUM_THRESH_VOXELS
        break; % Not enough voxels above threshold
    end
    
    cumulative_L_vertices_above_thresh = cumulative_L_vertices_above_thresh + l_count;
    cumulative_R_vertices_above_thresh = cumulative_R_vertices_above_thresh + r_count;
    
    % Compute bootstrap distribution for this threshold
    TB_LIs = bootstrapLI(l_above_thresh, r_above_thresh, n_resampling, RESAMPLE_RATIO);
    
    % Accumulate weighted LI
    weighted_li = weighted_li + mean(TB_LIs) * thresh_idx;
    num_threshvals = thresh_idx;
    
    % Store the weighted TB_LIs for CI calculation
    % We store TB_LIs * thresh_idx so that later we can sum them all and then divide by sum(1:num_threshvals).
    All_TB_LIs_weighted(:, thresh_idx) = TB_LIs * thresh_idx;
end

if num_threshvals > 0
    weighted_li = (weighted_li / sum(1:num_threshvals)) * 100;
end

L_vertices_above_thresh = round(cumulative_L_vertices_above_thresh / max(num_threshvals,1));
R_vertices_above_thresh = round(cumulative_R_vertices_above_thresh / max(num_threshvals,1));

% Compute the combined bootstrap distribution for the final weighted LI
if num_threshvals > 0
    CombinedLI = (sum(All_TB_LIs_weighted(:,1:num_threshvals), 2) / sum(1:num_threshvals)) * 100;
    % Calculate 95% CI from CombinedLI
    CI = prctile(CombinedLI, [2.5 97.5]);
else
    CI = [NaN NaN]; % No CI if no thresholds were processed
end
end

function TB_LIs = bootstrapLI(Lvals, Rvals, n_samples, resample_ratio)
% Bootstrap LI for logical suprathreshold observations.
%
% Flatten explicitly so that the bootstrap population consists of individual
% vertex-time observations. MATLAB datasample samples rows from a matrix by
% default, which previously created extremely large k-by-time arrays.
Lvals = Lvals(:);
Rvals = Rvals(:);

if isempty(Lvals) || isempty(Rvals)
    TB_LIs = NaN(n_samples, 1);
    return;
end

l_n = max(round(numel(Lvals) * resample_ratio), 1);
r_n = max(round(numel(Rvals) * resample_ratio), 1);

% For a logical population sampled with replacement, the number of true
% observations in each sample is binomial. This is distributionally
% equivalent to explicitly resampling every logical element, without
% allocating the large temporary sample vectors.
l_means = binornd(l_n, mean(Lvals), [n_samples, 1]) ./ l_n;
r_means = binornd(r_n, mean(Rvals), [n_samples, 1]) ./ r_n;

denominator = l_means + r_means;
TB_LIs = (l_means - r_means) ./ denominator;
TB_LIs(denominator == 0) = NaN;
end

function wi = generateTimeWindows(cfg_main)

strt = cfg_main.strt;
spt = cfg_main.spt;
overlap = cfg_main.overlap;       % Fraction in the range [0, 1)
linterval = cfg_main.linterval;   % Window length in seconds

if linterval <= 0
    error('Window length must be greater than zero.');
end
if overlap < 0 || overlap >= 1
    error('Window overlap must be at least 0 and less than 100 percent.');
end

duration = spt - strt;

% Allow for binary floating-point roundoff at exact boundaries such as
% 0.700 - 0.400 versus a requested 0.300-second window.
timeScale = max([1, abs(strt), abs(spt), abs(linterval)]);
timeTolerance = 100 * eps(timeScale);
if duration + timeTolerance < linterval
    wi = zeros(0, 2);
    return;
end

step = linterval * (1 - overlap);
nWindows = floor((duration - linterval + timeTolerance) / step) + 1;
starts = strt + (0:nWindows-1)' * step;
stops = starts + linterval;
stops(abs(stops - spt) <= timeTolerance) = spt;
wi = [starts, stops];

end

function [Summ_LI, LI_label_out, L_count, R_count, threshold] = computeLI_Svalue_Count(cfg_LI)
% computeLI_Svalue_Count Computes the Laterality Index (LI) for each ROI defined in cfg_LI.
%
% INPUT:
%   cfg_LI: A struct containing configuration and data.
%       .RoiIndices      : Cell array, each cell contains indices of vertices for a given ROI (left/right).
%       .RoiLabels       : Cell array of ROI names (strings).
%       .sScout          : Struct with field .Scouts, each containing a list of vertices for a subregion.
%       .ImageGridAmp    : Matrix of source amplitudes (Vertices x Time).
%       .Threshtype      : Integer (1,2,3) specifying how threshold is determined.
%       .Ratio4Threshold : Float ratio (0-1) used to compute threshold.
%       .GlobalMax       : Max amplitude over all vertices and times.
%       .AllMax          : Max amplitude in the time interval of interest.
%       .method          : Integer specifying LI calculation method (e.g. 1=Power-based, 2=Count-based).
%       .t1, .t2         : Time indices defining the interval of interest.
%       .Tinterval       : Interval type (1,2,3).
%       (Optional) .globmax_rois : If Tinterval=3, per-ROI max values.
%
% OUTPUT:
%   Summ_LI     : Array of computed LI values for each ROI.
%   LI_label_out: Cell array of ROI labels corresponding to Summ_LI.
%   L_count     : Array of left hemisphere counts (or power) used in LI computation.
%   R_count     : Array of right hemisphere counts (or power).
%   threshold   : Computed threshold value used for LI calculation.

% Initialize output variables
RoiIndices = cfg_LI.RoiIndices;
RoiLabels = cfg_LI.RoiLabels;
sScout = cfg_LI.sScout;
ImageGridAmp = cfg_LI.ImageGridAmp;
Threshtype = cfg_LI.Threshtype;
Ratio4Threshold = cfg_LI.Ratio4Threshold;
GlobalMax = cfg_LI.GlobalMax;
AllMax = cfg_LI.AllMax;
method = cfg_LI.method;
Summ_LI = zeros(1, length(RoiIndices));
LI_label_out = cell(1, length(RoiIndices));

t1 = cfg_LI.t1;
t2 = cfg_LI.t2;

L_count = zeros(1, length(RoiIndices));
R_count = zeros(1, length(RoiIndices));

% Process each ROI
for ii = 1:length(RoiIndices)
    s2 = RoiLabels{ii};
    %     hemi_roi_num = length(RoiIndices{ii});
    curr_subregion = sScout.Scouts(RoiIndices{ii});
    
    % Split indices into left and right Rois
    Ltemp_region = [];
    for iSub = 1 : 2 : length(curr_subregion)
        % Turn each .Vertices into a column, then append
        Ltemp_region = [Ltemp_region; curr_subregion(iSub).Vertices(:)];
    end
    
    Rtemp_region = [];
    for iSub = 2 : 2 : length(curr_subregion)
        Rtemp_region = [Rtemp_region; curr_subregion(iSub).Vertices(:)];
    end    
    % Get values for the current ROI
    switch cfg_LI.Tinterval % modified by VY
        case 2
            %First parse the maps into separate space-times maps for each side
            LHvals = ImageGridAmp(Ltemp_region, :);
            LH_max = max(max(LHvals));
            RHvals = ImageGridAmp(Rtemp_region, :);
            RH_max = max(max(RHvals));
            ROIMax = max(LH_max,RH_max);
        otherwise
            %First parse the maps into separate space-times maps for each side
            LHvals = ImageGridAmp(Ltemp_region,t1:t2);
            LH_max = max(max(LHvals));
            RHvals = ImageGridAmp(Rtemp_region,t1:t2);
            RH_max = max(max(RHvals));
            ROIMax = max(LH_max,RH_max);
    end
    
    % For window-based approach (Tinterval == 3), use the precomputed global max for each ROI:
    if cfg_LI.Tinterval == 3
        ROIMax = cfg_LI.globmax_rois(ii);
    end
    
    % Select threshold based on type
    switch Threshtype
        case 1
            threshold = Ratio4Threshold * GlobalMax;
        case 2
            threshold = Ratio4Threshold * AllMax;
        case 3
            threshold = Ratio4Threshold * ROIMax;
    end
    
    switch method
        case 1
            % Power-based LI calculation
            pow_left  = sum(LHvals(LHvals(:) > threshold));
            pow_left = pow_left/size(LHvals,1);
            
            pow_right = sum(RHvals(RHvals(:) > threshold));
            pow_right = pow_right/size(RHvals,1);
            
            LI_ROIval = 100 * ((pow_left - pow_right) / (pow_left + pow_right));
            
            L_count(ii) = pow_left;
            R_count(ii) = pow_right;
            
            % Store the results
            Summ_LI(ii) = LI_ROIval;
            
        otherwise
            % Counting-based LI calculation
            L_ROIcount = sum(LHvals(:) > threshold);
            R_ROIcount = sum(RHvals(:) > threshold);
            
            LI_ROIcount = 100 * ((L_ROIcount - R_ROIcount) / (L_ROIcount + R_ROIcount));
            
            L_count(ii) = L_ROIcount;
            R_count(ii) = R_ROIcount;
            
            % Store the results
            Summ_LI(ii) = LI_ROIcount;
    end
    LI_label_out{ii} = s2;
end
end

function plot_LI(cfg_LI)
% PLOT_LI Plot window-based LI time courses using the actual window centers.

if ~isfield(cfg_LI, 'windows') || isempty(cfg_LI.windows)
    warning('No time windows are available for plotting.');
    return;
end
if ~isfield(cfg_LI, 'final_LI') || isempty(cfg_LI.final_LI)
    warning('No LI values are available for plotting.');
    return;
end

windows = cfg_LI.windows;
LI = cfg_LI.final_LI;
[W, R] = size(LI);

if size(windows,1) ~= W || size(windows,2) ~= 2
    warning('The LI matrix and time-window matrix have incompatible dimensions.');
    return;
end

% Plot each estimate at the center of its analysis window.
timeCenters = mean(windows, 2);
xLeft = min(windows(:,1));
xRight = max(windows(:,2));
if ~isfinite(xLeft) || ~isfinite(xRight)
    warning('The time-window limits contain non-finite values.');
    return;
end
if xRight <= xLeft
    xLeft = timeCenters(1) - 0.05;
    xRight = timeCenters(1) + 0.05;
end

switch cfg_LI.method
    case 1
        methodTitle = 'LI over time - Source magnitude';
    case 2
        methodTitle = 'LI over time - Vertex count';
    case 3
        methodTitle = 'LI over time - Bootstrap';
    otherwise
        methodTitle = 'Lateralization index over time';
end

windowLengthMs = 1000 * median(windows(:,2) - windows(:,1));
if W > 1
    windowStepMs = 1000 * median(diff(windows(:,1)));
    overlapPct = 100 * (1 - windowStepMs / windowLengthMs);
    overlapPct = min(max(overlapPct, 0), 99);
    detailTitle = sprintf('%d windows | %.0f ms length | %.0f%% overlap', ...
        W, windowLengthMs, overlapPct);
else
    detailTitle = sprintf('1 window | %.0f ms length', windowLengthMs);
end

hFig = figure( ...
    'Color', 'w', ...
    'Name', methodTitle, ...
    'NumberTitle', 'off', ...
    'Position', [100, 100, 1080, 540], ...
    'InvertHardcopy', 'off');
ax = axes('Parent', hFig);
hold(ax, 'on');

colors = lines(max(R, 7));
colors = colors(1:R,:);
lineStyles = {'-', '--', '-.', ':'};

% LI is expressed as a percentage and is bounded by -100 and +100.
xlim(ax, [xLeft, xRight]);
ylim(ax, [-100, 100]);

% Bootstrap confidence intervals are drawn first so curves remain prominent.
hasCI = cfg_LI.method == 3 && ...
    isfield(cfg_LI, 'final_CI') && ~isempty(cfg_LI.final_CI) && ...
    size(cfg_LI.final_CI,1) == W && size(cfg_LI.final_CI,2) >= R && ...
    size(cfg_LI.final_CI,3) >= 2;

if hasCI
    for r = 1:R
        lowerCI = reshape(cfg_LI.final_CI(:,r,1), [], 1);
        upperCI = reshape(cfg_LI.final_CI(:,r,2), [], 1);
        valid = isfinite(timeCenters) & isfinite(lowerCI) & isfinite(upperCI);
        if nnz(valid) >= 2
            fill(ax, ...
                [timeCenters(valid); flipud(timeCenters(valid))], ...
                [lowerCI(valid); flipud(upperCI(valid))], ...
                colors(r,:), ...
                'FaceAlpha', 0.13, ...
                'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
    end
end

% Reference lines: no lateralization and, when visible, stimulus time zero.
line(ax, [xLeft, xRight], [0, 0], ...
    'Color', [0.25, 0.25, 0.25], ...
    'LineStyle', '-', ...
    'LineWidth', 0.9, ...
    'HandleVisibility', 'off');
if xLeft < 0 && xRight > 0
    line(ax, [0, 0], [-100, 100], ...
        'Color', [0.45, 0.45, 0.45], ...
        'LineStyle', ':', ...
        'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
end

hLines = gobjects(R, 1);
for r = 1:R
    hLines(r) = plot(ax, timeCenters, LI(:,r), ...
        'Color', colors(r,:), ...
        'LineStyle', lineStyles{mod(r-1, numel(lineStyles)) + 1}, ...
        'LineWidth', 2.0);
end

xlabel(ax, 'Time relative to event (s)');
ylabel(ax, {'Lateralization index (LI, %)', 'Positive = left | Negative = right'});
title(ax, {methodTitle, detailTitle}, 'FontWeight', 'bold');

set(ax, ...
    'FontSize', 10, ...
    'LineWidth', 1.0, ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'Layer', 'top', ...
    'XGrid', 'on', ...
    'YGrid', 'on', ...
    'GridAlpha', 0.16);

if isfield(cfg_LI, 'RoiLabels') && ~isempty(cfg_LI.RoiLabels)
    roiLabels = cfg_LI.RoiLabels(:);
    if numel(roiLabels) < R
        missingLabels = arrayfun(@(x) sprintf('ROI %d', x), ...
            numel(roiLabels)+1:R, 'UniformOutput', false)';
        roiLabels = [roiLabels; missingLabels];
    end
else
    roiLabels = arrayfun(@(x) sprintf('ROI %d', x), 1:R, ...
        'UniformOutput', false)';
end

legend(ax, hLines, roiLabels(1:R), ...
    'Location', 'eastoutside', ...
    'Interpreter', 'none', ...
    'Box', 'off');

if hasCI
    text(ax, 0.99, 0.02, 'Shading: 95% bootstrap CI', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', ...
        'Color', [0.35, 0.35, 0.35], ...
        'FontSize', 9);
end

hold(ax, 'off');
end

function final_table = report_tLI(cfg_LI)
% REPORT_TLI evaluates each ROI separately, finding the max LI value, the corresponding
% time interval, the median LI, and optionally the 95% CI at the max LI time point if
% bootstrapping was used.

[W, R] = size(cfg_LI.final_LI);

final_table = table();

for r = 1:R
    % Extract LI values for the current ROI
    roi_li = cfg_LI.final_LI(:, r);
    
    % Find max LI and corresponding window index
    [mx, max_win] = max(roi_li);
    md = nanmedian(roi_li);
    
    % Extract the interval for the max LI window
    max_interval = cfg_LI.windows(max_win, :);
    interval_str = sprintf('[%.2f, %.2f]', max_interval(1), max_interval(2));
    
    % Basic columns: ROI, Max_LI, Median_LI, and merged interval column
    ROI_label = cfg_LI.RoiLabels{r};
    a = table({ROI_label}, 'VariableNames', {'ROI'});
    b = table(mx, 'VariableNames', {'Max_LI'});
    d = table(md, 'VariableNames', {'Median_LI'});
    c = table({interval_str}, 'VariableNames', {'Max_Interval'});
    
    new_row = [a, b, d, c];
    
    % If bootstrapping was used, add CI information for this ROI at max_win
    if cfg_LI.method == 3
        roi_ci = squeeze(cfg_LI.final_CI(max_win, r, :));
        ci_str = sprintf('[%.2f - %.2f]', roi_ci(1), roi_ci(2));
        
        ci_t = table({ci_str}, 'VariableNames', {'CI_95'});
        
        new_row = [new_row, ci_t];
    end
    final_table = [final_table; new_row];
end

disp(' ')
disp('Summary of LI Results per ROI:')
disp(final_table)

if cfg_LI.method == 3
    disp('Note: The CI_95 represents the 95% confidence interval of the LI estimate.')
    disp('It indicates the range within which the true LI value would lie in about 95% of repeated bootstrap samples.')
end

end

function report_LI(cfg_LI)
switch cfg_LI.method
    case {1,2}
        % Optional: Display results
        disp('                 ')
        a = table(cfg_LI.RoiLabels'); a.Properties.VariableNames{'Var1'} = 'ROI';
        b = table(cfg_LI.Summ_LI'); b.Properties.VariableNames{'Var1'} = 'LI';
        c = table([cfg_LI.L_count;cfg_LI.R_count]'); c.Properties.VariableNames{'Var1'} = 'Left_vs_right';
        d = [a,b,c];
        disp(d)
    case 3
        
        % Convert each to a column (4x1)
        Summ_LI = cfg_LI.Summ_LI(:);
        RoiLabels = cfg_LI.RoiLabels(:);
        L_vertices_total = cfg_LI.L_count(:);
        R_vertices_total = cfg_LI.R_count(:);
        
        % Combine L and R vertex counts into a single string per row
        Vertices_total = arrayfun(@(l, r) sprintf('%d  %d', l, r), L_vertices_total, R_vertices_total, 'UniformOutput', false);
        
        % Assuming all other arrays are also 1x4, just transpose them:
        CI_strings = cfg_LI.CI_strings(:);
        CI_widths = cfg_LI.CI_widths(:);
        
        % Create a table with the combined vertices column
        T = table(RoiLabels, Summ_LI, CI_strings, CI_widths, Vertices_total, ...
            'VariableNames', {'ROI', 'LI', 'CI_95', 'CI_Width', 'Left_vs_right'});
        disp(T);
end
end

function export_LI_results(cfg_LI)
% EXPORT_LI_RESULTS Save any interval mode as MAT or tab-delimited text.

folderPath = strtrim(char(cfg_LI.savedir));
if isempty(folderPath)
    warning('No Results folder was specified. LI results were not saved.');
    return;
end
if ~exist(folderPath, 'dir')
    [ok, msg] = mkdir(folderPath);
    if ~ok
        warning('Could not create Results folder "%s": %s', folderPath, msg);
        return;
    end
end

baseName = strtrim(char(cfg_LI.sname));
if isempty(baseName)
    baseName = 'LI_results';
end
[~, baseName, ~] = fileparts(baseName);
baseName = regexprep(baseName, '[<>:"/\\|?*]', '_');
baseName = regexprep(baseName, '\s+', '_');

switch cfg_LI.method
    case 1
        methodLabel = 'Magnitude';
    case 2
        methodLabel = 'Count';
    case 3
        methodLabel = 'Bootstrap';
    otherwise
        methodLabel = sprintf('Method%d', cfg_LI.method);
end

switch cfg_LI.Tinterval
    case 1
        intervalLabel = 'Single';
    case 2
        intervalLabel = 'Average';
    case 3
        intervalLabel = 'Windows';
    otherwise
        intervalLabel = sprintf('Interval%d', cfg_LI.Tinterval);
end

fileStem = sprintf('%s_%s_%s', baseName, intervalLabel, methodLabel);

if ~isfield(cfg_LI, 'saveFormat') || isempty(cfg_LI.saveFormat)
    saveFormat = 1;
else
    saveFormat = cfg_LI.saveFormat;
end

switch saveFormat
    case 1
        filename = fullfile(folderPath, [fileStem, '.mat']);
        sLI = build_LI_export_struct(cfg_LI, methodLabel, intervalLabel);
        save(filename, '-struct', 'sLI');
    case 2
        filename = fullfile(folderPath, [fileStem, '.tsv']);
        if ~write_LI_tsv(cfg_LI, filename, methodLabel, intervalLabel)
            return;
        end
    otherwise
        warning('Unknown save format code %g. LI results were not saved.', saveFormat);
        return;
end

fprintf('LI results saved to:\n%s\n', filename);
end


function sLI = build_LI_export_struct(cfg_LI, methodLabel, intervalLabel)
% Create a compact result structure without copying the full source matrix.

sLI = struct();
sLI.format_version = 2;
sLI.method_code = cfg_LI.method;
sLI.method = methodLabel;
sLI.interval_code = cfg_LI.Tinterval;
sLI.interval_method = intervalLabel;
sLI.analysis_range_s = cfg_LI.timerange;
sLI.ROI_labels = cfg_LI.RoiLabels(:);
sLI.source_comment = cfg_LI.Comment;
sLI.threshold_type = cfg_LI.Threshtype;
sLI.threshold_ratio = cfg_LI.Ratio4Threshold;

if cfg_LI.Tinterval == 3 && isfield(cfg_LI, 'final_LI')
    sLI.LI = cfg_LI.final_LI;
    sLI.windows_s = cfg_LI.windows;
    sLI.window_centers_s = mean(cfg_LI.windows, 2);
    if isfield(cfg_LI, 'LI_summary')
        sLI.summary = cfg_LI.LI_summary;
    end
    if isfield(cfg_LI, 'final_CI')
        sLI.CI_95 = cfg_LI.final_CI;
    end
else
    sLI.LI = cfg_LI.Summ_LI(:);
    if isfield(cfg_LI, 'Summ_CI')
        sLI.CI_95 = cfg_LI.Summ_CI;
    end
    if isfield(cfg_LI, 'L_count')
        sLI.left_measure = cfg_LI.L_count(:);
        sLI.right_measure = cfg_LI.R_count(:);
    end
    if isfield(cfg_LI, 'threshold')
        sLI.threshold = cfg_LI.threshold;
    end
end

if cfg_LI.method == 3
    sLI.bootstrap_samples = cfg_LI.n_resampling;
    sLI.threshold_divisions = cfg_LI.divs;
    sLI.resample_ratio = cfg_LI.RESAMPLE_RATIO;
end
end


function success = write_LI_tsv(cfg_LI, filename, methodLabel, intervalLabel)
% Write a genuine tab-delimited text file for any interval mode.

success = false;
fid = fopen(filename, 'wt');
if fid == -1
    warning('Could not open output file for writing: %s', filename);
    return;
end
closeFile = onCleanup(@() fclose(fid));

fprintf(fid, '# LI method\t%s\n', methodLabel);
fprintf(fid, '# Interval method\t%s\n', intervalLabel);
if numel(cfg_LI.timerange) >= 2
    fprintf(fid, '# Analysis range (s)\t%.10g\t%.10g\n', ...
        cfg_LI.timerange(1), cfg_LI.timerange(2));
end

if cfg_LI.Tinterval == 3 && isfield(cfg_LI, 'final_LI')
    [nWindows, nRois] = size(cfg_LI.final_LI);
    hasCI = cfg_LI.method == 3 && isfield(cfg_LI, 'final_CI') && ...
        size(cfg_LI.final_CI,1) == nWindows && ...
        size(cfg_LI.final_CI,2) >= nRois && size(cfg_LI.final_CI,3) >= 2;
    
    if hasCI
        fprintf(fid, 'Window\tStart_s\tEnd_s\tCenter_s\tROI\tLI\tCI_Lower\tCI_Upper\n');
    else
        fprintf(fid, 'Window\tStart_s\tEnd_s\tCenter_s\tROI\tLI\n');
    end
    
    for w = 1:nWindows
        centerTime = mean(cfg_LI.windows(w,:));
        for r = 1:nRois
            roiLabel = regexprep(char(cfg_LI.RoiLabels{r}), '[\t\r\n]', ' ');
            if hasCI
                fprintf(fid, '%d\t%.10g\t%.10g\t%.10g\t%s\t%.10g\t%.10g\t%.10g\n', ...
                    w, cfg_LI.windows(w,1), cfg_LI.windows(w,2), centerTime, ...
                    roiLabel, cfg_LI.final_LI(w,r), ...
                    cfg_LI.final_CI(w,r,1), cfg_LI.final_CI(w,r,2));
            else
                fprintf(fid, '%d\t%.10g\t%.10g\t%.10g\t%s\t%.10g\n', ...
                    w, cfg_LI.windows(w,1), cfg_LI.windows(w,2), centerTime, ...
                    roiLabel, cfg_LI.final_LI(w,r));
            end
        end
    end
else
    roiLabels = cfg_LI.RoiLabels(:);
    if cfg_LI.method == 3 && isfield(cfg_LI, 'Summ_CI')
        fprintf(fid, 'ROI\tLI\tCI_Lower\tCI_Upper\tCI_Width\tLeft_Measure\tRight_Measure\n');
        for r = 1:numel(roiLabels)
            roiLabel = regexprep(char(roiLabels{r}), '[\t\r\n]', ' ');
            fprintf(fid, '%s\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n', ...
                roiLabel, cfg_LI.Summ_LI(r), ...
                cfg_LI.Summ_CI(r,1), cfg_LI.Summ_CI(r,2), ...
                cfg_LI.Summ_CI(r,2) - cfg_LI.Summ_CI(r,1), ...
                cfg_LI.L_count(r), cfg_LI.R_count(r));
        end
    else
        fprintf(fid, 'ROI\tLI\tLeft_Measure\tRight_Measure\n');
        for r = 1:numel(roiLabels)
            roiLabel = regexprep(char(roiLabels{r}), '[\t\r\n]', ' ');
            fprintf(fid, '%s\t%.10g\t%.10g\t%.10g\n', ...
                roiLabel, cfg_LI.Summ_LI(r), cfg_LI.L_count(r), cfg_LI.R_count(r));
        end
        if isfield(cfg_LI, 'threshold')
            fprintf(fid, '# Threshold\t%.10g\n', cfg_LI.threshold);
        end
    end
end

success = true;
clear closeFile
end

function globmax_rois = compute_globmax_rois(cfg_LI)

%- Global max, to avoid LI biases
globmax_rois = [];
for ii = 1:length(cfg_LI.RoiIndices)
    curr_subregion = cfg_LI.sScout.Scouts(cfg_LI.RoiIndices{ii});
    idx_region = [curr_subregion.Vertices];
    % Calculate mean across windows
    mdwin = [];
    for j = 1:size(cfg_LI.windows,1)
        [~,timind1] = min(abs(cfg_LI.Time - cfg_LI.windows(j,1)));
        [~,timind2] = min(abs(cfg_LI.Time - cfg_LI.windows(j,2)));
        dwin = cfg_LI.ImageGridAmp(idx_region,timind1:timind2);
        mdwin(j,:) = mean(dwin,2);
    end
    globmax_rois(ii) = max(mdwin(:));
end
end


function [sScout, ProtocolInfo] = convertDesikanKillianyScout(sResultP)
% Convert the Desikan-Killiany scout to select scouts

ProtocolInfo = bst_get('ProtocolInfo');
SurfaceFile = load(fullfile(ProtocolInfo.SUBJECTS, sResultP.SurfaceFile));

Scouts = [];
sScout = [];
for i = 1:length(SurfaceFile.Atlas)
    if contains(SurfaceFile.Atlas(i).Name, {'Desikan-Killiany'})
        Scouts = SurfaceFile.Atlas(i).Scouts;
    end
end
sScout.Scouts = Scouts;

expectedRegions = {...
    'bankssts L', 'bankssts R', ...
    'caudalanteriorcingulate L', 'caudalanteriorcingulate R', ...
    'caudalmiddlefrontal L', 'caudalmiddlefrontal R', ...
    'cuneus L', 'cuneus R', ...
    'entorhinal L', 'entorhinal R', ...
    'frontalpole L', 'frontalpole R', ...
    'fusiform L', 'fusiform R', ...
    'inferiorparietal L', 'inferiorparietal R', ...
    'inferiortemporal L', 'inferiortemporal R', ...
    'insula L', 'insula R', ...
    'isthmuscingulate L', 'isthmuscingulate R', ...
    'lateraloccipital L', 'lateraloccipital R', ...
    'lateralorbitofrontal L', 'lateralorbitofrontal R', ...
    'lingual L', 'lingual R', ...
    'medialorbitofrontal L', 'medialorbitofrontal R', ...
    'middletemporal L', 'middletemporal R', ...
    'paracentral L', 'paracentral R', ...
    'parahippocampal L', 'parahippocampal R', ...
    'parsopercularis L', 'parsopercularis R', ...
    'parsorbitalis L', 'parsorbitalis R', ...
    'parstriangularis L', 'parstriangularis R', ...
    'pericalcarine L', 'pericalcarine R', ...
    'postcentral L', 'postcentral R', ...
    'posteriorcingulate L', 'posteriorcingulate R', ...
    'precentral L', 'precentral R', ...
    'precuneus L', 'precuneus R', ...
    'rostralanteriorcingulate L', 'rostralanteriorcingulate R', ...
    'rostralmiddlefrontal L', 'rostralmiddlefrontal R', ...
    'superiorfrontal L', 'superiorfrontal R', ...
    'superiorparietal L', 'superiorparietal R', ...
    'superiortemporal L', 'superiortemporal R', ...
    'supramarginal L', 'supramarginal R', ...
    'temporalpole L', 'temporalpole R', ...
    'transversetemporal L', 'transversetemporal R'...
    };

% Handle case when number of anatomical regions are not identical to atlas regions
actualRegions = {sScout.Scouts.Label};
missingRegions = setdiff(expectedRegions, actualRegions);

% Assuming sScout.Scouts is not empty and has at least one scout
if ~isempty(sScout.Scouts)
    % Identify all fields from the first scout as a template
    fieldNames = fieldnames(sScout.Scouts(1));
    % Prepare an empty scout template with all fields
    emptyScout = cell2struct(cell(length(fieldNames), 1), fieldNames, 1);
    
    % Default empty values for known fields
    emptyScout.Label = ''; % Update as necessary
    emptyScout.Vertices = [];
    emptyScout.Seed = 0; % Or any appropriate 'empty' value
    
    % For any other fields in your structure, set an appropriate 'empty' value
    % For example, if there's a 'Color' field, you might do:
    % emptyScout.Color = [0, 0, 0]; % Assuming color is RGB
    % Adjust the above line for each additional field with a sensible 'empty' value
    
    % Now, handle missing regions with this updated emptyScout
    for i = 1:length(missingRegions)
        emptyScout.Label = missingRegions{i};
        % Insert the empty scout at the correct position
        idx = find(strcmp(expectedRegions, missingRegions{i}));
        sScout.Scouts = [sScout.Scouts(1:idx-1), emptyScout, sScout.Scouts(idx:end)];
    end
else
    warning('sScout.Scouts is empty, cannot determine structure fields.');
end
end

function [RoiLabels, RoiIndices] = defineROIs_DK()
% Define regions of interest (ROIs)

AngSmg   = [15,16,63,64];
Front    = [3,4,5,6,11,12,25,26,29,30,33,34,37,38,39,40,41,42,49,50,53,54,55,56,57,58];
LatFront = [5,6,11,12,37,38,39,40,41,42,55,56,57,58];
LatTemp  = [1,2,17,18,31,32,61,62,65,66,67,68];
PeriSyl  = [15,16,37,38,41,42,61,62,63,64];
Tanaka   = [37,38,41,42,61,62,63,64];
Temp     = [1,2,9,10,13,14,17,18,19,20,27,28,31,32,35,36,61,62,65,66,67,68];
Whole    = 1:68;

RoiLabels = {'AngSmg', 'Front','LatFront','LatTemp', 'PeriSyl', 'Tanaka','Temp','Whole'};
RoiIndices = {AngSmg, Front, LatFront, LatTemp, PeriSyl, Tanaka, Temp, Whole};

end