function varargout = process_computeLI_HCP(varargin )

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
% Author: Vahab YoussofZadeh, 2024
% Last Update: 12/13/24
% Changes: Window-based LI analysis added.

eval(macro_method);

end

function sProcess = GetDescription() %#ok<DEFNU>
% PROCESS_COMPUTELI_HCP: Compute the Lateralization Index (LI) on HCP atlas-based source data.
% Provides several methods (source magnitude, counting, bootstrapping) and time interval strategies (specific, averaged, window-based).
% Supports thresholding and optional saving of results as .mat and plotting.

sProcess.Comment     = 'Compute LI (HCP atlas, surface-based)';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 337;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/CoregisterSubjects';
sProcess.InputTypes  = {'results'};
sProcess.OutputTypes = {'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

%% Time Window Parameters
% Heading as a label
sProcess.options.timeParamsLabel.Comment = '<B>Time Window Parameters:</B>';
sProcess.options.timeParamsLabel.Type    = 'label';
sProcess.options.timeParamsLabel.Value   = {};

% Group to hold related options (no special formatting here)
sProcess.options.timeParams.Comment = '';
sProcess.options.timeParams.Type    = 'group';
sProcess.options.timeParams.Value   = [];

sProcess.options.twindow.Comment    = 'Window length (ms): For window-based methods, defines the length of each time window.';
sProcess.options.twindow.Type       = 'value';
sProcess.options.twindow.Value      = {300, 'ms', 100, 1000, 1};

sProcess.options.toverlap.Comment   = 'Overlap between windows (%): For window-based analysis, the percentage overlap between consecutive windows.';
sProcess.options.toverlap.Type      = 'value';
sProcess.options.toverlap.Value     = {50, '%', 0, 100, 1};

%% LI Computation Methods
% Heading as a label
sProcess.options.liMethodsLabel.Comment = '<B>Lateralization Index Methods:</B>';
sProcess.options.liMethodsLabel.Type    = 'label';
sProcess.options.liMethodsLabel.Value   = {};

% Group
sProcess.options.liMethods.Comment  = '';
sProcess.options.liMethods.Type     = 'group';
sProcess.options.liMethods.Value    = [];

sProcess.options.methodSource.Comment = 'Source Magnitude Method: Compute LI based on source amplitude values.';
sProcess.options.methodSource.Type    = 'checkbox';
sProcess.options.methodSource.Value   = 0;

sProcess.options.methodCounting.Comment = 'Counting Method: Compute LI based on counting suprathreshold vertices.';
sProcess.options.methodCounting.Type    = 'checkbox';
sProcess.options.methodCounting.Value   = 0;

sProcess.options.methodBootstrap.Comment = 'Bootstrapping Method: Compute LI with confidence intervals via resampling.';
sProcess.options.methodBootstrap.Type    = 'checkbox';
sProcess.options.methodBootstrap.Value   = 0;

%% Bootstrap Parameters
% Heading as a label
sProcess.options.bootstrapParamsLabel.Comment = '<B>Bootstrap Parameters:</B>';
sProcess.options.bootstrapParamsLabel.Type    = 'label';
sProcess.options.bootstrapParamsLabel.Value   = {};

% Group
sProcess.options.bootstrapParams.Comment = '';
sProcess.options.bootstrapParams.Type    = 'group';
sProcess.options.bootstrapParams.Value   = [];

sProcess.options.divs.Comment = 'Number of divisions: For bootstrapping LI, how many threshold divisions to consider.';
sProcess.options.divs.Type    = 'value';
sProcess.options.divs.Value   = {10, '', 1, 100, 1};

sProcess.options.n_resampling.Comment = 'Number of resampling iterations: More iterations yield more stable CI at the cost of computation time.';
sProcess.options.n_resampling.Type    = 'value';
sProcess.options.n_resampling.Value   = {200, '', 1, 1000, 1};

sProcess.options.RESAMPLE_RATIO.Comment = 'Resample ratio (%): Percentage of data to resample in each bootstrap iteration.';
sProcess.options.RESAMPLE_RATIO.Type    = 'value';
sProcess.options.RESAMPLE_RATIO.Value   = {75, '%', 0, 100, 1, 1};

%% Time Interval Selection
% Label for interval analysis
sProcess.options.timeIntervalLabel.Comment = '<B>Time Interval Analysis:</B> Choose how time is handled.';
sProcess.options.timeIntervalLabel.Type    = 'label';
sProcess.options.timeIntervalLabel.Value   = {};

sProcess.options.window.Comment = '';
sProcess.options.window.Type = 'combobox';
sProcess.options.window.Value = {1, {'Specific Time Interval', 'Averaged Time Interval', 'Window based'}};

% Specific time interval setting
sProcess.options.poststim_custom.Comment = 'Enter specific time interval: Only used if "Specific Time Interval" is selected.';
sProcess.options.poststim_custom.Type    = 'poststim';
sProcess.options.poststim_custom.Value   = [];

%% Effect Type
% Label for effect type
sProcess.options.effectLabel.Comment = '<B>Effect Type:</B> Choose how to process the data (positive, negative, absolute values).';
sProcess.options.effectLabel.Type    = 'label';
sProcess.options.effectLabel.Value   = {};

sProcess.options.effect.Comment = '';
sProcess.options.effect.Type = 'combobox';
sProcess.options.effect.Value = {1, {'Positive values', 'Negative values', 'Absolute values'}};

%% Threshold Settings
% Label for threshold settings
sProcess.options.thresholdLabel.Comment = '<B>Threshold Settings:</B>';
sProcess.options.thresholdLabel.Type    = 'label';
sProcess.options.thresholdLabel.Value   = {};

sProcess.options.threshold.Comment = '';
sProcess.options.threshold.Type    = 'group';
sProcess.options.threshold.Value   = [];

sProcess.options.threshtype.Comment = ['Threshold type: Select how to determine the max amplitude for thresholding. ' ...
                                       'In window-based approaches, global maxima (compute_globmax_rois) are used.'];
sProcess.options.threshtype.Type    = 'combobox';
sProcess.options.threshtype.Value   = {1, {'Global-max', 'Time-max', 'Region-max'}};

sProcess.options.threshold_for_windows.Comment = ['<I>Note:</I> For the window-based approach, thresholds rely on global maxima ' ...
                                                  '(from compute_globmax_rois) to ensure consistency across windows.'];
sProcess.options.threshold_for_windows.Type = 'label';
sProcess.options.threshold_for_windows.Value = {};

sProcess.options.ratio4threshold.Comment = 'Threshold ratio (%): Percentage of the determined max value used as threshold.';
sProcess.options.ratio4threshold.Type = 'value';
sProcess.options.ratio4threshold.Value = {20, '%', 0, 100, 1, 1};

%% Output Settings
% Label for output settings
sProcess.options.outputLabel.Comment = '<B>Output Settings:</B>';
sProcess.options.outputLabel.Type    = 'label';
sProcess.options.outputLabel.Value   = {};

sProcess.options.savedir.Comment = 'Saving Dir.: Directory where results will be saved.';
sProcess.options.savedir.Type    = 'text';
sProcess.options.savedir.Value   = '';

sProcess.options.sname.Comment = 'Saving filename: Base name for output files.';
sProcess.options.sname.Type    = 'text';
sProcess.options.sname.Value   = '';

% Additional options
sProcess.options.saveMat.Comment = 'Save results as .mat file (window-based only).';
sProcess.options.saveMat.Type    = 'checkbox';
sProcess.options.saveMat.Value   = 0;

sProcess.options.plotResults.Comment = 'Enable Plotting (window-based only). If selected, temporal evolution of LI is plotted.';
sProcess.options.plotResults.Type    = 'checkbox';
sProcess.options.plotResults.Value   = 1;

sProcess.options.noteApplicability.Comment = ['<I>Note:</I> The .mat saving and plotting options apply only when ' ...
                                              'the "Window based" method is selected.'];
sProcess.options.noteApplicability.Type = 'label';
sProcess.options.noteApplicability.Value = {};

end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput)

% Stop any Brainstorm progress bar that may be running
bst_progress('stop');

% Extract user-defined options
plotResults = sProcess.options.plotResults.Value;
saveMat     = sProcess.options.saveMat.Value;
savedir     = sProcess.options.savedir.Value;

OutputFiles = {};

% ===== LOAD INPUT DATA =====
sResultP = in_bst_results(sInput.FileName, 1);
if isempty(sResultP)
    error('Could not load the selected results file.');
end

% ===== TIME INTERVAL SETTINGS =====
Tinterval = selectTimeInterval(sProcess.options.window.Value{1});

% Determine the time range based on user selection
if Tinterval == 1 || Tinterval == 3
    timerange = sProcess.options.poststim_custom.Value{1};
else
    timerange = 1; % For averaged or other intervals, default/legacy behavior
end

% ===== EFFECT TYPE =====
effect = selectEffectType(sProcess.options.effect.Value{1});

% ===== THRESHOLD SETTINGS =====
Ratio4Threshold = sProcess.options.ratio4threshold.Value{1}/100;
Threshtype = selectThresholdType(sProcess.options.threshtype.Value{1});

% ===== PROCESS IMAGE DATA =====
ImageGridAmp = processImageGridAmp(sResultP.ImageGridAmp, effect);

% Determine max values for various time windows
[AllMax, GlobalMax, t1, t2] = determineMaxValues(Tinterval, ImageGridAmp, sResultP, timerange);

[sScout, ~] = convertHCPScout(sResultP);

% If "Averaged Time Interval" is selected (value=2), and enough data points are available
if sProcess.options.window.Value{1} == 2 && length(sResultP.Time) > 10
    timerange = sProcess.options.poststim.Value{1};
    [~, t1] = min(abs(sResultP.Time - timerange(1)));
    [~, t2] = min(abs(sResultP.Time - timerange(2)));
    ImageGridAmp = mean(ImageGridAmp(:, t1:t2), 2);
end

% If "Window based" method is selected (value=3), set up time windows
if sProcess.options.window.Value{1} == 3
    winLength = sProcess.options.twindow.Value{1};
    overlap   = sProcess.options.toverlap.Value{1};
    
    cfg = [];
    cfg.strt    = timerange(1);
    cfg.spt     = timerange(2);
    cfg.overlap = overlap/1000;
    cfg.linterval = winLength/1000;
    wi = generateTimeWindows(cfg);
else
    wi = [];
end

disp('Selected intervals:')
if sProcess.options.window.Value{1} == 3
    disp(['Interval: [', num2str(timerange(1)), ', ', num2str(timerange(2)), '] sec, Window length: ', ...
          num2str(winLength), ' ms, Overlap: ', num2str(overlap), '%']);
    disp(wi);
else
    disp(['Interval: [', num2str(timerange(1)), ', ', num2str(timerange(2)), '] sec']);
end

% ===== DEFINE ROIs =====
[RoiLabels, RoiIndices] = defineROIs_HCP(sScout);

% ===== CONFIGURE LI COMPUTATION =====
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
cfg_LI.saveMat         = saveMat;
cfg_LI.savedir         = savedir;
cfg_LI.Comment         = sResultP.Comment;
cfg_LI.windows         = wi; % Only meaningful if window-based
cfg_LI.sname           = sProcess.options.sname.Value;
cfg_LI.Time            = sResultP.Time;

% ===== COMPUTE LI BASED ON SELECTED METHODS =====
% At least one method should be selected; if not, warn the user
methodsSelected = [sProcess.options.methodSource.Value, sProcess.options.methodCounting.Value, sProcess.options.methodBootstrap.Value];
if ~any(methodsSelected)
    warning('No LI computation method selected. Please choose at least one method (Source, Counting, or Bootstrapping).');
    return;
end

if sProcess.options.methodSource.Value == 1
    disp('Computing LI using Source Magnitude Method...')
    cfg_LI.method = 1;
    computeLI(cfg_LI); 
    pause(0.2);
end

if sProcess.options.methodCounting.Value == 1
    disp('Computing LI using Counting Method...')
    cfg_LI.method = 2;
    computeLI(cfg_LI); 
    pause(0.2);
end

if sProcess.options.methodBootstrap.Value == 1
    disp('Computing LI with Bootstrapping...')
    cfg_LI.method = 3;
    cfg_LI.divs           = sProcess.options.divs.Value{1};
    cfg_LI.n_resampling   = sProcess.options.n_resampling.Value{1};
    cfg_LI.RESAMPLE_RATIO = sProcess.options.RESAMPLE_RATIO.Value{1} / 100;
    computeLI(cfg_LI);  
    pause(0.2);
end

disp(['LI assessed for interval: ', num2str(timerange(1)), '-', num2str(timerange(2)), ' sec using the HCP-MMP1 atlas.']);
disp('If needed, edit process_computeLI.m after ensuring Brainstorm is running.')
disp('LI analysis is completed!')

OutputFiles = {}; % Adjust if your code generates output files

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
    if contains(SurfaceFile.Atlas(i).Name, {'mmp_in_mni_symmetrical_1.nii_05'})
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
    
    cfg_LI.Summ_LI = Summ_LI;
    cfg_LI.L_count = L_count;
    cfg_LI.R_count = R_count;
    cfg_LI.threshold = threshold;
    cfg_LI.LI_label_out = LI_label_out;
    
    export_LI(cfg_LI)
    report_LI(cfg_LI)
    
elseif isempty(windows) && cfg_LI.method == 3
    
    cfg_LI.report = 0;
    [Summ_LI, Summ_CI, LI_label_out, L_count, R_count, CI_strings, CI_widths]  = computeLI_bootstrap(cfg_LI);
    
    cfg_LI.Summ_LI = Summ_LI;
    cfg_LI.Summ_CI = Summ_CI;
    cfg_LI.L_count = L_count;
    cfg_LI.R_count = R_count;
    cfg_LI.R_count = R_count;
    cfg_LI.CI_strings = CI_strings;
    cfg_LI.CI_widths = CI_widths;
    cfg_LI.LI_label_out = LI_label_out;
    
    export_LI(cfg_LI)
    report_LI(cfg_LI)
    
elseif ~isempty(windows) && cfg_LI.Tinterval == 3
    
    globmax_rois = compute_globmax_rois(cfg_LI);
    
    final_LI = [];
    final_CI = [];
    for j = 1:length(windows)
        [~, cfg_LI.t1] = min(abs(cfg_LI.Time - windows(j,1)));
        [~, cfg_LI.t2] = min(abs(cfg_LI.Time - windows(j,2)));
        cfg_LI.globmax_rois = globmax_rois;
        if cfg_LI.method == 1 || cfg_LI.method == 2
            [Summ_LI, ~, ~, ~, ~] = computeLI_Svalue_Count(cfg_LI);
        elseif cfg_LI.method == 3
            cfg_LI.report = 0;
            disp(['Interval:', num2str(j), '/', num2str(length(windows)), ', [', num2str(windows(j,1)), ',', num2str(windows(j,2)), '] sec'])
            [Summ_LI, Summ_CI, ~] = computeLI_bootstrap(cfg_LI);
            final_CI(j,:,:) = Summ_CI;
        end
        final_LI = [final_LI;Summ_LI];
    end
    
    cfg_LI.final_LI = final_LI;
    if cfg_LI.method == 3, cfg_LI.final_CI = final_CI; end
    
    if cfg_LI.plotResults == 1
        plot_LI(cfg_LI);
    end
    LI_summary = report_tLI(cfg_LI);
    
    cfg_LI.LI_summary = LI_summary;
    
    if cfg_LI.saveMat == 1
        export_tLI(cfg_LI);
    end
    
end
end

function [weighted_li, num_threshvals, L_vertices_above_thresh, R_vertices_above_thresh, CI] = do_LI_bootstrap(cfg_main)
divs = cfg_main.divs;
n_resampling = cfg_main.n_resampling;
RESAMPLE_RATIO = cfg_main.RESAMPLE_RATIO;
RoiIndices = cfg_main.RoiIndices;
MIN_NUM_THRESH_VOXELS = divs / RESAMPLE_RATIO; % Adjust as needed

% Extract the data in the time interval
if cfg_main.Tinterval == 2 || size(cfg_main.ImageGridAmp,2) == 1
    ImageGridAmp = mean(cfg_main.ImageGridAmp,2);
else
    ImageGridAmp = cfg_main.ImageGridAmp(:, cfg_main.t1:cfg_main.t2);
end

% Optional downsampling if too large
if size(ImageGridAmp,2) > 500
    ImageGridAmp = ImageGridAmp(:,1:10:end);
end

LHscout = horzcat(cfg_main.atlas.Scouts(RoiIndices(1:2:end)).Vertices);
RHscout = horzcat(cfg_main.atlas.Scouts(RoiIndices(2:2:end)).Vertices);

LHvals = ImageGridAmp(LHscout, :);
RHvals = ImageGridAmp(RHscout, :);

ROIMax = max([LHvals(:); RHvals(:)]);
threshvals = linspace(0, ROIMax, divs);

weighted_li = 0;
num_threshvals = 0;

cumulative_L_vertices_above_thresh = 0;
cumulative_R_vertices_above_thresh = 0;

% Initialize a matrix to store all bootstrap LI distributions across thresholds
All_TB_LIs_weighted = [];

for thresh_idx = 1:numel(threshvals)
    thresh = threshvals(thresh_idx);
    l_above_thresh = LHvals >= thresh;
    r_above_thresh = RHvals >= thresh;
    
    if sum(l_above_thresh(:)) < MIN_NUM_THRESH_VOXELS && sum(r_above_thresh(:)) < MIN_NUM_THRESH_VOXELS
        break; % Not enough voxels above threshold
    end
    
    cumulative_L_vertices_above_thresh = cumulative_L_vertices_above_thresh + sum(l_above_thresh(:));
    cumulative_R_vertices_above_thresh = cumulative_R_vertices_above_thresh + sum(r_above_thresh(:));
    
    % Compute bootstrap distribution for this threshold
    TB_LIs = bootstrapLI(l_above_thresh, r_above_thresh, n_resampling, RESAMPLE_RATIO);
    
    % Accumulate weighted LI
    weighted_li = weighted_li + mean(TB_LIs) * thresh_idx;
    num_threshvals = thresh_idx;
    
    % Store the weighted TB_LIs for CI calculation
    % We store TB_LIs * thresh_idx so that later we can sum them all and then divide by sum(1:num_threshvals).
    All_TB_LIs_weighted = [All_TB_LIs_weighted, TB_LIs * thresh_idx];
end

if num_threshvals > 0
    weighted_li = (weighted_li / sum(1:num_threshvals)) * 100;
end

L_vertices_above_thresh = round(cumulative_L_vertices_above_thresh / max(num_threshvals,1));
R_vertices_above_thresh = round(cumulative_R_vertices_above_thresh / max(num_threshvals,1));

% Compute the combined bootstrap distribution for the final weighted LI
if num_threshvals > 0
    CombinedLI = (sum(All_TB_LIs_weighted, 2) / sum(1:num_threshvals)) * 100;
    % Calculate 95% CI from CombinedLI
    CI = prctile(CombinedLI, [2.5 97.5]);
else
    CI = [NaN NaN]; % No CI if no thresholds were processed
end
end

function TB_LIs = bootstrapLI(Lvals, Rvals, n_samples, resample_ratio)
% Bootstrap Laterality Index computation for given left and right hemisphere values.
l_n = max(round(numel(Lvals) * resample_ratio), 1);
r_n = max(round(numel(Rvals) * resample_ratio), 1);

TB_LIs = zeros(n_samples, 1);
for i = 1:n_samples
    l_sample = datasample(Lvals, l_n);
    r_sample = datasample(Rvals, r_n);
    TB_LIs(i) = (mean(l_sample) - mean(r_sample)) / (mean(l_sample) + mean(r_sample));
end
end

function wi  = generateTimeWindows (cfg_main)

strt = cfg_main.strt; % strt = -0.5 sec.
spt = cfg_main.spt; % spt = 0 sec.
overlap = cfg_main.overlap; % overlap = 0;
linterval = cfg_main.linterval; % linterval = 0.25 sec.

wi = [];
w1 = strt;
l = linterval;
ov = overlap;
j = 1;

while w1 + l <= spt
    wi(j,:) = [w1, w1 + l];
    j = j + 1;
    w1 = w1 + (ov > 0) * ov + (ov <= 0) * l; % Adjust w1 based on overlap presence
end

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
    Ltemp_region = [curr_subregion(1:2:end).Vertices];
    Rtemp_region = [curr_subregion(2:2:end).Vertices];
    
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
% PLOT_LI Visualizes the Lateralization Index (LI) across multiple time windows.
% This function is intended for window-based analyses, where LI values are computed
% for each defined time window. It plots LI values per window for all ROIs, allowing
% users to observe temporal dynamics in laterality.

% Ensure that windows and final_LI are available
if ~isfield(cfg_LI, 'windows') || isempty(cfg_LI.windows)
    warning('No windows defined in cfg_LI. This plot is primarily intended for the window-based approach.');
    return;
end

if ~isfield(cfg_LI, 'final_LI') || isempty(cfg_LI.final_LI)
    warning('No LI data available to plot (cfg_LI.final_LI is empty).');
    return;
end

% Extract window start times and compute mean of start times for labeling
val = mean(cfg_LI.windows(:,1),2);
val_rounded = round(val,2);

% Check if final_LI and windows have compatible dimensions
[W, R] = size(cfg_LI.final_LI);
if W ~= size(cfg_LI.windows,1)
    warning('Number of windows in cfg_LI.final_LI does not match cfg_LI.windows. Cannot plot LI per window.');
    return;
end

figure('Color', 'w'); % White background figure
plot(cfg_LI.final_LI, 'LineWidth', 1.5);
ylabel('Lateralization Index (LI)');

% Set x-axis to show window start times. Label every other window for readability.
set(gca, 'XTick', 1:2:W, 'XTickLabel', val_rounded(1:2:end));
set(gca, 'FontSize', 8, 'XTickLabelRotation', 90);

% Adjust figure and axes appearance
set(gcf, 'Position', [400, 400, 700, 300]);
xlim([1 length(val_rounded)])
box off
set(gca, 'color', 'none'); % Transparent background

xlabel('Mean Temporal Windows (sec)');

% Set the title based on the method used
switch cfg_LI.method
    case 1
        title('LI by Source Magnitude Method');
    case 2
        title('LI by Vertex Counting Method');
    case 3
        title('LI by Bootstrapping Method');
    otherwise
        title('Lateralization Index');
end

% If ROI labels are available, add a legend
if isfield(cfg_LI, 'RoiLabels') && ~isempty(cfg_LI.RoiLabels)
    legend(cfg_LI.RoiLabels', 'Location', 'best');
else
    legend('off');
    disp('No ROI labels found to display in legend.');
end

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

function export_tLI(cfg_LI)

folderPath = fullfile(cfg_LI.savedir);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
    disp('Folder created successfully.');
else
    disp('Folder already exists.');
end

sname = cfg_LI.sname;

sLI = [];
sLI.final_LI      = cfg_LI.final_LI;
sLI.windows       = cfg_LI.windows;
sLI.RoiLabels     = cfg_LI.RoiLabels;
sLI.LI_summary    = cfg_LI.LI_summary;
sLI.filename      = cfg_LI.Comment;
sLI.Threshtype    = cfg_LI.Threshtype;
sLI.globmax_rois  = cfg_LI.globmax_rois;

switch cfg_LI.method
    case 1
        mtd_label = 'S';
    case 2
        mtd_label = 'C';
    case 3
        mtd_label = 'B';
        sLI.ConfInt     = cfg_LI.final_CI;
        sLI.n_resampling = cfg_LI.n_resampling;
        sLI.divs         = cfg_LI.divs;
end

save(fullfile(folderPath,[sname, '_',mtd_label, '.mat']), '-struct', 'sLI');

end

function export_LI(cfg_LI)
% EXPORT_LI saves the computed LI results to a file, adapting format based on the chosen method.
% For non-bootstrapping methods (1 or 2), it saves a simple table with LI and L/R counts.
% For the bootstrapping method (3), it includes CI and CI width as well.
%
% The function also handles directory creation and clearly reports where the file is saved.

% Ensure savedir is defined
folderPath = fullfile(cfg_LI.savedir);
if isempty(folderPath)
    warning('No saving directory specified (cfg_LI.savedir is empty). Results will not be saved.');
    return;
end

% Create the folder path if it does not exist
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
    disp(['Folder created successfully at: ', folderPath]);
else
    disp(['Using existing directory: ', folderPath]);
end

% Determine method label
switch cfg_LI.method
    case 1
        mtd_label = 'S'; % Source Magnitude
    case 2
        mtd_label = 'C'; % Counting
    case 3
        mtd_label = 'B'; % Bootstrapping
    otherwise
        warning('Unknown LI method code. Cannot determine filename label.');
        mtd_label = 'X';
end

% Construct filename based on Tinterval and method
if cfg_LI.Tinterval == 2
    % Averaged Time Interval case
    if cfg_LI.method == 3
        % Bootstrapping method for averaged interval
        filename = fullfile(folderPath, ['/LI_', mtd_label, '_', cfg_LI.sname, '.xls']);
    else
        % Source or Counting method for averaged interval
        filename = fullfile(folderPath, ['/LI_', mtd_label, '_', cfg_LI.sname, '_Th', num2str(cfg_LI.Ratio4Threshold), '.xls']);
    end
else
    % Specific Time Interval (1) or Window based (3)
    if cfg_LI.method == 3
        % Bootstrapping method for window-based or specific interval
        filename = fullfile(folderPath, ['/LI', mtd_label, '_', cfg_LI.sname, ' T ', num2str(cfg_LI.timerange(1)), '-', num2str(cfg_LI.timerange(2)), '.xls']);
    else
        % Source or Counting method for specific or window intervals
        filename = fullfile(folderPath, ['/LI', mtd_label, '_', cfg_LI.sname, ' T ', num2str(cfg_LI.timerange(1)), '-', num2str(cfg_LI.timerange(2)), '_Thre ', num2str(cfg_LI.Ratio4Threshold), '.xls']);
    end
end

% Open file for writing
tempfile = fopen(filename, 'w');
if tempfile == -1
    warning(['Could not open file for writing: ', filename]);
    return;
end

% Print results depending on the method
switch cfg_LI.method
    case {1, 2}
        % Simple header for Source or Counting method
        fprintf(tempfile, 'ROI\tLI\tLeft_count\tRight_count\n');
        
        L_count1 = cfg_LI.L_count(:);
        R_count1 = cfg_LI.R_count(:);
        
        for i = 1:length(cfg_LI.LI_label_out)
            fprintf(tempfile, '%s\t%f\t%d\t%d\n', ...
                cfg_LI.LI_label_out{i}, cfg_LI.Summ_LI(i), L_count1(i), R_count1(i));
        end
        
        % Print threshold used
        fprintf(tempfile, 'Threshold\t%f\n', cfg_LI.threshold);
        
    case 3
        % Bootstrapping method: Include CI and CI width
        % Print headers including CI and vertex counts
        fprintf(tempfile, 'ROI\tLI\tCI\tCI_Width\tL_Vertices\tR_Vertices\n');
        
        L_count1 = cfg_LI.L_count(:);
        R_count1 = cfg_LI.R_count(:);
        for i = 1:length(cfg_LI.LI_label_out)
            fprintf(tempfile, '%s\t%f\t%s\t%f\t%d\t%d\n', ...
                cfg_LI.LI_label_out{i}, cfg_LI.Summ_LI(i), cfg_LI.CI_strings{i}, cfg_LI.CI_widths(i), L_count1(i), R_count1(i));
        end
end

% Add a newline for separation
fprintf(tempfile, '\n');
fclose(tempfile);

% Display where the results are saved
disp('Results saved to:');
disp(filename);

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

