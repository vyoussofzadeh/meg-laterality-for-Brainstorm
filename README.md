# MEG Laterality for Brainstorm

This repository hosts a Brainstorm-compatible pipeline for analyzing brain laterality using magnetoencephalography (MEG) task responses. The pipeline leverages [Brainstorm](https://neuroimage.usc.edu/brainstorm/), a popular MATLAB toolbox, to facilitate robust analysis of hemispheric differences and lateralization indices (LI) in brain activity.

![Tool Interface](figures/fig1_tool.png)

*Figure 1: Example interface of the LI computation process in Brainstorm, showing the options for time intervals, LI methods, thresholding, and output preferences.*

## LI methods

The pipeline provides several methods for computing a Lateralization Index (LI):

- **Source Magnitude Method**: Compares source amplitude across hemispheres.
- **Counting Method**: Counts the number of suprathreshold vertices in left vs. right hemispheres.
- **Bootstrapping Method**: Uses resampling to obtain LI confidence intervals (95% CI), providing estimates of uncertainty around the LI measurement.

## Time intervals

You can also select different time intervals for analysis:
- **Specific Time Interval**: Analyze a predefined time range.
- **Averaged Time Interval**: Compute average activity across a chosen time segment.
- **Window-based Analysis**: Segment data into overlapping or non-overlapping windows, computing LI for each segment to capture temporal evolution of lateralization.

## Atlas (ROIs)

This pipeline currently supports the **Human Connectome Project (HCP) MMP1.0 atlas** and the **Desikan-Killiany (DK) atlas** for defining Regions of Interest (ROIs). Each atlas provides a distinct level of anatomical detail, allowing you to choose the granularity that best fits your research questions.

### HCP Atlas (MNI space only)
- **Description**: Utilizes the HCP MMP1.0 multimodal atlas aligned to a symmetrical MNI template.
- **Resolution**: High-resolution cortical parcellation, enabling a fine-grained assessment of hemispheric dominance.
- **ROIs**: Bilateral ROIs are grouped into functional categories (e.g., Angular, Frontal, Temporal, Lateral).
- **Reference**: See [Glasser et al. (2016), *Nature*](https://www.nature.com/articles/nature18933) for more details.

**Workflow:**
1. **Co-Registration**: Ensure your subject’s anatomy is projected onto the default Brainstorm anatomy in MNI space. To ensure accurate lateralization analysis, it's advisable to first project individual source maps into the default anatomy (e.g., right-click on the source map and select `Default_anatomy > Cortex_15002V`). Then, import the atlas by opening the source file, navigating to the Scout panel, going to the atlas tab, and loading the atlas. Finally, run the LI pipeline on the projected map. This approach helps maintain consistency in spatial alignment across subject anatomies and the atlas.
2. **Load the HCP Atlas**: In the Brainstorm Scout panel, load the HCP MMP1.0 atlas.
3. **Configure Pipeline**: The pipeline functions (`defineROIs_HCP`, `convertHCPScout`) are preconfigured to handle HCP-based ROIs.

### Desikan-Killiany (DK) Atlas (Individual/Native or MNI/Default space)
- **Description**: A widely used, lower-resolution cortical parcellation suitable for broader, more anatomically coarse analyses.
- **Resolution**: Coarser grained, beneficial for initial screenings or simpler ROI definitions.
- **ROIs**: DK-based ROIs are grouped into categories (e.g., AngSmg, Front, LatFront, LatTemp, PeriSyl, Tanaka, Temp, Whole).
- **No Extra Atlas Import Required**: If your anatomy is already segmented using FreeSurfer, the DK atlas labels are typically integrated into the default parcellation available in Brainstorm. This means you do not need to manually import a separate DK atlas—just select the DK scouts within Brainstorm.

**Workflow:**
1. **Co-Registration**: If using individual/native space, ensure that data is appropriately registered. For MNI space, follow the same steps as for HCP.
2. **Default Availability**: DK parcellation is typically included by default after FreeSurfer segmentation. Simply open the Scout panel in Brainstorm to access DK ROIs.
3. **Configure Pipeline**: The pipeline functions (`defineROIs_DK`, `convertDesikanKillianyScout`) can be applied directly to these default DK scouts.

## Thresholding Mechanism

Thresholding plays a crucial role in determining which brain regions are active and contributes to the calculation of LI. The pipeline provides three threshold types:

- **Global-max**: The threshold is set based on the maximum amplitude observed across all time points and all regions.
- **Time-max**: The threshold is based on the maximum amplitude during a specific time of interest across all regions.
- **Region-max**: The threshold is determined by the maximum amplitude within each region of interest during the time of interest.

Users can adjust the threshold ratio to control the sensitivity of region inclusion in the LI calculation. This ratio is applied to the determined maximum value according to the selected threshold type.

## Prerequisites

- **MATLAB**: Tested on MATLAB version 2024.
- **Brainstorm**: Download and install from [the official website](https://neuroimage.usc.edu/brainstorm). Ensure Brainstorm is on your MATLAB path.
- **MEG Data**: Source-level MEG data processed in Brainstorm (containing `ImageGridAmp`).
- **HCP Atlas**: Ensure your subject’s anatomy is co-registered with the HCP MMP1.0 atlas or a compatible symmetrical MNI atlas. A compatible version `scout_mmp_in_mni_symmetrical.mat` is included in this repository. Simply import it into your Brainstorm using the atlas panel. Make sure to import it to the source file that is projected to the default MNI space (see below).
 
## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/your-username/meg-laterality-for-Brainstorm.git

2. **Add to MATLAB Path:**:
   ```bash
   addpath('path_to_meg-laterality-for-Brainstorm');
   savepath;
   
### Step-by-Step in Brainstorm

1. **Open Brainstorm** and load your MEG protocol containing source-level results.

2. **Select a Source File**: in BS, drag and drop the source file to the 'Process1' window

3. **Run the Custom Process**:
    - In the Brainstorm **Process** tab, select **Compute LI (HCP or DK atlas, surface-based)**.
    - Adjust the process parameters as needed:
      - **Time interval method**: Specific, Averaged, or Window-based.
      - **Window parameters**: Window length and overlap (if using Window-based intervals).
      - **LI computation method(s)**: Source magnitude, Counting, and/or Bootstrapping.
      - **Bootstrap parameters**: If bootstrapping is enabled, define divisions, resampling iterations, and resample ratio.
      - **Thresholding**: Choose the threshold type and ratio.
      - **Output**: Specify the output directory and filename.
    - Click **Run** to execute the pipeline.

4. **Results**:
    - The pipeline computes LI for all defined ROIs and time intervals.
    - It generates output files (e.g., `.xls`) containing LI values, vertex counts, and (if bootstrapping is used) 95% confidence intervals.
    - Summaries are displayed, and if window-based analysis is used, LI evolution over time can be plotted.

5. **Check Outputs**:
    - Review the MATLAB command window for logs.
    - Check the specified output directory for result files.
    - If bootstrapping was performed, examine the CI columns in the output files to assess the precision and stability of the LI estimates.
### Example Code Snippet
If you prefer scripting, call the process directly in MATLAB:

```matlab
% Process: Compute LI (HCP or DK atlas, surface-based)
sFiles = bst_process('CallProcess', 'process_computeLI', sFiles, [], ...
    'atlas',                 1, ...  % HCP (MMP1.0)
    'atlas_note',            {}, ...
    'timeParamsLabel',       {}, ...
    'timeParams',            [], ...
    'twindow1',              300, ...
    'toverlap',              50, ...
    'liMethodsLabel',        {}, ...
    'liMethods',             [], ...
    'methodSource',          1, ...
    'methodCounting',        1, ...
    'methodBootstrap',       1, ...
    'bootstrapParamsLabel',  {}, ...
    'bootstrapParams',       [], ...
    'divs',                  10, ...
    'n_resampling',          20, ...
    'RESAMPLE_RATIO',        50, ...
    'timeIntervalLabel',     {}, ...
    'window',                3, ...  % Window based
    'poststim_custom',       [0, 1.95], ...
    'effectLabel',           {}, ...
    'effect',                3, ...  % Absolute
    'thresholdLabel',        {}, ...
    'threshold',             [], ...
    'threshtype1',           3, ...  % Region-max
    'ratio4threshold',       50, ...
    'threshold_for_windows', {}, ...
    'outputLabel',           {}, ...
    'savedir',               'SELECTED PATH', ...
    'sname',                 'NAME OF SAVED FILE', ...
    'saveMat',               0, ...
    'plotResults',           1, ...
    'noteApplicability',     {});

### Interpreting Results
- **LI**: A positive LI value suggests left-hemisphere dominance, while a negative LI suggests right-hemisphere dominance.
- **CI_95 (Bootstrapping)**: The 95% confidence interval provides an estimate of the uncertainty in LI. A narrower CI indicates more stable results.
- **Window-based Analysis**: By examining LI across multiple time windows, you can understand how hemispheric dominance changes over time or in response to different tasks or stimuli.

### Support and Contributions
If you have questions, encounter issues, or would like to contribute improvements:
- **Issues**: Open an issue on this repository for troubleshooting or feature requests.
- **Pull Requests**: Contributions are welcome via pull requests. Please ensure your code is well-documented and tested.
