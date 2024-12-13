# MEG Laterality for Brainstorm

This repository hosts a Brainstorm-compatible pipeline for analyzing brain laterality using magnetoencephalography (MEG) task responses. The pipeline leverages [Brainstorm](https://neuroimage.usc.edu/brainstorm/), a popular MATLAB toolbox, to facilitate robust analysis of hemispheric differences and lateralization indices (LI) in brain activity.

## Overview

The pipeline provides several methods for computing a Lateralization Index (LI):

- **Source Magnitude Method**: Compares source amplitude across hemispheres.
- **Counting Method**: Counts the number of suprathreshold vertices in left vs. right hemispheres.
- **Bootstrapping Method**: Uses resampling to obtain LI confidence intervals (95% CI), providing estimates of uncertainty around the LI measurement.

You can also select different time intervals for analysis:
- **Specific Time Interval**: Analyze a predefined time range.
- **Averaged Time Interval**: Compute average activity across a chosen time segment.
- **Window-based Analysis**: Segment data into overlapping or non-overlapping windows, computing LI for each segment to capture temporal evolution of lateralization.

## ROIs and Atlas

This pipeline uses the **Human Connectome Project (HCP) atlas** to define Regions of Interest (ROIs). The HCP MMP1.0 atlas is a high-resolution cortical parcellation providing detailed anatomical and functional organization of the cerebral cortex. Key points:

- **HCP Atlas**: A multimodal atlas aligned to symmetrical MNI space, ensuring fair hemisphere-to-hemisphere comparisons.
- **ROIs**: The code groups bilateral ROIs into categories (e.g., Angular, Frontal, Temporal, Lateral). Users can focus on specific ROI sets relevant to their research questions.
- **Missing or Incomplete Regions**: The pipeline checks for expected HCP-based ROIs and creates empty placeholders for missing regions, maintaining consistent indexing.

By leveraging the HCP atlas, the pipeline offers a comprehensive and fine-grained approach to assessing hemispheric dominance across various cortical regions.

For more information on the HCP atlas, see [Glasser et al. (2016), *Nature*](https://www.nature.com/articles/nature18933).

## Prerequisites

- **MATLAB**: Tested on recent MATLAB versions.
- **Brainstorm**: Download and install from [the official website](https://neuroimage.usc.edu/brainstorm). Ensure Brainstorm is on your MATLAB path.
- **MEG Data**: Source-level MEG data processed in Brainstorm (containing `ImageGridAmp`).
- **HCP Atlas**: Ensure your subject’s anatomy is co-registered with the HCP MMP1.0 atlas or a compatible symmetrical MNI atlas.

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

2. **Select a Results File**: Choose the `...results...mat` file you wish to analyze.

3. **Run the Custom Process**:
    - In the Brainstorm **Process** tab, select **Compute LI, surface-based, HCP atlas**.
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

### Interpreting Results
- **LI**: A positive LI value suggests left-hemisphere dominance, while a negative LI suggests right-hemisphere dominance.
- **CI_95 (Bootstrapping)**: The 95% confidence interval provides an estimate of the uncertainty in LI. A narrower CI indicates more stable results.
- **Window-based Analysis**: By examining LI across multiple time windows, you can understand how hemispheric dominance changes over time or in response to different tasks or stimuli.

### Support and Contributions
If you have questions, encounter issues, or would like to contribute improvements:
- **Issues**: Open an issue on this repository for troubleshooting or feature requests.
- **Pull Requests**: Contributions are welcome via pull requests. Please ensure your code is well-documented and tested.

### Example Code Snippet
If you prefer scripting, call the process directly in MATLAB:

```matlab
   % sFiles: cell array of Brainstorm result filenames
   sProcess = bst_process('CallProcess', 'process_computeLI_HCP', sFiles, [], ...
    'twindow',     300, ...
    'toverlap',    50, ...
    'methodSource', 0, ...
    'methodCounting', 1, ...
    'methodBootstrap', 1, ...
    'divs',        10, ...
    'n_resampling', 200, ...
    'RESAMPLE_RATIO', 75, ...
    'window',       3, ...
    'poststim',     [0.1, 0.5], ...
    'effect',       1, ...
    'threshtype',   1, ...
    'ratio4threshold', 20, ...
    'savedir',      '/path/to/save/dir', ...
    'sname',        'analysis_results');
