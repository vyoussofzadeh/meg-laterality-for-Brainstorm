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
addpath('path_to_meg-laterality-for-Brainstorm');
savepath;


