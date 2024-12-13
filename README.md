# MEG Laterality for Brainstorm

This repository hosts a Brainstorm-compatible pipeline designed for analyzing brain laterality using magnetoencephalography (MEG) task responses. The pipeline leverages [Brainstorm](https://neuroimage.usc.edu/brainstorm/), a popular MATLAB toolbox, to facilitate robust analysis of hemispheric differences and lateralization indices (LI) in brain activity.

## Overview

The pipeline provides several methods for computing a Lateralization Index (LI), including:
- **Source Magnitude Method**: Directly compares source amplitude across hemispheres.
- **Counting Method**: Counts the number of suprathreshold vertices in left vs. right hemispheres.
- **Bootstrapping Method**: Uses resampling to obtain LI confidence intervals (95% CI), providing estimates of uncertainty around the LI measurement.

Additionally, the pipeline allows for flexible time interval selection:
- **Specific Time Interval**: Analyze within a predefined time range.
- **Averaged Time Interval**: Compute average activity across a chosen time segment.
- **Window-based Analysis**: Segment the data into overlapping or non-overlapping time windows and compute LI for each.

This flexibility helps users characterize how laterality indices evolve over time or differ under various conditions.

## Prerequisites

- **MATLAB**: The pipeline is developed and tested on recent MATLAB versions.
- **Brainstorm**: Download and install Brainstorm from [the official website](https://neuroimage.usc.edu/brainstorm). Ensure you have a working Brainstorm environment and that Brainstorm is on your MATLAB path.
- **MEG Data**: You should have processed MEG source-level data in Brainstorm results format. The pipeline expects Brainstorm "results" files that contain `ImageGridAmp` and associated time vectors.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/your-username/meg-laterality-for-Brainstorm.git
