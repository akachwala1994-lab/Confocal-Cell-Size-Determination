# Confocal-Cell-Size-Determination
# Chlamydomonas Cell Size Analysis Using MATLAB

A MATLAB-based image analysis tool for measuring the cell size of *Chlamydomonas reinhardtii* from Leica Image File (LIF) microscopy datasets. The program combines manual cell boundary identification, binary-mask segmentation, and morphological analysis to determine cell dimensions across multiple imaging series and Z-stacks.

## Overview

This program was developed to quantify cell size using grayscale images, with chlorophyll autofluorescence serving as a visual reference for identifying individual cells. The user manually defines the approximate cell center and outer boundary, after which the program constructs a circular region of interest (ROI). The ROI can be adjusted into an ellipse or replaced with a manually traced boundary to account for variations in cell morphology.

For each imaging series, the program identifies the accepted Z-plane with the largest projected cell area and extracts the corresponding cell dimensions.

### Features

- Direct import of Leica `.lif` files using Bio-Formats.
- Selection of multiple imaging series and Z-stacks.
- Visualization of chlorophyll fluorescence over grayscale images.
- Manual selection of cell center and outer edge.
- Interactive circular, elliptical, and freehand ROI adjustment.
- Calculation of major-axis length, minor-axis length, equivalent diameter, cell area, aspect ratio, and orientation.
- Automatic selection of the largest-area accepted Z-plane per series.
- Export of binary masks, annotated microscopy images, CSV measurements, and MATLAB figures.

## Requirements

- MATLAB
- MATLAB Image Processing Toolbox
- Bio-Formats MATLAB Toolbox

Bio-Formats can be obtained from the [OME Bio-Formats website](https://www.openmicroscopy.org/bio-formats/).

Add the Bio-Formats MATLAB directory to the MATLAB search path before running the script.

## Image Channel Configuration

The program assumes the following channel arrangement:

| Channel | Image type | Purpose |
|---|---|---|
| 1 | GFP fluorescence | Available in the LIF file but not used for cell-size measurement |
| 2 | Grayscale | Identification of the outer cell boundary |
| 3 | Chlorophyll autofluorescence | Cell localization and Z-stack selection |

## Analysis Workflow

### 1. LIF File Import and Z-Stack Selection

The program imports a selected LIF file and identifies its available imaging series, channels, time points, and Z-planes.

For each selected series, the chlorophyll images are displayed as a Z-stack preview. The user selects the optical sections to be analyzed.

The physical pixel dimensions in the X and Y directions are retrieved from the LIF metadata and converted to µm/pixel. Manual calibration values can be supplied if the metadata is unavailable.

### 2. Image Normalization and Channel Overlay

The grayscale and chlorophyll images are independently normalized using the 1st and 99.5th percentiles of their intensity distributions.

Chlorophyll fluorescence is represented in magenta and superimposed on the grayscale image. This overlay provides a visual reference for locating the cell and identifying its outer boundary.

The chlorophyll signal is not used directly to calculate cell dimensions, since the chloroplast does not necessarily occupy the entire cell.

### 3. Manual Cell Boundary Identification

The program allows the user to identify two points on the overlaid image:

1. **Cell center:** The approximate geometric center of the cell.
2. **Outer edge:** A point along the cell's outer boundary.

The physical distance between these points is used to estimate the cell radius:

\[
r=\sqrt{[(x_2-x_1)p_x]^2+[(y_2-y_1)p_y]^2}
\]

An initial circular ROI is generated using this radius.

The proposed boundary is displayed on the grayscale and overlaid images for visual verification.

The user can then accept the boundary, adjust its shape using an interactive elliptical ROI, manually trace the outer boundary using a freehand ROI, repeat the center-and-edge selection, or skip the image.

Once accepted, the selected region is converted into a binary mask, with cell pixels assigned a value of 1 and background pixels assigned a value of 0.

### 4. Cell Size Calculation

Morphological parameters are determined from the accepted binary mask.

**Projected cell area**

\[
A=N_{\mathrm{pixels}}p_xp_y
\]

where \(N_{\mathrm{pixels}}\) is the number of pixels within the binary mask and \(p_x,p_y\) are the calibrated pixel dimensions.

**Major and minor axes**

The coordinates of pixels within the binary mask are converted to physical units. The second central-moment matrix is constructed, and its eigenvalues are used to calculate the axes of the moment-equivalent ellipse:

\[
D_{\mathrm{major}}=4\sqrt{\lambda_1}
\]

\[
D_{\mathrm{minor}}=4\sqrt{\lambda_2}
\]

where \(\lambda_1\) and \(\lambda_2\) are the larger and smaller eigenvalues, respectively.

**Equivalent circular diameter**

\[
D_{\mathrm{eq}}=\sqrt{\frac{4A}{\pi}}
\]

The equivalent diameter represents the diameter of a circle having the same projected area as the segmented cell.

**Aspect ratio**

\[
AR=\frac{D_{\mathrm{major}}}{D_{\mathrm{minor}}}
\]

The program also calculates cell orientation using the principal eigenvector of the second central-moment matrix.

### 5. Selection of Representative Z-Plane

Each selected Z-plane is analyzed independently. For each series, the program identifies the accepted optical section with the largest projected cell area:

\[
Z_{\mathrm{representative}}
=\underset{Z}{\operatorname{arg\,max}}\;A(Z)
\]

Cell masks intersecting the image boundary are excluded from this selection.

The largest-area optical section is used as an approximation of the central cell cross-section. The corresponding major-axis length, minor-axis length, and equivalent diameter are retained for comparison across imaging series.

### 6. Statistical Analysis and Visualization

The cell-size measurements can be imported into a separate MATLAB statistical analysis script or entered manually as arrays corresponding to different experimental conditions.

The following statistics are calculated:

- Mean cell size
- Median cell size
- Sample standard deviation (SD)
- Standard error of the mean (SEM)
- Minimum and maximum cell size
- Number of independent cell measurements

The measurements are visualized using bar graphs showing mean cell size with standard deviation error bars. Individual measurements are superimposed to illustrate variation among cells.

## Output Files

The program creates an output directory containing the following files.

| File | Description |
|---|---|
| `All_Series_All_Z.csv` | Morphological measurements for all accepted Z-planes |
| `All_Series_Central_Cell_Size.csv` | Measurements from the largest-area accepted Z-plane of each series |
| `Largest_Area_Central_Z.csv` | Representative cell measurements for an individual series |
| `*_Cell_Size.png` | Annotated grayscale image with cell boundary, fitted ellipse, major/minor axes, and equivalent-diameter circle |
| `*_Cell_Size.fig` | Editable MATLAB figure |
| `*_Mask.png` | Binary cell mask |
| `*_Mask.mat` | Binary mask, calibration, and measurement data |
| `Cell_Size_All_Series_BarGraph.png` | Comparison of cell dimensions across imaging series |
| `Cell_Size_All_Series_BarGraph.fig` | Editable MATLAB bar graph |

## Limitations

The program determines cell dimensions from two-dimensional optical sections rather than reconstructing the complete three-dimensional cell geometry.

The accuracy of the measurements depends on correct identification of the outer cell boundary. Manual verification is therefore incorporated into the segmentation procedure.

The largest-area Z-plane is treated as a representative cell cross-section but is not guaranteed to coincide with the exact geometric center of the cell. Additionally, multiple Z-plane measurements of the same cell should not be treated as independent biological replicates.

## Summary

This MATLAB program provides an interactive, reproducible workflow for measuring *Chlamydomonas reinhardtii* cell morphology from LIF microscopy datasets. By combining chlorophyll-guided cell localization with manually verified grayscale boundary selection, the program avoids relying solely on chloroplast fluorescence for cell-size determination.

The resulting binary masks are used to calculate physically calibrated morphological parameters, while the largest-area optical section from each series provides a consistent basis for comparing cell dimensions across experimental conditions.
