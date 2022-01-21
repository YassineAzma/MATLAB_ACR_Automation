# MATLAB ACR Automation
## _Introduction_

This repository contains scripts for automation of key [ACR tests](https://www.acraccreditation.org/-/media/acraccreditation/documents/mri/largephantomguidance.pdf), with additional functionality included for calculation of both single-image and subtraction Signal-to-Noise Ratio along with Modulation Transfer Function.
In its current form, all functions can be run readily from the ACR_Main.m script, which also contains additional options that can be customised as needed.

If anything is not clear please leave a message in [Discussions](https://github.com/YassineRMH/MATLAB_ACR_Automation/discussions). If you've discovered a bug or would like to make a feature request, please submit this in [Issues](https://github.com/YassineRMH/MATLAB_ACR_Automation/issues).

Documentation will be added for scripts that are nominally finished.

## _Uniformity_

The centroid is first determined.

<p align="center">
<img src="https://user-images.githubusercontent.com/96583432/150523210-9f1959f8-1a01-4748-bb01-c4487fae9a73.png" width="500"> 
</p>

This is used to determine the location of the 195-205cm<sup>2</sup> ROI within the image:

<p align="center">
<img src="https://user-images.githubusercontent.com/96583432/150523943-f08567e7-7d15-4cbb-8ee5-7dee97528b80.png" width="500">
</p>

Within this large ROI, small ROIs of roughly 1cm<sup>2</sup> are placed within it. The ROIs with maximum and minimum mean value are used to calculate the percent integral uniformity and are then displayed in an image:

<p align="center">
<img src="https://user-images.githubusercontent.com/96583432/150523525-cd16f82e-3ea4-42a8-aafb-ddd851d59bc2.png" width="500">
</p>

The function is also capable of calculating the Peak Deviation Non-Uniformity (PDNU) and the Normalised Absolute Average Deviation Uniformity (NAADU) as defined by [NEMA](https://www.nema.org/standards/view/determination-of-image-uniformity-in-diagnostic-magnetic-resonance-images). These are not reported by the function, but can by modified to do so trivially.
