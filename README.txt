This matlab package contains basic code for extracting time traces from different regions in photoactivated images.

Currently (2/28/2020), only randomly scattered regions, within a particular cellular region, and within a preset distance from the activation center, have been implemented.

To run, go into example_randregions.m and step through it using Matlab's cell mode (hitting ctrl-enter or the equivalent to run one chunk of code at a time).

The code will ask for interactive input to select the activation center and the overall cellular region of interest.

It will then generate a number of randomly located circular regions within the cell, and within some maximum cutoff distance from the activated center.

Final output:
regionTraces = matrix of intensity per area for each region, for each timepoint
regcent = list of center positions for the regions
regdist = list of distances for each region relative to the activation center.
