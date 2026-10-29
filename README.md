# Sleeping Beauty transposition manuscript
This repository contains the scripts and pipelines used in the manuscript "Local chromatin, 3D genome folding and transposon cargo shape Sleeping Beauty transposition" by Moene and Eder et al.

## Introduction
In this study we investigated the features that influence Sleeping Beauty transposition. We combine previously published and unpublished transposition ('hopping') datasets and evaluate the effects of chromatin context, cargo and 3D genome on transposition efficiency and integration location.

## Repository guide
This repository contains information about
- Figures: R markdown files and some input tables to create the manuscript figures.
- Chromatin context: Config, input and output files for the chromatin context analysis. This analysis is performed using the TRIP_local_context pipeline, which can be found here: [TRIP_local_context](https://github.com/vansteensellab/TRIP_local_context)
- Tagmentation: All config files for the Tn5 based mapping of transposon integrations in the paper.
- Transposed fraction: The scripts used to determine the transposed fraction from the intermediate files of the Tagmentation analysis pipeline, and output files for all samples. 
