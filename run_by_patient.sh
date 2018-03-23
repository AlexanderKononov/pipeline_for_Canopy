#!/bin/bash

SNA = $1
CNA = $2
work_dir = $3
project_name = $4
norm_sampl = $5

python file_adapter_to_Canopy.py $SNA $CNA $work_dir $norm_sampl

/apps/statistics2/R-3.4.3/bin/Rscript pipeline_Canopy_with_clustering.R $work_dir $project_name