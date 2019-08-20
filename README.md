# TuningMetrics
 Repository containing several tuning curve metrics.
 
The metrics included here address the following issues: 
1) Detection of cells that show a non-trivial modulation by a stimulus, in the sense that their temporal pattern of activity may be reliable, but their average firing rate is no different during the presence or absence of a stimulus. (getVisualResponsiveness.m)
2) A parameter-free, information-based metric of stimulus tuning that is more sensitive than common metrics such as (1 – circular variance), or the orientation selectivity index. (getDeltaPrime.m)
3) A parameter-free smoothness-based metric for quantifying tuning bandwidth/sparseness. (getTuningRho.m)
