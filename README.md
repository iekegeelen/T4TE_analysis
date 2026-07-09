Study:
 TMS-EEG investigation of whether pre-stimulus µ-alpha oscillation phase modulates cortical and 
 Corticospinal excitability in left primary motor cortex (M1).


Stage 0: Preprocessing (Pipeline runs on raw data)
T4TE_preprocessing_pipeline_v9_BEL_S01.m - Complete TMS-EEG preprocessing pipeline. Produces trial-level EEG data with artifact removal.


Stage 1: Individual Alpha Frequency Estimation
T4TE_RS_preprocessing_BEL_S01.m - Preprocesses resting-state EEG recordings (eyes-open, pre TMS blocks).

Stage 2: MEP Analysis
T4TE_MEP_analysis_BEL_S01.m - Interactive trial-by-trial validation of motor evoked potentials (MEPs) from FDI EMG.

Stage 3: Phase Extraction
T4TE_phase_MEP_analysis_v7.m - Analyzes relationship between pre-stimulus µ-alpha phase and MEP amplitude.

Stage 4: TEP Extraction & Analysis
T4TE_TEP_extraction_v4.m - Extracts trial-averaged TEP waveforms per condition per subject.

Stage 5: Figure Generation
T4TE_TEP_waveform_v6.m - Generates waveform plots overlaying phase-binned TEP responses.
T4TE_TEP_topomap_v6.m - Generates topographic maps of TEP components across phase bins.
T4TE_TEP_grandaverage_IAF.m - Generates grand-averaged TEP waveforms collapsed across phase bins, separated by intensity.
T4TE_S5_confound_table.m - Generates supplementary table quantifying pre-stimulus alpha power as a potential confound.
T4TE_power_phase_figure.m - Visualizes the relationship between pre-stimulus alpha power and alpha phase.
T4TE_prestim_power.m - Analyzes pre-stimulus alpha power as a potential predictor of MEP/TEP amplitude (independent of phase).


*Phase Binning*
Trough: −45° to +45°
Rising: +45° to +135°
Peak: +135° to −135°
Falling: −135° to −45°



*Statistical Testing*
Circular-linear correlation (circ_corrcl): uses all trials without discretization; Wilcoxon signed-rank test at group level
One-way RM-ANOVA: tests phase bin effect; Bonferroni-corrected post-hoc for multiple comparisons
Confound checks: pre-stimulus power is not a significant predictor of MEP/TEP amplitude

