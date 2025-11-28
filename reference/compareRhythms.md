# Run differential rhythmicity analysis

The differential rhythmicity analysis is run with a call to this
function. To execute this function, the three necessary ingredients are
the timeseries data, the experimental design and parameters to choose
and tune the method.

## Usage

``` r
compareRhythms(
  data,
  exp_design,
  lengths = NULL,
  method = "mod_sel",
  period = 24,
  rhythm_fdr = 0.05,
  compare_fdr = 0.05,
  amp_cutoff = 0.5,
  criterion = "bic",
  schwarz_wt_cutoff = 0.6,
  just_classify = TRUE,
  robust = TRUE,
  outliers = FALSE,
  longitudinal = FALSE,
  just_rhythms = TRUE
)
```

## Arguments

- data:

  A matrix of log2 expression values (if microarray), expression counts
  (RNA-seq) or normalized data (see Details).

- exp_design:

  A data.frame of the experimental design with at least two columns:
  "time" and "group".

- lengths:

  A data.frame of average transcript lengths. Only used with methods
  "deseq" and "edgeR".

- method:

  The method of analysis. It should be one of "mod_sel" for model
  selection, "dodr" for analysis using
  [DODR::dodr](https://rdrr.io/pkg/DODR/man/dodr.html), "limma" for
  linear-modeling approach based on limma, "voom" for linear-modeling
  approach for RNA-Seq using
  [limma::voom](https://rdrr.io/pkg/limma/man/voom.html), "deseq2" for
  RNA-seq analysis using DESeq2, "edger" for RNA-seq analysis using
  edgeR, and "cosinor" for simple cosinor-based analysis both for
  independent samples and repeated samples.

- period:

  The period of rhythm being tested (default = 24)

- rhythm_fdr:

  The false discovery cutoff for finding rhythmic time series (default =
  0.05)

- compare_fdr:

  The false discovery cutoff for the comparison of rhythms (default =
  0.05)

- amp_cutoff:

  The minimum peak-to-trough amplitude in log2 scale considered
  biologically relevant (default = 0.5)

- criterion:

  The criterion used for model selection. These can be "aic" or "bic"
  (default = "bic"). Only used for method = "mod_sel".

- schwarz_wt_cutoff:

  The conditional probability that the best model is the true model.
  Genes with a conditional probability smaller than this cutoff are
  deemed unclassifiable. This is only used for method = "mod_sel".
  (default = 0.4)

- just_classify:

  Boolean specifying whether genes must only be classified (TRUE) or if
  the amplitude and phases of fits should also be returned (FALSE)

- robust:

  Boolean to turn on robust computation of statistics in different
  methods (default = TRUE).

- outliers:

  Boolean specifying if weights must be computed for each sample to
  account for outliers. Only used by method = "voom".

- longitudinal:

  Boolean specifying if repeated samples from one experimental unit.
  Only used by method = "cosinor".

- just_rhythms:

  Boolean specifying whether only rhythm analysis (True, default) or
  differential expression/magnitude analysis should also be performed.

## Value

A *data.frame* with the names of the differentially rhythmic features,
the category it is classified under and optionally the rhythm parameters
of the features in each group. The differential rhythmicity categories
are **gain** of, **loss** of, **change** of, or **same** rhythms (with
respect to the reference/control group).
