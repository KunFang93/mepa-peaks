# MEPA (Multi-Experiment Peak Aggregation)

*Shift-calibrated pooling of broad peaks across labs.*

Build a pooled broad-peak BED from three labs profiling the same mark in comparable samples. Step 3 builds within-lab consensus with MSPC v6.1.0. Step 4 combines lab-level peak calls, rescues loci with reproducible signal, and selects rescue thresholds using a within-chromosome shift test. A lab with one replicate contributes its full primary MACS3 peak set.

![MEPA workflow from MACS3 broadPeak and BAM or midpoint-count BigWig through MSPC, shift calibration, and pooled BED export](docs/workflow.png)

## Install

Use Python 3.11 or newer. Install this repository in a virtual environment. Download the MSPC v6.1.0 executable for your platform and put it on your path or pass its location with `--mspc`. [MSPC v6.1.0 release](https://github.com/Genometric/MSPC/releases/tag/v6.1.0)

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
mepa-peaks --help
```

## Inputs

All coordinates must use the same genome assembly and chromosome names. Supply a tab-delimited sample sheet with one row per treatment replicate or IgG control. Row order defines the lab order used for shift pairs. File paths are absolute or relative to the sample sheet. [Example with two paired labs and one singleton](examples/samples.tsv)

| Column | Treatment row (`target`) | Control row (`igg`) |
| --- | --- | --- |
| `lab`, `sample`, `role` | Lab ID, unique sample ID, and `target` | Matched lab ID, unique sample ID, and `igg` |
| `broadpeak` | Primary nine-column MACS3 broadPeak, plain text or gzip. MSPC reads column 8. | Leave blank |
| `permissive_broadpeak` | Optional second MACS3 call set used for the ±500 bp check. A blank field uses the primary set. | Leave blank |
| `bam` | Paired-end BAM for fragment counting | Matched paired-end IgG BAM |
| `bigwig` | Alternative 250 bp raw fragment-midpoint count BigWig | Alternative count BigWig for an IgG control |

Step 3 needs the primary peaks. Step 4 needs one signal source per sample. If both sources are supplied, **BAM takes precedence**. Use an empty `bam` field to select BigWig. BAM is read sequentially, so a BAM index is optional. Multiple IgG rows for a lab are pooled by summing midpoint counts and fragment depths. A lab without an IgG row uses its treatment flanks for background.

A count BigWig stores the **integer number of fragment midpoints per 250 bp tile**, anchored at chromosome position 0. The terminal tile can be shorter and omitted tiles count as zero. The library depth is the sum over the selected chromosomes. Use the following command to generate a compatible track, keeping MAPQ 0. Ordinary CPM coverage tracks contain a different quantity and need the original BAM for this conversion.

```bash
mepa-peaks export-midpoints \
  --genome genome.tsv --bam data/A_rep1.bam --output data/A_rep1.midpoints.bw
```

The genome TSV has columns `chrom`, `length`, and `split`. Order rows as desired for the BED. Only the listed chromosomes contribute to counts, normalization, and output. At least one chromosome needs the `calibration` label for Step 4.

| `split` | Role |
| --- | --- |
| `calibration` | Select the rescue thresholds |
| `holdout` | Report a separate shift-test result with the selected thresholds |
| `output_only` | Call peaks using the selected thresholds |

The included mm39 example uses odd autosomes for calibration, even autosomes other than chr2 for holdout, and chr2, X, and Y for output only. The largest default shift is 875 kb and calibration excludes 875 kb at either end. For shorter chromosomes, pass six distinct pairs with repeated `--shift-kb LEFT RIGHT` arguments. Use displacements large enough to break the local alignment of peak signals. [mm39 chromosome table](examples/mm39.genome.tsv)

## Run

Run an input check, then stop after Step 3 to inspect the per-lab BEDs. Use a new output directory for each stage.

```bash
mepa-peaks check --samples samples.tsv --genome genome.tsv --mspc /path/to/mspc

mepa-peaks step3 \
  --samples samples.tsv --genome genome.tsv --mspc /path/to/mspc \
  --outdir results/step3_mspc

mepa-peaks step4 \
  --samples samples.tsv --genome genome.tsv --step3-dir results/step3_mspc \
  --outdir results/step4_pool --workers 4 --target-shift-ratio 0.10
```

To run both stages in one invocation, give `run` a parent directory. The stage outputs go into `step3_mspc/` and `step4_pool/` within it.

```bash
mepa-peaks run \
  --samples samples.tsv --genome genome.tsv --mspc /path/to/mspc \
  --outdir results --workers 4
```

## Method

### Step 3 · Within-lab consensus

For a replicated lab, convert each primary broadPeak to BED5 with `−log10(p)` from column 8 as the score. Cap the score at 300 to avoid underflow during MSPC's p-value conversion. Run MSPC in biological-replicate mode with overlap required in all replicates (`-c 2` for two replicates). The weak p-value cutoff is 0.01, the stringent cutoff and combined Fisher cutoff are both $10^{-5}$, and the Benjamini–Hochberg alpha is 0.05. When there are several overlapping peaks in a partner replicate, use the lowest p-value.

```bash
mspc -i rep1.bed rep2.bed -r bio -c 2 -w 1e-2 -s 1e-5 -g 1e-5 \
  -a 0.05 -m lowest -d 4 --excludeHeader -o consensus_output
```

Export the MSPC consensus coordinates as the lab-level BED. A singleton bypasses MSPC and uses its primary MACS3 coordinates and peak IDs. [MSPC input and arguments](https://genometric.github.io/MSPC/docs/cli/args/)

### Step 4 · Direct support, signal rescue, and extension

**Direct cores.** Form the union of exact overlaps of calls from at least two labs. Retain continuous shared spans of at least 500 bp. Also accept a 250–499 bp overlap when it covers at least half of the shorter lab-level peak.

**Fragment counts and enrichment.** Count one fragment per primary proper read 1 with a mapped mate on the same chromosome and nonzero template length. Exclude duplicate, QC-failed, secondary, and supplementary read 1 records. Retain every MAPQ, including 0. Place each fragment at the leftmost mate position plus half the absolute template length, rounded down. Count midpoints in seven 250 bp tiles centered on each focal tile (1,750 bp).

For each replicate, estimate the expected count from the combined per-bin mean in the flanks (tiles 11–42 on each side), scaled to the size of the count window. For optional IgG, sum the matched IgG window counts, add one pooled pseudocount, and multiply by the ratio of treatment depth to pooled IgG depth. The expected count $\mu$ is the maximum of the flank expectation, the IgG expectation if supplied, and 0.5. Near chromosome ends use available bins. With zero available flank bins the flank estimate is zero. The signed Poisson deviance score for observed count $k$ is

$$
z = s\sqrt{2[k\ln(k/\mu) - (k-\mu)]}
$$

The sign $s$ is $+1$ when $k > \mu$, $-1$ when $k < \mu$, and $0$ when $k = \mu$. Use $k\ln(k/\mu)=0$ when $k=0$. The **lowest replicate score** is the lab score. For example, 12 observed fragments against 10 expected gives $z\approx 0.61$.

**Rescue a one-lab call.** The center of a candidate tile must have exactly one lab-level call and the calling lab must meet the weak score cutoff. A supporting lab needs a permissive MACS3 peak within 500 bp of that center in **every replicate**. The lab's shape must also correlate with the calling lab's shape. For shape, sum counts in a 750 bp window, convert to fragment CPM per replicate, average replicate profiles equally within each lab, and take the square root. Compute Pearson $r$ on 17 points (4,250 bp in tiles). Keep the best $r$ across shifts of −500, −250, 0, +250, and +500 bp. Each qualifying supporting lab contributes one vote. Accept either one strong supporting lab or both other labs at the weak cutoff. Require at least two adjacent passing tiles with the same calling lab. Clip the rescued runs to the union of the primary replicate peaks and retain clipped spans of at least 500 bp.

**Extend boundaries.** Two or more labs must reach the weak enrichment cutoff in each tile, with evidence from a lab-level call or permissive calls near the center in every replicate. Grow through runs of at least two such tiles that overlap a core after clipping to the primary replicate peak union. A gap in support or in that union breaks extension.

### Auto cutoff · Shift test

Evaluate the grid on the calibration chromosomes.

| Parameter | Candidate values |
| --- | --- |
| Weak signed Poisson score | 0.50, 0.75, 1.00, 1.25 |
| Strong signed Poisson score | 2.00, 2.25, 2.50, 2.75, 3.00, 3.50 |
| Minimum Pearson $r$ | 0.55, 0.60, 0.65, 0.70, 0.75, 0.80 |

1. Count observed passing seed tiles after the continuity check and before clipping. Candidate centers start with a single lab-level call and lie outside the retained 500 bp direct cores.
2. Hold the calling lab fixed. Move each other lab's score, permissive-call mask, shape profile, and lab-level call mask together within that chromosome. For each anchor, assign the two offsets to the other two labs in sample-sheet order. Use six offset pairs in kb `(211, -389)`, `(-227, 457)`, `(517, 875)`, `(-571, -821)`, `(383, -679)`, and `(-419, 733)`. Compare the same interior positions in all seven tests to avoid wrapping across chromosome ends. After shifting, require one lab-level call at the focal tile under the shifted masks as well.
3. For each shift, multiply the count of passing seed tiles by $\text{observed eligible tiles}/\text{shifted eligible tiles}$. Divide the mean of these six adjusted counts by the observed passing tile count.
4. Among settings with a defined ratio at or below 0.10, choose the setting with the most observed passing tiles. Break ties by higher weak score, higher strong score, then higher correlation. Evaluate the selected setting independently on the holdout chromosomes when supplied.

For example, 26,190 observed seed tiles and a mean of 1,814.53 adjusted shifted tiles give a ratio of 0.0693. The reference run selected weak 0.5, strong 2.0, and $r\ge 0.55$. The ratio was 6.9% on calibration and 6.7% on holdout. Interpret the ratio as an accidental-support diagnostic for **signal-rescue tiles**. Estimating peak-level biological error requires independent truth data. If no grid setting qualifies, Step 4 reports `direct_only` and exports the direct overlap cores.

## Outputs

| File | Contents |
| --- | --- |
| Step 3 `peak_beds/*.bed` | One lab-level BED4 per lab |
| Step 3 `peak_summary.tsv` and `step3.json` | Lab counts, MSPC version, and input and output peak checksums |
| Step 3 `mspc/` and `mspc_inputs/` | Native MSPC files, commands, logs, and scored inputs |
| Step 4 `pooled_peaks.bed` | Final 0-based half-open pooled peaks |
| Step 4 `pooled_peak_evidence.tsv` | Per-peak bases from direct overlap, short overlap, signal rescue, and shoulder extension |
| Step 4 `thresholds.json` | Selected cutoffs, shift configuration, and calibration and holdout metrics |
| Step 4 `signal_threshold_calibration.tsv` | All grid scores, grouped by split and anchor lab |
| Step 4 `count_summary.tsv`, `peak_inputs.tsv`, `pooling_summary.json` | Signal source and depth, peak inputs, and final peak counts and spans |

## Upstream MACS3 settings for the reference run

**Applicability.** MEPA was originally developed for broad H3K9me3 peaks in mouse embryonic stem (mES) cells. For transcription factors and other histone modifications with different peak widths or signal shapes, test assay-specific MACS3 settings on your own data. Use replicate agreement and matched controls, where available, to check that MEPA's overlap lengths and signal windows are appropriate before interpreting pooled peaks.

The Step 3 and Step 4 commands expect existing peak files. The reference study used MACS3 3.0.2 on properly paired, primary, nonduplicate nuclear alignments, including MAPQ 0. The primary call used `-f BAMPE -g mm --broad --keep-dup all --llocal 20000 --max-gap 500 -q 0.05 --broad-cutoff 0.10`. The separate permissive call used `-f BAMPE -g mm --broad --keep-dup all -p 0.001 --broad-cutoff 0.05` with the caller's default local window and gap. In each call, supply the matching IgG BAMs via `-c` where available. The shift test calibrates the **Step 4 score and shape cutoffs**. MACS3 and MSPC use their specified peak-calling settings. [MACS3 callpeak documentation](https://macs3-project.github.io/MACS/docs/callpeak.html)

## Validation

A replay of the reference mm39 analysis reproduced all 43,309 pooled peak coordinates, the per-peak evidence lengths, and all 1,152 calibration grid rows. The included synthetic test runs MSPC and the shift calibration on three labs and compares the BED, evidence and cutoffs obtained from BAMs with those from their exported BigWigs.

```bash
MEPA_MSPC=/path/to/mspc python -m unittest discover -s tests -v
```
