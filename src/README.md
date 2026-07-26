# Aleph

Version of the FCCAnalyses code that supports command line arguments, to be able to process both data and MC with the same script. Nightlies version of the key4hep stack is required for this. 

## Setup

Folow the steps described in the general [README](../README.md) at the top level to setup the code, then just `cd src`. 
 
<!-- ```bash
git clone https://github.com/Apranikstar/Aleph.git
cd Aleph
git submodule update --init --recursive
cd FCCAnalyses
fccanalysis build -j 8
cd ../src
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh #or compile and source the FCCAnalyses module.
```

Need nightlies because updated FCCAnalyses version after this commit is needed: https://github.com/HEP-FCC/FCCAnalyses/pull/474 -->

## Stage1: Produce ntuples

Note: Change the data fraction based on your needs.

### Run on MC:
```bash
fccanalysis run stage1.py -- --tag <version_tag>  --MCflavour <flavour_index>
```

Output files will be in: 
`/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/<year>/<mc_type>/stage1/<version_tag>/<flavour_name>.root`. 

Fraction of events to process can be set via `--fraction <val>`, default is to process all events. 

`--year <year>` and `--MCtype <type>` are also supported command line arguments, currently we only have `1994` and `zqq` here. 

### Run on data:
```bash
fccanalysis run stage1.py -- --tag <version_tag> --doData 
```

Output files will be in: `/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/<year>/stage1/<version_tag>/`

`--year` and `--fraction` is also supported as an argument here. 



### The two-tier V0 module (`--newV0`)

Standalone V0 (Ks/Λ) reconstruction in [`analyzer_v0new.h`](analyzer_v0new.h), enabled with `--newV0`. It is *V0-first* by design: its track claims are meant to be consumed downstream (secondary-vertex finding runs on the unclaimed tracks), so the module optimises the correctness of each claim, not just the candidate list.

**Candidate building.** All opposite-charge secondary track pairs are vertexed with a single consistent fit (`VertexFitter_Tk`); every downstream quantity (momenta, invariant masses under both hypotheses, Armenteros–Podolanski (AP) variables, pointing) is derived from the refitted momenta at that vertex — there is no second fit.

**Two selection tiers, evaluated per hypothesis (Ks and Λ):**

- **Tight** — the adopted physics selection: mass window; momentum-tiered pointing cut (separate Ks and Λ ladders, `ksPointThr` / `lamPointThr`); a qT veto against photon conversions (Λ only); and a resolution-scaled AP-band cut around the exact kinematic locus — the band half-width follows the measured σ_ell(p) of each species (`ksBandThr`, `lamBandThrTight`; the Λ band is floored at low p and capped at the loose-band edge), plus common fit-quality (χ²) and displacement requirements.
- **Loose** — the ML-training tier: same windows/χ²/displacement, but flat pointing, widened AP bands and a relaxed Λ qT veto. It is a strict superset of tight and is what gets stored, so *any* tighter selection (including the historical ones) can be re-derived offline from any production.

**Hypothesis arbitration.** A pair passing both hypotheses is booked as the one whose invariant mass is closer to its window centre (normalised by the window half-width).

**Exclusive claiming, tight first.** Candidates claim their tracks exclusively in quality order: all tight candidates claim before any loose one, and within a tier the best-χ² candidate claims first. A track is claimed once; later candidates using it are dropped. This preserves the tight-only output exactly regardless of the loose tier, and defines the track set left to the SV finder.

**Stored flags and ML inputs.** `v0n_tight` re-derives the adopted tight package offline from the same single-source helpers (booking a candidate ≠ selecting it — variant productions still carry the adopted-package flag). Note the adopted package changed on 2026-07-27 (C′ pointing + σ-scaled band): the flag is not comparable across pre-/post-change productions, but any tighter historical selection is re-derivable offline from the stored loose tier. `v0n_bandSig` and `v0n_massSig` store the AP-band and mass cut variables as signed pulls in resolution units for training; all other cut variables (cosPointing, pointSig, qT, χ², displacement, p, invM) are stored raw.

**Study variants.** `--v0nLamPointKsTiers` (Λ pointing fully aligned to the Ks ladder) and `--v0nWideLamLoose` (widened loose-Λ AP band, for measuring the band tail) select wrapper configurations for systematic studies; the stored `v0n_tight` flag always encodes the adopted package.

### STAGE 2:

Don't touch stage2.py!
Open up stage2_all.py and change the desired input and output directories. 
Set the number of cpus.
Now you can decide if you want to divide each flavor into multiple files, then you can change the argument `n_final_files`.
Run it with nightlies.





