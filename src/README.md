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
`/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/<year>/<mc_type>/stage1/<version_tag>/<flavour_name>.root`; set `ALEPH_OUT_DIR=<dir>` to write to `<dir>/wp1_stage1/<version_tag>/` instead (e.g. when `/eos/experiment` is read-only).

Fraction of events to process can be set via `--fraction <val>`, default is to process all events. 

`--year <year>` and `--MCtype <type>` are also supported command line arguments, currently we only have `1994` and `zqq` here. 

### Run on data:
```bash
fccanalysis run stage1.py -- --tag <version_tag> --doData 
```

Output files will be in the working directory (`--batch`: `/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/<year>/stage1/<version_tag>/`); set `ALEPH_OUT_DIR=<dir>` to write to `<dir>/wp2_data/<version_tag>/`.

`--year` and `--fraction` is also supported as an argument here. 

### Run on batch:
```
fccanalysis submit stage1.py -- --tag VXX-XX --MCflavour X --batch --chunks X
```

### Reconstruction modules: defaults and opt-outs

A flag-less `stage1.py` runs the standalone primary-vertex fitter ([`analyzer_pvnew.h`](analyzer_pvnew.h)), the two-tier V0 module ([`analyzer_v0new.h`](analyzer_v0new.h)) and the V0-first secondary-vertex module ([`analyzer_svnew.h`](analyzer_svnew.h)), writing the `pv_*`, `v0n_*`/`v0njet_*` and `svn_*`/`svm_*` branches. On MC the V0 truth-matching branches (`truev0_*`, `v0c_*`, `v0n_class`, ...) are added automatically; they are skipped under `--doData`.

The legacy code paths remain available as opt-outs:

| flag | meaning |
| --- | --- |
| `--oldPV` | legacy PV chain exactly as implemented in `FCCAnalyses`: `get_PrimaryTracks` + `VertexFitter_Tk`, and the origin-referenced track pre-selection instead of the beamspot-referenced one. No `pv_*` flag branches. Note that the beamspot constraint of the `get_PrimaryTracks` selection fit is passed in mm while its track parameters are read in cm, so that constraint is off by a factor 1000 and is effectively absent; the final `VertexFitter_Tk` fit is unaffected. |
| `--oldSV` | drop the standalone SV module: no `svn_*`/`svm_*` branches. |
| `--oldV0` | drop the two-tier V0 module: no `v0n_*`/`v0njet_*` and no V0 truth branches. Implies `--oldSV`, since the SV finder consumes the tight-V0 track veto. |

The selection is not configurable from the command line: every tuned value is a named `constexpr` in the header that uses it, and that declaration is its only definition. The paragraphs below describe the selections in words; the numbers quoted are those declarations.

### The primary-vertex selection

Tracks enter the fit through a pre-selection window on the impact parameters, `|D0| < 0.75 cm` and `|Z0| < 2 cm` (`PVN_D0_MAX`, `PVN_Z0_MAX`), referenced to the run beamspot — or to the origin under `--oldPV`, which leaves it off-centre by the beamspot offset in data. Track/vertex compatibility is then judged at `chi2max = 5` (`PVN_CHI2_MAX`); lowering it claims fewer tracks as primary and so leaves more to the secondary finders.

Both PV fits — the selection fit that prunes the track list and the final position fit — share one beamspot constraint, of Gaussian widths 200 um transverse in x, 100 um in y and 2 cm along the beam (`PVN_BS_SIGMA_X/Y/Z`, declared in cm). All of these live in [`analyzer_pvnew.h`](analyzer_pvnew.h), which is loaded for either PV chain because the legacy chain reads the same numbers, rescaled to its own unit convention.

### The secondary-vertex selection

The SV values live in [`analyzer_svnew.h`](analyzer_svnew.h) as the `SVN_*` constants, which are also the default arguments of `findSVs`, so `stage1.py` passes none of them. A seed or a growth step is accepted at normalised vertex `chi2 < 10` (`SVN_CHI2`) with a per-track chi2 contribution capped at 5 (`SVN_TRK_CHI2`; <=0 disables the cap). A candidate is kept if it sits between 0.03 and 3 cm from the PV (`SVN_DIS_LO`, `SVN_DIS_HI`), its longitudinal vertex sigma stays under 0.10 cm (`SVN_SIGL_MAX`), and it points back to the PV with `cosPointing > 0.7` (`SVN_COS_POINT`). Growth stops at 8 tracks per candidate (`SVN_MAX_TRK`); the per-step fitted-vertex displacement guard (`SVN_GROW_SHIFT`) is off. Seeds are claimed best-chi2 first (`SVN_CLAIM_MODE`).

Two SV collections are written from the same event: `svn_*` runs on the tracks left after the tight-claimed V0 tracks are masked out (`SVN_MASK_MODE`), and `svm_*` is the unmasked control twin (`SVN_MASK_NONE`), for studying the V0/SV interplay.

### Other options

`--excludeRuns RUN [RUN ...]` (data only) vetoes the listed run numbers before any selection; `eventsProcessed` still counts the raw input.

Environment overrides for local (non-`--batch`) runs: `ALEPH_RECLUS_DIR=<dir>` reads the input files from `<dir>` instead of `/eos/experiment` (a re-clustered copy lifts the RDataFrame thread cap set by the few TTree clusters of the raw files); `ALEPH_OUT_DIR=<dir>` writes the output under `<dir>` (see above).

### The two-tier V0 module

Standalone V0 (Ks/Λ) reconstruction in [`analyzer_v0new.h`](analyzer_v0new.h). It is *V0-first* by design: its track claims are meant to be consumed downstream (secondary-vertex finding runs on the unclaimed tracks), so the module optimises the correctness of each claim, not just the candidate list.

**Candidate building.** All opposite-charge secondary track pairs are vertexed with a single consistent fit (`VertexFitter_Tk`); every downstream quantity (momenta, invariant masses under both hypotheses, Armenteros–Podolanski (AP) variables, pointing) is derived from the refitted momenta at that vertex — there is no second fit.

**Two selection tiers, evaluated per hypothesis (Ks and Λ):**

- **Tight** — the adopted physics selection: mass window; momentum-tiered pointing cut (separate Ks and Λ ladders, `ksPointThr` / `lamPointThr`); a qT veto against photon conversions (Λ only); and a resolution-scaled AP-band cut around the exact kinematic locus — the band half-width follows the measured σ_ell(p) of each species (`ksBandThr`, `lamBandThrTight`; the Λ band is floored at low p and capped at the nominal ramp edge), plus common fit-quality (χ²) and displacement requirements.
- **Loose** — the ML-training tier: same windows/χ²/displacement, but flat pointing, a wider Ks AP band, a Λ AP band equal to a fixed fraction (0.8) of the ramp half-width floored at the tight band, and a relaxed Λ qT veto. It is a strict superset of tight and is what gets stored, so *any* tighter selection (including the historical ones) can be re-derived offline from any production.

**Hypothesis arbitration.** A pair passing both hypotheses is booked as the one whose invariant mass is closer to its window centre (normalised by the window half-width).

**Exclusive claiming, tight first.** Candidates claim their tracks exclusively in quality order: all tight candidates claim before any loose one, and within a tier the best-χ² candidate claims first. A track is claimed once; later candidates using it are dropped. This preserves the tight-only output exactly regardless of the loose tier, and defines the track set left to the SV finder.

**Stored flags and ML inputs.** `v0n_tight` re-derives the adopted tight package offline from the same single-source helpers (booking a candidate ≠ selecting it — variant productions still carry the adopted-package flag). Note the adopted package has changed over the module's history (p-tiered Λ pointing, σ-scaled AP band): the flag is not comparable across module versions, but any tighter selection is re-derivable offline from the stored loose tier. `v0n_bandSig` and `v0n_massSig` store the AP-band and mass cut variables as signed pulls in resolution units for training; all other cut variables (cosPointing, pointSig, qT, χ², displacement, p, invM) are stored raw.

**Study variants.** `--v0nLamPointKsTiers` (Λ pointing fully aligned to the Ks ladder) and `--v0nWideLamLoose` (loose-Λ AP band ramp edges doubled, for measuring the band tail) select wrapper configurations for systematic studies; the stored `v0n_tight` flag always encodes the adopted package.

**Output branches.** Two groups, both written by default and both dropped by `--oldV0`:

- `n_v0n_event`, `v0n_pdg`, `v0n_invM`, `v0n_alpha`, `v0n_qt`, `v0n_chi2`, `v0n_dxyz`, `v0n_p`, `v0n_px/py/pz`, `v0n_cosPointing`, `v0n_pointSig`, `v0n_tight`, `v0n_bandSig`, `v0n_massSig`, `v0n_vx/vy/vz`, the vertex-fit covariance `v0n_cov_*`, the daughter joins `v0n_trk1_origIdx`/`v0n_trk2_origIdx` and their `v0n_trk{1,2}_dEdx_{pads,wires}_{value,error}` — event-order candidate quantities, independent of jet assignment.
- `n_v0njet_jets`, `n_v0njet_ks`, `n_v0njet_lambda` and `v0njet_*` — the new candidates pushed through the same per-jet assignment and jet-relative getters as the existing `v0_*` block, so old vs new is an apples-to-apples comparison at jet level.

Two conditional additions on top: the V0→nearest-SV pointing features `v0n_svnCosPoint`, `v0n_svnPointSig`, `v0n_svnIdx` need the SV module (dropped by `--oldSV`), and the truth-matching branches (`v0n_class`, `v0n_trueidx`, `v0n_pairmult`, `v0n_trackshared`, `v0n_trk1`, `v0n_trk2`, `truev0_foundnew_*`, `truev0_*`, `v0c_*`) are MC only — they are skipped under `--doData`.

**Further utilities in the headers.** Beyond the alternative wrapper entry points reached by the study-variant flags above, [`analyzer_truth.h`](analyzer_truth.h) carries the MC truth-matching utilities (true-V0 finding, track↔MC index recovery, candidate truth classification) used to derive the tunings above; it is loaded unconditionally, because its truth-FREE candidate accessors (`candChi2`, `candDxyz`, `candP`, `candPcomp`, `candCosPointing`, `candVtxPos`) are also used on data.

### STAGE 2:

Don't touch stage2.py!
Open up stage2_all.py and change the desired input and output directories. 
Set the number of cpus.
Now you can decide if you want to divide each flavor into multiple files, then you can change the argument `n_final_files`.
Run it with nightlies.





