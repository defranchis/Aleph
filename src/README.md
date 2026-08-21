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
`/eos/user/m/mdefranc/aleph_vertex/wp1_stage1/<version_tag>/<flavour_name>.root` (`--batch` writes to `/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/<year>/<mc_type>/stage1/<version_tag>/` instead).

Fraction of events to process can be set via `--fraction <val>`, default is to process all events. 

`--year <year>` and `--MCtype <type>` are also supported command line arguments, currently we only have `1994` and `zqq` here. 

### Run on data:
```bash
fccanalysis run stage1.py -- --tag <version_tag> --doData 
```

Output files will be in: `/eos/user/m/mdefranc/aleph_vertex/wp2_data/<version_tag>/` (`--batch` writes to `/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedData/<year>/stage1/<version_tag>/` instead).

`--year` and `--fraction` is also supported as an argument here. 



### Reconstruction modules: defaults and opt-outs

A flag-less `stage1.py` runs the standalone primary-vertex fitter ([`analyzer_pvnew.h`](analyzer_pvnew.h)), the two-tier V0 module ([`analyzer_v0new.h`](analyzer_v0new.h)) and the V0-first secondary-vertex module ([`analyzer_svnew.h`](analyzer_svnew.h)), writing the `pv_*`, `v0n_*`/`v0njet_*` and `svn_*`/`svm_*` branches. On MC the V0 truth-matching branches (`truev0_*`, `v0c_*`, `v0n_class`, ...) are added automatically; they are skipped under `--doData`.

The legacy code paths remain available as opt-outs:

| flag | meaning |
| --- | --- |
| `--oldPV` | legacy PV chain exactly as implemented in `FCCAnalyses`: `get_PrimaryTracks` + `VertexFitter_Tk`, and the origin-referenced track pre-selection instead of the beamspot-referenced one. No `pv_*` flag branches. Note that the beamspot constraint of the `get_PrimaryTracks` selection fit is passed in mm while its track parameters are read in cm, so that constraint is off by a factor 1000 and is effectively absent; the final `VertexFitter_Tk` fit is unaffected. |
| `--oldSV` | drop the standalone SV module: no `svn_*`/`svm_*` branches. |
| `--oldV0` | drop the two-tier V0 module: no `v0n_*`/`v0njet_*` and no V0 truth branches. Implies `--oldSV`, since the SV finder consumes the tight-V0 track veto. |

### Primary-vertex options

| flag | default | meaning |
| --- | --- | --- |
| `--pvchi2` | `5.0` | chi2max for track compatibility in the PV fit (lower = fewer tracks claimed primary = more secondaries). |
| `--pvIPWindow D0MAX Z0MAX` | `0.75 2.0` | PV pre-selection `\|D0\|`, `\|Z0\|` upper bounds [cm]. |
| `--pvBSWidth SX SY SZ` | `200 100 2` | beamspot-constraint widths for BOTH PV fits: `SX`, `SY` in um, `SZ` in cm. |
| `--excludeRuns RUN [RUN ...]` | none | data only: veto these runs before any selection (`eventsProcessed` still counts the raw input). |

### Secondary-vertex options

The SV cut values live only in the `stage1.py` argparse defaults, which are the production values; [`analyzer_svnew.h`](analyzer_svnew.h) takes them all as explicit arguments.

| flag | default | meaning |
| --- | --- | --- |
| `--svMaskMode` | `1` | V0-track masking: 0 none, 1 tight-claimed, 2 all claimed. |
| `--svChi2` | `10` | normalised vertex chi2 cut (seed and growth). |
| `--svDisLo` / `--svDisHi` | `0.03` / `3` | PV displacement window [cm]. |
| `--svSigL` | `0.10` | longitudinal vertex sigma guard [cm]. |
| `--svTrkChi2` | `5` | per-track chi2 contribution cap (<=0 off). |
| `--svCosPoint` | `0.7` | minimum SV cosPointing. |
| `--svMaxTrk` | `8` | maximum tracks per SV candidate. |
| `--svGrowShift` | `0` | max fitted-vertex displacement per growth step [cm]; 0 = off. |
| `--svClaimMode` | `0` | seed ordering: 0 best-chi2 first, 1 densest first, 2 grow-all-then-claim. |

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

### STAGE 2:

Don't touch stage2.py!
Open up stage2_all.py and change the desired input and output directories. 
Set the number of cpus.
Now you can decide if you want to divide each flavor into multiple files, then you can change the argument `n_final_files`.
Run it with nightlies.





