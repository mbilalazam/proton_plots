# NDLAr-2x2 Proton Sample Distributions

A focused study of proton track kinematics in the ND-LAr 2×2 detector using MiniRun6.5 MC and sandbox data CAFs collected with the NuMI beam at Fermilab. Angular and track-length distributions of reconstructed proton candidates are compared between simulation and data.
All plots produced by the macros listed below are documented in the technote at docDB:12345.

---

## Input File Lists

| File | Description |
|------|-------------|
| `MiniRun6p5.txt` | List of MiniRun6.5 MC CAF files used as input to the MC analysis script |
| `sandbox_v11.txt` | List of sandbox data CAF files used as input to the data analysis script |

---

## Analysis Scripts

| File | Description |
|------|-------------|
| `purity_dlp_multiplicity.cc` | Processes MC CAF files and writes output histograms and trees to `output.root` |
| `data.cc` | Processes data CAF files and writes output to `sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root` |

---

## Output ROOT Files

| File | Description |
|------|-------------|
| `output.root` | Output produced by `purity_dlp_multiplicity.cc` |
| `sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root` | Output produced by `data.cc` |

---

## Plotting Macros

### Style
| File | Description |
|------|-------------|
| `protoDUNEStyle.C` | ROOT plotting style macro for consistent formatting used across all plots |

### Proton Selection Efficiency
| File | Description |
|------|-------------|
| `plot_MC_proton_cosTheta_efficiency.C` | Plots proton selection efficiency as a function of cos θ, comparing the baseline selection (Nominal-MINERvA) and the final selection; includes raw counts and shape-normalized versions with a bin-by-bin efficiency ratio pad |

### Proton Kinematics — Data/MC Comparisons
| File | Description |
|------|-------------|
| `plot_dataMC_proton_cosTheta_GENIE.C` | Plots proton cos θ data/MC comparison with GENIE systematic uncertainty band |
| `plot_dataMC_proton_length_GENIE.C` | Plots proton track length data/MC comparison with GENIE systematic uncertainty band |
| `plot_dataMC_proton_g4rw_band.C` | Plots proton kinematic data/MC comparison with Geant4Reweighting (G4RW) systematic uncertainty band |
| `plot_dataMC_proton_detector_band.C` | Plots proton kinematic data/MC comparison with detector response and TPC physics systematic uncertainty band |
| `plot_dataMC_proton_combined_band.C` | Plots proton kinematic data/MC comparison with the combined systematic uncertainty band (GENIE + Geant4 + Detector) |
| `plot_dataMC_proton_cosTheta_GENIE_proposed.C` | Plots the proposed Neutrino 2026 poster plot: proton cos θ data/MC comparison with GENIE systematic uncertainty band, with both absolute counts and shape-normalized versions |

### Proton Kinematics — Stacked MC Distributions
| File | Description |
|------|-------------|
| `plot_proton_stacked.C` | Plots stacked MC distributions of proton cos θ and track length, broken down by interaction type (QE, RES, DIS, MEC, COH) and final-state particle composition (True Proton, True Pion, True Muon, True Other) |
