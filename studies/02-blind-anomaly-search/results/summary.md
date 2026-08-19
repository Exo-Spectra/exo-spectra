# Phase 5 — model-independent anomaly search over the archive

## Tier A — per-spectrum structure statistics
- Transmission: 615 spectra (>=5 usable points), structured @FDR0.01: 288 (47%), median chi2_red_flat 1.60
- Eclipse: 126 spectra (>=5 usable points), structured @FDR0.01: 83 (66%), median chi2_red_flat 3.62

Top 15 by structure amplitude (chi2_red vs flat):

    pl_name    spec_type                                                      instrument               authors  n_used  chi2_red_flat  chi2_red_slope  acf_lag1_snr
HD 189733 b      Eclipse                            Mid-Infrared Instrument (MIRI) - LRS    Inglis et al. 2024      33     740.138832       22.536626      4.002765
   KELT-9 b Transmission Potsdam Echelle Polarimetric & Spectroscopic Instrument (PEPSI)    Cauley et al. 2019      19     588.799349      421.551626      0.200762
HD 189733 b      Eclipse                            Mid-Infrared Instrument (MIRI) - LRS    Inglis et al. 2024      33     444.352538       15.371160      2.932624
HD 189733 b      Eclipse                                   Near Infrared Camera (NIRCam)     Zhang et al. 2025      53     359.607811       47.202548      6.764625
  WASP-43 b      Eclipse                                  Mid-Infrared Instrument (MIRI)      Bell et al. 2024      14     227.630245       18.267957      0.655961
  WASP-43 b      Eclipse                                  Mid-Infrared Instrument (MIRI)      Bell et al. 2024      14     164.339897       15.841488      0.984706
  WASP-43 b      Eclipse                                  Mid-Infrared Instrument (MIRI)      Bell et al. 2024      14     153.538539       17.502897      1.329649
  WASP-43 b      Eclipse                                  Mid-Infrared Instrument (MIRI)      Bell et al. 2024      14     144.487438       15.296891      0.939964
  WASP-43 b      Eclipse                                  Mid-Infrared Instrument (MIRI)      Bell et al. 2024      14     135.119059       15.056662      0.991290
 WASP-121 b      Eclipse         Near Infrared Imager and Slitless Spectrograph (NIRISS) Pelletier et al. 2026     480      92.291519        3.080485     12.193311
 HD 80606 b      Eclipse                    Near Infrared Spectrograph (NIRSpec) - G395H    Sikora et al. 2025      24      60.081860        1.567216     -0.154061
  WASP-69 b      Eclipse                                   Near Infrared Camera (NIRCam)  Schlawin et al. 2024      51      49.619927       12.754213      5.730861
  WASP-18 b      Eclipse                                             Wide Field Camera 3 Arcangeli et al. 2018      14      49.252853        3.132027      0.970502
WASP-77 A b      Eclipse                    Near Infrared Spectrograph (NIRSpec) - G395H    August et al. 2023     150      48.391576        9.520204     10.891867
 WASP-107 b Transmission                                   Near Infrared Camera (NIRCam)    Murphy et al. 2024      30      45.943004       12.319188      3.244116

## Tier C — point anomalies (|z_local| > 4)
- 243 anomalous points in 93 spectra; 149 confirmed by >=1 other spectrum of the same planet, 67 contradicted
- instrument hotspots (>=3 anomalies from >=2 planets in a 1% wavelength bin): 12 — recurring bins are suspected instrument/reduction systematics, not astrophysics

Top 15 point anomalies by |z|:

     pl_name    spec_type                                                      instrument             authors    wave    z_local  n_other_specs  n_confirming  n_contradicting
    KELT-9 b Transmission Potsdam Echelle Polarimetric & Spectroscopic Instrument (PEPSI)  Cauley et al. 2019 0.65628  97.052632              0             0                0
    KELT-9 b Transmission Potsdam Echelle Polarimetric & Spectroscopic Instrument (PEPSI)  Cauley et al. 2019 0.48613  68.315789              0             0                0
   WASP-17 b Transmission         Near Infrared Imager and Slitless Spectrograph (NIRISS)   Louie et al. 2025 2.05820 -25.464576              1             1                0
    KELT-9 b Transmission Potsdam Echelle Polarimetric & Spectroscopic Instrument (PEPSI)  Cauley et al. 2019 0.49576 -23.440000              0             0                0
    KELT-9 b Transmission Potsdam Echelle Polarimetric & Spectroscopic Instrument (PEPSI)  Cauley et al. 2019 0.51690  20.896552              0             0                0
  TOI-5205 b Transmission                    Near Infrared Spectrograph (NIRSpec) - PRISM   Cañas et al. 2026 1.59334 -13.373494              7             3                2
  TOI-5205 b Transmission                    Near Infrared Spectrograph (NIRSpec) - PRISM   Cañas et al. 2026 1.10667 -10.614583              7             7                0
    KELT-9 b Transmission Potsdam Echelle Polarimetric & Spectroscopic Instrument (PEPSI)  Cauley et al. 2019 0.48242  -9.655172              0             0                0
    KELT-9 b Transmission Potsdam Echelle Polarimetric & Spectroscopic Instrument (PEPSI)  Cauley et al. 2019 0.53166   9.379310              0             0                0
   WASP-43 b      Eclipse                                  Mid-Infrared Instrument (MIRI)    Bell et al. 2024 5.25000  -9.088161              4             4                0
   WASP-43 b      Eclipse                                  Mid-Infrared Instrument (MIRI)    Bell et al. 2024 5.25000  -8.639037              4             4                0
   WASP-43 b      Eclipse                                  Mid-Infrared Instrument (MIRI)    Bell et al. 2024 5.25000  -8.552995              4             4                0
TRAPPIST-1 b Transmission                    Near Infrared Spectrograph (NIRSpec) - PRISM Rathcke et al. 2025 1.58096  -8.507617              1             0                0
 V1298 Tau b Transmission                            Near Infrared Spectrograph (NIRSpec)   Barat et al. 2025 5.12040  -8.484032              0             0                0
   WASP-43 b      Eclipse                                  Mid-Infrared Instrument (MIRI)    Bell et al. 2024 5.25000  -8.450855              4             4                0

Instrument hotspots:

                                             instrument  n_anomalies  n_planets  n_spectra  wave_lo  wave_hi    mean_z
                   Near Infrared Spectrograph (NIRSpec)            5          2          2  5.02637  5.12695 -0.124803
                   Mid-Infrared Instrument (MIRI) - LRS            3          2          3  6.62500  6.72000 -2.636800
Near Infrared Imager and Slitless Spectrograph (NIRISS)            3          2          2  0.98140  0.98414 -1.790690
Near Infrared Imager and Slitless Spectrograph (NIRISS)            3          2          3  0.76294  0.76685  1.545232
Near Infrared Imager and Slitless Spectrograph (NIRISS)            3          2          3  1.08175  1.08682  5.493555
Near Infrared Imager and Slitless Spectrograph (NIRISS)            3          2          2  1.12389  1.13470 -1.652611
Near Infrared Imager and Slitless Spectrograph (NIRISS)            3          2          2  2.37615  2.38604 -1.480810
Near Infrared Imager and Slitless Spectrograph (NIRISS)            3          2          2  2.17008  2.17991 -1.804002
                   Near Infrared Spectrograph (NIRSpec)            3          3          3  3.03088  3.04576  1.820433
           Near Infrared Spectrograph (NIRSpec) - G395H            3          2          3  3.31269  3.34367 -4.133788
                                    Wide Field Camera 3            3          3          3  1.15600  1.16300  1.101992
                                    Wide Field Camera 3            3          3          3  1.36050  1.38005  1.331009

## Tier B — cohort shape oddballs (7 cohorts, 435 spectra)
- **E_WFC3_G141** (n=41): CoRoT-1 b (Changeat et al. 2022, amp_snr 1.6); WASP-4 b (Changeat et al. 2022, amp_snr 1.4); WASP-79 b (Changeat et al. 2022, amp_snr 1.3)
- **T_G395** (n=141): TOI-776 c (Teske et al. 2025, amp_snr 1.1); L 168-9 b (Alam et al. 2025, amp_snr 1.7); L 168-9 b (Alam et al. 2025, amp_snr 1.3)
- **T_MIRI_LRS** (n=24): WASP-43 b (Bell et al. 2024, amp_snr 3.1); GJ 1214 b (Kempton et al. 2023, amp_snr 2.2); LTT 3780 c (Rigby et al. 2025, amp_snr 1.5)
- **T_NIRISS** (n=31): LHS 1140 b (Cadieux et al. 2024, amp_snr 2.4); GJ 357 b (Taylor et al. 2025, amp_snr 0.8); TRAPPIST-1 c (Radica et al. 2025, amp_snr 1.4)
- **T_PRISM** (n=19): WASP-39 b (Rustamkulov et al. 2023, amp_snr 13.5); TRAPPIST-1 d (Piaulet-Ghorayeb et al. 2025, amp_snr 1.0); TRAPPIST-1 b (Rathcke et al. 2025, amp_snr 1.9)
- **T_STIS** (n=17): WASP-79 b (Rathcke et al. 2021, amp_snr 1.6); WASP-62 b (Alam et al. 2021, amp_snr 3.4); WASP-101 b (Rathcke et al. 2023, amp_snr 0.6)
- **T_WFC3_G141** (n=162): GJ 1132 b (Swain et al. 2021, amp_snr 3.1); GJ 1132 b (Libby-Roberts et al. 2022, amp_snr 0.9); Kepler-9 c (Edwards et al. 2023, amp_snr 1.3)

## Cross-check vs phase 3 pair test
- spectra in both analyses: 395; Spearman(oddball percentile, log max pair chi2_red) = 0.10 (p = 5.8e-02)
- plot: phase5_vs_phase3.png