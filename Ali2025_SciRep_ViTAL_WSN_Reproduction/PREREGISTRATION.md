# Pre-registration (written 2026-10-12T14:10Z, before any simulation of ViTAL was run)

Target: S. J. Ali, J. Iqbal, F. Khan, Y.-K. Moon, "Vigorous technique for augmented lifetime in WSNs",
Scientific Reports 15:15400 (2025), https://doi.org/10.1038/s41598-025-96830-w (CC BY 4.0).
No code or data are public. The paper says data are available from the corresponding author on request.

These choices are fixed now and will not be changed after seeing results. Alternatives are only run as
labelled sensitivity analyses.

| ID | Item | Paper | Primary choice |
|---|---|---|---|
| P1 | Area, sink, N, E0, Rtx, Rs | 100x100 m, sink (5,5), N=100..500, 0.5 J, 15 m, 5 m | as stated |
| P2 | Packet size | 800 bits | 800 bits; 4000 bits as sensitivity |
| P3 | Eelec, EDA | not given | 50 nJ/bit, 5 nJ/bit/signal (standard LEACH values) |
| P4 | eps_fs, eps_mp | "10 pJ/bit", "0.0013 pJ/bit" (labelled receive energy) | 10 pJ/bit/m^2, 0.0013 pJ/bit/m^4 |
| P5 | p (CH fraction) | not given | 0.1 |
| P6 | w1..w4 (Eq. 20) | not given | 0.25 each; sweep as sensitivity |
| P7 | Eq. 20 distance term | Dist/Avg_dist, as printed | as printed; inverted Avg/Dist as sensitivity |
| P8 | SNR | log-normal shadowing, no parameters | Pt 0 dBm, PL(1 m) 40 dB, n=3, sigma 4 dB static per node, noise -100 dBm; SNR in dB floored at 0 |
| P9 | Eq. 22 use | "ViTAL = NV x T(n)"; Alg. 2 picks highest candidate | used as a LEACH-style threshold: node is CH if rand <= NV*T(n) |
| P10 | CDS | each node builds CDS; CHs chosen from the CDS | greedy connected dominating set on alive nodes (Rtx 15 m), recomputed when a node dies; ViTAL CH candidates restricted to it |
| P11 | HSA deployment | Eq. 9 objective, HMCR/PAR/BW from ref. 22 (not given) | maximise cr*dmin with s=N fixed; binary disk coverage on a 1 m grid; HMS 10, HMCR 0.9, PAR 0.3, BW 5, 2000 improvisations |
| P12 | Data path | single-hop LEACH | member -> nearest CH (or BS if nearer), CH -> BS single hop; no CH in a round -> nodes send directly |
| P13 | EE-LEACH, I-LEACH | described in one sentence each; originals not accessible | EE-LEACH*: T(n)*Eres/E0. I-LEACH*: only nodes with Eres >= mean alive energy; T(n)*(neigh/mean neigh)*(mean dist/dist). Both labelled as reconstructions |
| P14 | Mobile sink | random waypoint, 1-10 m/s, 10 s pause | 1 round = 1 s; same model |
| P15 | Packets to BS | not defined precisely | every packet that reaches the BS (CH aggregate = 1, direct = 1) |
| P16 | Average residual energy | time not stated | mean over all N nodes at round 1000 (paper: "each run lasts 1000 s"); round 500 also reported |
| P17 | Repetitions | 30 runs | 30 seeds per configuration |
| P18 | Improvement % | paper's convention appears to be (V-B)/V | report both (V-B)/V and (V-B)/B |
