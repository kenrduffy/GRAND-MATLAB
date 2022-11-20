# GRAND-MATLAB
Guessing Random Additive Noise Decoding (GRAND)

Subject to license:
"GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf"

Non-parallelized MATLAB implementations of: GRAND (hard detection); basic ORBGRAND (soft detection); 1-line ORBGRAND (soft detection).

Simulation is setup and run with GRAND_Code/driver_GRAND.m

Sample output is in RESULTS, and sample plots from those results can be made with MAKE_FIGS/driver_sample_figs.m

Note that for an [n,k] code, where k information bits become n coded bits, GRAND algorithms accurately and efficiently decode codes where n-k is moderate. This MATLAB implementation is solely intended to be instructive and is not parallelised, even though highly-parallelised implementations are possible. As a result, obtaining the full performance of a code with n-k>16 may prove time-consuming with the present implementation.

The following should be cited in association with results from this code.

K. R. Duffy, J. Li, and M. Medard, "Capacity-achieving guessing random additive noise decoding," IEEE Trans. Inf. Theory, vol. 65, no. 7, pp. 4023–4040, 2019.

A. Riaz, V. Bansal, A. Solomon, W. An, Q. Liu, K. Galligan, K. R. Duffy, M. Medard and R. T. Yazicigil, "Multi-code multi-rate universal maximum likelihood decoder using GRAND," Proceedings of IEEE ESSCIRC, 2021.

K. R. Duffy, “Ordered reliability bits guessing random additive noise decoding," Proceedings of IEEE ICASSP, 2021, pp. 8268–8272.

K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing random additive noise decoding,” IEEE Trans. Signal Process., vol. 70, pp. 4528-4542, 2022.

For further details on GRAND, see: https://www.granddecoder.mit.edu/
