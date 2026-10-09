# Gauss+ Three-Factor Affine Term Structure Model

Implementation of the Gauss+ model from Tuckman & Serrat, *Fixed Income Securities*, Chapter 9 and Appendix A9.2.

This is not an exact replication because some details for the fitting of the curve are not disclosed in the book. I also added Kalman filter as an alternative to the exact extraction of the factors, but note that this is not an appropriate implementation becaause it uses Q-dynamics in the transition, i.e. that risk premia are zero.

Finally, the book has multiple typos in the annex and the main chapter. I tried to higlight all of them here.