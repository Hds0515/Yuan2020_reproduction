# V15 final decision report

1. **A model is delivered:** `models_v15.Fixed60HighResolutionReference` has
   4800 thermal states and a frozen scope guard.
2. **COMSOL converged:** the 1200-state independent implementation compiled in
   COMSOL 6.4 and completed its transient BDF solve.
3. **Implementation consistency passed:** maximum COMSOL–Python regional
   difference was 0.00739 °C; RMSE was 0.00135 °C.
4. **Energy and mesh checks passed:** COMSOL mid-hold energy residual was
   0.00996%; the Python medium–fine sensor difference was 0.121 °C.
5. **Experimental accuracy is scoped:** independent spatial RMSE was 1.197 °C
   and fixed-60% regional dynamic RMSE was 1.485 °C.
6. **The strict high-fidelity gate is not rewritten:** V13 maximum error and
   ordering still fail, and V14 variable-PWM validation fails.
7. **Control studies remain gated:** variable-duty ROM, MM-EKF, and decisive
   Oracle MPC are not authorized from this plant.

The honest deliverable is therefore a high-resolution, experimentally
constrained fixed-duty reference, plus a converged COMSOL implementation—not a
general high-fidelity digital twin.  The 55 °C reference is not a safety limit.
