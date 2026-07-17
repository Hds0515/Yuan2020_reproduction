# COMSOL 6.4 implementation confirmation

`V15Fixed60NetworkConfirmation.java` independently assembles the frozen
thermal equations in COMSOL Global Equations.  It does not import Python
predictions.  The implementation contains 600 surface states and 600 matched
core states and includes:

- in-plane and through-stack conduction;
- nonuniform fan flow with plenum mixing;
- conservative air-enthalpy marching through every cathode path;
- linearized radiation plus edge convection;
- two thermal masses with surface–core exchange;
- thermoneutral heat generation with the frozen spatial heat pattern.

COMSOL 6.4.0.293 compiled the Java source and solved the 4 A ramp with its BDF
transient solver.  The solve converged in 7.09 s of reported core solver time.
Against the independently executed Python model on the same 1200-state mesh,
the maximum regional-temperature difference was 0.00739 °C and the RMSE was
0.00135 °C.  Both are below the frozen 0.05 °C threshold.

The energy ledger evaluates heat input, heat rejection, and stored-energy rate
as separate expressions.  Exact 100 s switching instants have two legitimate
one-sided state derivatives and are excluded from residual scoring.  At the
mid-hold experimental sample times, the maximum residual was 0.00996%, below
the 0.5% limit.

Reproduction commands:

```powershell
& 'E:\COMSOL\COMSOL64\Multiphysics\bin\win64\comsolcompile.exe' `
  comsol\v15_reference\V15Fixed60NetworkConfirmation.java

& 'E:\COMSOL\COMSOL64\Multiphysics\bin\win64\comsolbatch.exe' `
  -inputfile comsol\v15_reference\V15Fixed60NetworkConfirmation.class `
  -batchlog comsol\v15_reference\v15_run.log -batchlogout
```

The `.mph` file is a discrete energy-network implementation for solver
confirmation.  It is not an explicit fan, plenum, manifold, membrane-water, or
electrochemical CAD reconstruction.
