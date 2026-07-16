# V9 COMSOL equivalent-model confirmation

`V9EquivalentConfirmation.java` is an independent COMSOL 6.4 implementation
of the frozen V9 one-dimensional equivalent heat-transfer model.

It is deliberately described as an equivalent finite-element model, not as a
reconstruction of the original paper's CAD, channels, or manifolds. The four
inverse parameters are copied unchanged from
`outputs_v9/frozen_equivalent_parameters.json`.

Run from the repository root:

```powershell
& 'E:\COMSOL\COMSOL64\Multiphysics\bin\win64\comsolcompile.exe' `
  'comsol\v9_confirmation\V9EquivalentConfirmation.java'

& 'E:\COMSOL\COMSOL64\Multiphysics\bin\win64\comsolbatch.exe' `
  -inputfile 'comsol\v9_confirmation\V9EquivalentConfirmation.class' `
  -batchlog 'comsol\v9_confirmation\comsol_confirmation.log'
```

The Java program uses COMSOL's default stationary solver. It writes raw mesh
metrics, profiles, images, and MPH files to this directory. The separate
Python postprocessor checks hashes and frozen values, compares COMSOL against
Python V9, evaluates the unchanged Fig. 7 validation split, and creates the
auditable final deliverables.

Run the postprocessor with a Python environment containing NumPy, pandas, and
Pillow:

```powershell
python 'comsol\v9_confirmation\postprocess_confirmation.py'
```

The required deliverables are written only to `outputs_v9_comsol/` and
`docs_v9_comsol/`. Generated `.class`, `.class.status`, and COMSOL automatic
error-recovery MPH files are build artifacts and are not retained.
