# V10 COMSOL literature-channel confirmation

This directory contains an independent COMSOL 6.4 implementation of the
frozen, literature-constrained representative-channel equations. It does not
recover full-cell CAD or a manifold.

From the repository root:

```powershell
& 'E:\COMSOL\COMSOL64\Multiphysics\bin\win64\comsolcompile.exe' `
  'comsol\v10_reference\V10LiteratureChannel.java'

& 'E:\COMSOL\COMSOL64\Multiphysics\bin\win64\comsolbatch.exe' `
  -inputfile 'comsol\v10_reference\V10LiteratureChannel.class' `
  -batchlog 'comsol\v10_reference\comsol_confirmation.log'

.\.venv\Scripts\python.exe `
  'comsol\v10_reference\postprocess_confirmation.py'
```

The Java model uses COMSOL's default stationary solver and runs the six
published cases on 50, 100, and 200 element meshes. The postprocessor compares
the 200-element COMSOL solution with the frozen Python implementation. It does
not fit or alter any parameter.
