# cascade-shock_buffet Project:
- Aim: Reproduce the experiment result in paper “High‐Speed PIV of shock boundary layer interactions in the transonic buffet flow of a compressor cascade”
- Ref: Andre Weiner/naca0012--shock_buffet 




## Test_cases:

### 01-LowInletVelocity-Acoustic-Duct-Cascade
Solver | InletBoundary |  OutletBoundary | Timescheme   | TurbulenceScheme |  Effeciency |   Info |
|-------------|--------------| -------------|  ---------------|----------------- | ---------- |-------|
|rhoPimpleFoam| Uinlet (1 0 0) | inletOutlet | Backward | IDDES,S-A |  69s,40cores,0.025s,1.0e-6s,97,664 cells  | [01.gif](https://github.com/jiaqiwang969/DLR-buffet/blob/main/Workspace/run/01-DLR-buffet/result/01.gif) |
|rhoPimpleFoam| (1 0 0) sinwave-3000Hz | inletOutlet | Backward | IDDES,S-A | ... | ... |
|rhoPimpleFoam| (1 0 0) sinwave-3000Hz | inletOutlet | Eluer | IDDES,S-A |  ... | ... |
|rhoPimpleFoam| (1 0 0) sinwave-3000Hz | NRBCx | Backward | IDDES,S-A |  ... | ... |
|rhoPimpleFoam| (1 0 0) sinwave-3000Hz | NRBCx | Eluer | IDDES,S-A |  ... | ... |
|LusgsFoam| (1 0 0) sinwave-3000Hz | inletOutlet | Backward | IDDES,S-A | ... | ... |
|LusgsFoam| (1 0 0) sinwave-3000Hz | inletOutlet | Eluer | IDDES,S-A |  ... | ... |
|LusgsFoam| (1 0 0) sinwave-3000Hz | NRBCx | Backward | IDDES,S-A |  ... | ... |
|LusgsFoam| (1 0 0) sinwave-3000Hz | NRBCx | Eluer | IDDES,S-A |  ... | ... |






