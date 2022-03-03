# cascade-shock_buffet Project:
- Aim: Reproduce the experiment result in paper “High‐Speed PIV of shock boundary layer interactions in the transonic buffet flow of a compressor cascade”
- Ref: Andre Weiner/naca0012--shock_buffet 
- 致谢：董广明，宋沫儒，[陈十七](https://blog.csdn.net/weixin_39124457/article/details/120152679?spm=1001.2014.3001.5502)



## Test_cases:

### 01-LowInletVelocity-Acoustic-Duct-Cascade
Solver | InletBoundary |  OutletBoundary | Timescheme   | TurbulenceScheme |divSchemes|  Effeciency |   Info |
|-------------|--------------| -------------|  ---------------|----------------- | -------- | ---------- |-------|
|rhoPimpleFoam| Uinlet (1 0 0) | inletOutlet | Backward | IDDES,S-A |  div(phi,U):Gauss MinmodV| 69s,40cores,0.025s,1.0e-6s,97,664 cells  | [01.gif](https://github.com/jiaqiwang969/DLR-buffet/blob/main/Workspace/run/01-DLR-buffet/result/01.gif) |
|rhoPimpleFoam| (1 0 0) sinwave-3000Hz | inletOutlet | Backward | IDDES,S-A | div(phi,U):Gauss MinmodV| ... | ... |
|rhoPimpleFoam| (1 0 0) sinwave-3000Hz | inletOutlet | Eluer | IDDES,S-A |div(phi,U):Gauss MinmodV|   ... | ... |
|rhoPimpleFoam| (1 0 0) sinwave-3000Hz | NRBCx | Backward | IDDES,S-A |div(phi,U):Gauss MinmodV|   ... | ... |
|rhoPimpleFoam| (1 0 0) sinwave-3000Hz | NRBCx | Eluer | IDDES,S-A | div(phi,U):Gauss MinmodV|  ... | ... |
|LusgsFoam| (1 0 0) sinwave-3000Hz | inletOutlet | Backward | IDDES,S-A | div(phi,U):Gauss MinmodV| ... | ... |
|LusgsFoam| (1 0 0) sinwave-3000Hz | inletOutlet | Eluer | IDDES,S-A |div(phi,U):Gauss MinmodV|   ... | ... |
|LusgsFoam| (1 0 0) sinwave-3000Hz | NRBCx | Backward | IDDES,S-A |div(phi,U):Gauss MinmodV|   ... | ... |
|LusgsFoam| (1 0 0) sinwave-3000Hz | NRBCx | Eluer | IDDES,S-A |div(phi,U):Gauss MinmodV|   ... | ... |


### 02-HighInletVelocity-NRBCxOutlet-Acoustic-Duct-Cascade-1.2-Mach
CaseName |  Effeciency |   Info |
|-------------| ---------- |-------|
|04-test02-rhoPimpleFoam-time1e9-Uinlet12-inletOutlet-backward-IDDES-MinmodV| | |
|04-test02-rhoPimpleFoam-time1e9-Uinlet12-inletOutlet-backward-IDDES-SuperBee| | |
|04-test02-rhoPimpleFoam-time1e9-Uinlet12-inletOutlet-backward-IDDES-vanLeer| | |

NEXT LIST

- Add surface force,y+,... postprocessing
- Add rhoCentralFoam, LusgsFoam for comparision
- Plan: how to compare with experiment data









