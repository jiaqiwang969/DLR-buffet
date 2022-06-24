# cascade-shock_buffet Project: STEP BY STEP
<img src="https://cdn.mathpix.com/snip/images/vDFwkqAIId0aSXhd9ji53H8UPbIBQkZkv60RMrlrlKk.original.fullsize.png" width="340px">

- Aim: Reproduce the experiment result in paper “High‐Speed PIV of shock boundary layer interactions in the transonic buffet flow of a compressor cascade”
- Ref: 1. [Andre Weiner/naca0012--shock_buffet](https://github.com/AndreWeiner/naca0012_shock_buffet) , 2. [turbulence model](https://github.com/jiaqiwang969/SSTtransition-turbulence-model) , 3. [AeroAcoustics solver](https://github.com/jiaqiwang969/Hybrid-Methods-in-Openfoam) 4. [NSCBC](https://github.com/jiaqiwang969/NSCBC-openfoam) 5. [dbns](https://github.com/ilyapopov/dbns-port/blob/259e46a7b4f96301f872498f950702aabef4ed92/dbnsFoam/dbnsFoam.C) 6. thesis-Linear aeroacoustic solver in OpenFOAM 7. [Axis-2Dbump](https://github.com/jiaqiwang969/Axis-2Dbump) 8. [WENO](https://github.com/jiaqiwang969/WENOEXT-project) 9.[RK4+LTS](https://github.com/jiaqiwang969/rhoCentralRKFoam) 10. [WENO-example](https://figshare.com/authors/Panagiotis_Tsoutsanis/2595451) 11. [thermophysical](https://cfd.direct/openfoam/user-guide/v6-thermophysical/) 12. Markus Zauner-phd-thesis
- Sponsor：Dong
- Acknowledge：Song Moru, Cheng Long，Du Lei


# Gallery

### [hisa solver](https://gitlab.com/hisa/hisa), with SAMU+plus+NSCBC (1 order space scheme)
```
    reconstruct(rho) wVanLeer;
    reconstruct(U)   wVanLeer;
    reconstruct(T)   wVanLeer;
```


Animation: [TTT02-hisa-1.05-0.71-AUSMPlusUp-5000-maxCo1-5.0e-8-long2-avi.avi](https://www.youtube.com/watch?v=ASdCueqCFpY)

- Cycle

<img src="https://cdn.mathpix.com/snip/images/CtJf8WQpi7ubGOaAumKUYiIy4vF2L4pP7iee_H8AjMA.original.fullsize.png" width="340px">

Simu_freq_buffet $\approx$ 118Hz, T_cycle= 0.00845 s;    (Experiments_freq=171Hz)

### [blastFoam solver](https://github.com/synthetik-technologies/blastfoam), with SAMU+plus+waveTransimision  (3 order space scheme)

```
    reconstruct(rho)               quadraticMUSCL Minmod;
    reconstruct(U)                 quadraticMUSCL Minmod;
    reconstruct(e)                 quadraticMUSCL Minmod;
    reconstruct(p)                 quadraticMUSCL Minmod;
    reconstruct(speedOfSound)      quadraticMUSCL Minmod;
```

- Early stages
<img src="https://cdn.mathpix.com/snip/images/AedspmaItmFj-sckjw2DFNKvCnb0IgZiZlBmCx9V6X8.original.fullsize.png" width="340px">

- Interaction
<img src="https://cdn.mathpix.com/snip/images/d5NlV0GG0Nd01WX2jLtbdtbfGj5wLIVIt4MOBcTqN_c.original.fullsize.png" width="340px">

- Cycle

waiting...



### Parameter

<img src="https://cdn.mathpix.com/snip/images/KSFuMM79-hnuki6cm_LVKSspXaCo1KNqRHBsfirzhMw.original.fullsize.png" width="340px">

```
mu              =7.1488e-07;  %kinematic viscosity  [m^2/s]
nu=mu/rho       =0.00000122;  %dynamic viscosity  [kg/(m*s)]
dy              =1.42e-7；
gamma           =1.4;
R               =287.05;
ptotal          =101325;
Mach            =1.2;
Tinlet          =300;
soundSpeedInlet =sqrt(gamma*R*Tinlet);
vinlet          =soundSpeedInlet*Mach;  % 416.6627m/s
pinlet          =ptotal/(1+0.5*(vinlet)^2/R/Tinlet); %5.0461e+04 pa
rhoinlet        =pinlet/R/Tinlet;       % 0.5860 kg/m^3
Ttotal          =Tinlet + Tinlet*Mach*Mach*(gamma-1)/2;
```






### Running Effeciency 
Parameter: T=0.0001s/foam_divSchemes="SuperBeeV"/NRBCx/yPlus_max=1.0680066/64core each node
| cores |  rhoPimpleFoam | rhoCentralFoam |myFoam| hisa |
|-------------| ---------- |  ---------- |---------- |  ---------- |
| ..| 8:1220 s,16:709 s,32:455 s,128:732 s,256:1223 s |..| ..| ..|
| 64| 395 s stable with small dt=2.0e-9| sensitive to dt  |error| 85 s:most stable dt=7.6e-07,shock wave capture|


### Exploring issues

- Instability

solution1: choose better limiter to alleviate dispersion effect

<img src="https://cdn.mathpix.com/snip/images/cSA0LK_DRCgxbaZwfoskt4tkcWjpIwyjpcwmE_sjL5E.original.fullsize.png" width="340px">

- testcase02: also instability in the shear layer （Is it a numerically similar instability with increased grid and denser stripes？）
<img src="https://cdn.mathpix.com/snip/images/jlel6iby3EftNmMgQbaRwX_sasyTB8jxKytpoEKPbVk.original.fullsize.png" width="340px">
<img src="https://cdn.mathpix.com/snip/images/wvuHTNBzri8Q1RcADlwxndsthXQqjBNUUjwM4MoyrGQ.original.fullsize.png" width="340px">

<img src="https://cdn.mathpix.com/snip/images/fplVunmCng73ChZb_SO8JAAOLQ6xp1wboAaoEsbUUYk.original.fullsize.png" width="340px">

- same effect also find in tailing region
<img src="https://cdn.mathpix.com/snip/images/5-CmkH4Jk3X_ZMMmvwS7JOLL7lc7mKp4j3UK24TsKVc.original.fullsize.png" width="240px">



## First Results: Oscillation！！

<!-- <img src="https://cdn.mathpix.com/snip/images/jIiNlXnmTouWWUWyQvveZ0r6STzL6axQXNV08WRW7gA.original.fullsize.png" width="640px">
 -->
<!-- - Steady vribraion in one cycle？？ -->
<img src="https://cdn.mathpix.com/snip/images/3TwqsU5mEdTTtdI-rreur58Cx0-ILZJjI8zB1vMuoRE.original.fullsize.png"  width="340px">
<img src="https://cdn.mathpix.com/snip/images/MjcNL-z5TpqnZEo8hJdfUb9hSYwATnmwUL5DBhhH2yc.original.fullsize.png" width="340px">
<img src="https://cdn.mathpix.com/snip/images/l6w45yANJMQ3WLqx4kIweDZB3Y049eWN9ILyfqo_v1k.original.fullsize.png" width="340px">

<!-- - Performace: not bad not good~
<div align="left">
  <img src="https://github.com/jiaqiwang969/DLR-buffet/blob/main/Performance.jpeg" width="400px">
</div>
 -->

<!-- - [Naca65](https://www.cfd-china.com/topic/2591/rhosimplefoam-%E6%B1%82%E8%A7%A3%E5%8F%AF%E5%8E%8B%E7%BC%A9%E6%B5%81%E5%8A%A8%E6%97%B6%E4%B8%80%E4%B8%AA%E5%A5%87%E6%80%AA%E7%9A%84%E7%8E%B0%E8%B1%A1/32?_=1654024267410)
<img src="https://cdn.mathpix.com/snip/images/xj3BPyQ0S0EsbJK94W4sptt8LC-qpDqti9PCr04qbyM.original.fullsize.png" width="640px">

- DLR-compressor: better to enlarge the upstream channel
<img src="https://cdn.mathpix.com/snip/images/EwGO2iP8OZrG_aLwgY2W7RzI1mLizkCCj2_9UJTCk_8.original.fullsize.png" width="640px">

- Enlarge the upstream
<img src="https://cdn.mathpix.com/snip/images/YwxuxBPoW6K1iEF1JKlGbKzZO_CtuimX5UTBu-niD0Y.original.fullsize.png" width="640px">

 -->

# Improv01: mesh
After learning the Lessons from [axis-2d bump](https://github.com/jiaqiwang969/Axis-2Dbump), we start to improve this case!!

### Remesh: the mesh is a natural Filter

<img src="https://cdn.mathpix.com/snip/images/nYrqq5zrX3B_iPKwCAuQ18KIPIXD_dwy8rS8WeY_6q0.original.fullsize.png" width="640px">

- 3D mesh with extrudeMesh & createPatch command
<img src="https://cdn.mathpix.com/snip/images/3oKYYcGNBwa8LSVYYePmkw1F1RpZOBajSgBnbTeGZEc.original.fullsize.png"  width="340px">



## Improve02: add [WENO](https://github.com/jiaqiwang969/WENOEXT-project)+[RK4+LTS](https://github.com/jiaqiwang969/rhoCentralRKFoam)+[NSCBC](https://github.com/jiaqiwang969/NSCBC-openfoam)
<img src="https://cdn.mathpix.com/snip/images/_Oi4PRJVOW4nogHdeQB8erf3A44_ALxHGeLf_cIRA98.original.fullsize.png" width="640px">

## Improve03: add [acousticdampingsource](https://caefn.com/openfoam/fvoptions-acousticdampingsource)
<img src="https://cdn.mathpix.com/snip/images/z72z1Vg0FoZNlVHsQcJPFqQ4G5oBzpFLQxoCLRAZ9os.original.fullsize.png" width="640px">
“acousticDampingSource” to damp acoustic waves generated from unsteady flow before they propagate to imperfect non-reflective inflow/outflow boundaries!

