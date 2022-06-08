# cascade-shock_buffet Project:
- Aim: Reproduce the experiment result in paper “High‐Speed PIV of shock boundary layer interactions in the transonic buffet flow of a compressor cascade”
- Ref: 1. Andre Weiner/naca0012--shock_buffet , 2. [turbulence model](https://github.com/jiaqiwang969/SSTtransition-turbulence-model) , 3. [AeroAcoustics solver](https://github.com/jiaqiwang969/Hybrid-Methods-in-Openfoam) 4. [NSCBC](https://github.com/jiaqiwang969/NSCBC-openfoam) 5. [dbns](https://github.com/ilyapopov/dbns-port/blob/259e46a7b4f96301f872498f950702aabef4ed92/dbnsFoam/dbnsFoam.C) 6. thesis-Linear aeroacoustic solver in OpenFOAM 7. [vertification with cases-axibump](https://github.com/jiaqiwang969/case_axibump)
- Sponsor：Dong
- Acknowledge：Song Moru, Cheng Long




### Frame：

<img src="https://cdn.mathpix.com/snip/images/FLoKLu0FjsDo5S0EtAbn2p4PrqaI1Vb70ONrhrDzlPY.original.fullsize.png" />

### Running Effeciency 
Parameter: T=0.0001s/foam_divSchemes="SuperBeeV"/NRBCx/yPlus_max=1.0680066/64core each node
| cores |  rhoPimpleFoam | rhoCentralFoam |myFoam| hisa |
|-------------| ---------- |  ---------- |---------- |  ---------- |
| ..| 8:1220 s,16:709 s,32:455 s,128:732 s,256:1223 s |..| ..| ..|
| 64| 395 s stable with small dt=2.0e-9| sensitive to dt  |error| 85 s:most stable dt=7.6e-07,shock wave capture|


### Exploring issues

- Instability

solution1: choose better limiter to alleviate dispersion effect
<img src="https://cdn.mathpix.com/snip/images/cSA0LK_DRCgxbaZwfoskt4tkcWjpIwyjpcwmE_sjL5E.original.fullsize.png" width="640px">

- testcase02: also instability in the shear layer （Is it a numerically similar instability with increased grid and denser stripes？）
<img src="https://cdn.mathpix.com/snip/images/jlel6iby3EftNmMgQbaRwX_sasyTB8jxKytpoEKPbVk.original.fullsize.png" width="340px">
<img src="https://cdn.mathpix.com/snip/images/wvuHTNBzri8Q1RcADlwxndsthXQqjBNUUjwM4MoyrGQ.original.fullsize.png" width="340px">

<img src="https://cdn.mathpix.com/snip/images/fplVunmCng73ChZb_SO8JAAOLQ6xp1wboAaoEsbUUYk.original.fullsize.png" width="340px">

- same effect also find in tailing region
<img src="https://cdn.mathpix.com/snip/images/5-CmkH4Jk3X_ZMMmvwS7JOLL7lc7mKp4j3UK24TsKVc.original.fullsize.png" width="240px">


### Gallery
- Performace
<div align="left">
  <img src="https://github.com/jiaqiwang969/DLR-buffet/blob/main/Performance.jpeg" width="400px">
</div>


- [Naca65](https://www.cfd-china.com/topic/2591/rhosimplefoam-%E6%B1%82%E8%A7%A3%E5%8F%AF%E5%8E%8B%E7%BC%A9%E6%B5%81%E5%8A%A8%E6%97%B6%E4%B8%80%E4%B8%AA%E5%A5%87%E6%80%AA%E7%9A%84%E7%8E%B0%E8%B1%A1/32?_=1654024267410)
<img src="https://cdn.mathpix.com/snip/images/xj3BPyQ0S0EsbJK94W4sptt8LC-qpDqti9PCr04qbyM.original.fullsize.png" width="640px">

- DLR-compressor: better to enlarge the upstream channel
<img src="https://cdn.mathpix.com/snip/images/EwGO2iP8OZrG_aLwgY2W7RzI1mLizkCCj2_9UJTCk_8.original.fullsize.png" width="640px">
<img src="https://cdn.mathpix.com/snip/images/YwxuxBPoW6K1iEF1JKlGbKzZO_CtuimX5UTBu-niD0Y.original.fullsize.png" width="640px">
<img src="https://cdn.mathpix.com/snip/images/vDFwkqAIId0aSXhd9ji53H8UPbIBQkZkv60RMrlrlKk.original.fullsize.png" width="640px">


- [TTT02-hisa-1.05-0.71-AUSMPlusUp-5000-maxCo1-5.0e-8-long2-avi.avi](https://www.youtube.com/watch?v=ASdCueqCFpY)
