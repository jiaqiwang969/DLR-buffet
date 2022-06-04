# cascade-shock_buffet Project:
- Aim: Reproduce the experiment result in paper “High‐Speed PIV of shock boundary layer interactions in the transonic buffet flow of a compressor cascade”
- Ref: Andre Weiner/naca0012--shock_buffet 
- Sponsor：Dong
- Acknowledge：Song Moru




### Frame：

<img src="https://cdn.mathpix.com/snip/images/FLoKLu0FjsDo5S0EtAbn2p4PrqaI1Vb70ONrhrDzlPY.original.fullsize.png" />

### Running Effeciency 
Parameter: T=0.0001s/foam_divSchemes="SuperBeeV"/NRBCx/yPlus_max=1.0680066/64core each node
| cores |  rhoPimpleFoam | rhoCentralFoam |myFoam| hisa |
|-------------| ---------- |  ---------- |---------- |  ---------- |
| ..| 8:1220 s,16:709 s,32:455 s,128:732 s,256:1223 s |..| ..| ..|
| 64| 395 s stable with small dt=2.0e-9| sensitive to dt  |error| 85 s:most stable dt=7.6e-07,shock wave capture|


### Result
<div align="left">
  <img src="https://github.com/jiaqiwang969/DLR-buffet/blob/main/Performance.jpeg" width="400px">
</div>

- Further, deep learning will be implemented to approximate the F_length curve. This will open a new project called [TransitionAI](https://github.com/jiaqiwang969/TransitionAI).



### Gallery

- [Naca65](https://www.cfd-china.com/topic/2591/rhosimplefoam-%E6%B1%82%E8%A7%A3%E5%8F%AF%E5%8E%8B%E7%BC%A9%E6%B5%81%E5%8A%A8%E6%97%B6%E4%B8%80%E4%B8%AA%E5%A5%87%E6%80%AA%E7%9A%84%E7%8E%B0%E8%B1%A1/32?_=1654024267410)
<img src="https://cdn.mathpix.com/snip/images/xj3BPyQ0S0EsbJK94W4sptt8LC-qpDqti9PCr04qbyM.original.fullsize.png" />

- DLR-compressor
<img src="https://cdn.mathpix.com/snip/images/EwGO2iP8OZrG_aLwgY2W7RzI1mLizkCCj2_9UJTCk_8.original.fullsize.png" />
