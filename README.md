# cascade-shock_buffet Project:
- Aim: Reproduce the experiment result in paper “High‐Speed PIV of shock boundary layer interactions in the transonic buffet flow of a compressor cascade”
- Ref: Andre Weiner/naca0012--shock_buffet 
- Sponsor：Dong
- Acknowledge：Song，[Cheng](https://blog.csdn.net/weixin_39124457/article/details/120152679?spm=1001.2014.3001.5502)


### Running Effeciency 
Parameter: T=0.0001s/foam_divSchemes="SuperBeeV"/NRBCx/yPlus_max=1.0680066/64core each node
| cores |  rhoPimpleFoam | rhoCentralFoam |myFoam| hisa |
|-------------| ---------- |  ---------- |---------- |  ---------- |
| ..| 8:1220 s,16:709 s,32:455 s,128:732 s,256:1223 s |..| ..| ..|
| 64| 395 s stable with small dt=2.0e-9| sensitive to dt  |error| 85 s:most stable dt=7.6e-07,shock wave capture|









