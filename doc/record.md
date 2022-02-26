## 06-cascade-shock_buffet

- Aim: Reproduce the experiment result in paper “High‐Speed PIV of shock boundary layer interactions in the transonic buffet flow of a compressor cascade”
- Ref: Andre Weiner/naca0012--shock_buffet

#### Step1：Add experimental data/geometry file for comparison.

- [pressure coeffient](https://en.wikipedia.org/wiki/Pressure_coefficient)

$$
C_{p}=\frac{p-p_{\infty}}{\frac{1}{2} \rho_{\infty} V_{\infty}^{2}}=\frac{p-p_{\infty}}{p_{0}-p_{\infty}}
$$

where:

![p](https://wikimedia.org/api/rest_v1/media/math/render/svg/81eac1e205430d1f40810df36a0edffdc367af36) is the [static pressure](https://en.wikipedia.org/wiki/Static_pressure#Static_pressure_in_fluid_dynamics) at the point at which pressure coefficient is being evaluated

![p_{\infty }](https://wikimedia.org/api/rest_v1/media/math/render/svg/d4b1d64d7d5bc5f0a8490d9844565f5cb2b83802) is the static pressure in the [freestream](https://en.wikipedia.org/wiki/Freestream) (i.e. remote from any disturbance)

![p_{0}](https://wikimedia.org/api/rest_v1/media/math/render/svg/2b969ada68a88e2aeba9a2d2096abaf1fd53c21d) is the [stagnation pressure](https://en.wikipedia.org/wiki/Stagnation_pressure)/total pressure in the [freestream](https://en.wikipedia.org/wiki/Freestream) (i.e. remote from any disturbance)

![\rho _{\infty }](https://wikimedia.org/api/rest_v1/media/math/render/svg/496bc807b3a0e1d4dc492b892944f61159aab11e) is the freestream [fluid density](https://en.wikipedia.org/wiki/Density) (Air at [sea level](https://en.wikipedia.org/wiki/Sea_level) and 15 °C is 1.225 ( ![{\displaystyle {\rm {kg/m^{3}}}}](https://wikimedia.org/api/rest_v1/media/math/render/svg/0e6a2bd14730d35d6495bc5d5f95ae7d7f1ee3d2))

![V_{\infty }](https://wikimedia.org/api/rest_v1/media/math/render/svg/0b9ee987969cd445eb595eb4b44367b42a8aeeb7) is the freestream velocity of the fluid, or the velocity of the body through the fluid



- [Isentropic Mach number](https://www.researchgate.net/publication/224993164_Numerical_investigation_of_the_effect_of_intake_integration_at_transonic_speeds_on_the_DLR-F17E_UCAV_model_for_the_TWG)

The isentropic Mach number is the Mach number you would have without any losses in the flow. This property is often used ==to investigate the ideal surface Mach number you would have without losses and walls without any friction== (slip conditions). The isentropic Mach number is for example often plotted on turbomachinery blades. The isentropic Mach number can be computed from the [[isentropic-flow-relations|thermo.isentropic-compression]] using the following formula: 
$$
\frac{p_0}{p} = \left( 1+\frac{\gamma -1}{2}\operatorname{Ma_s}^2\right)^{\frac{\gamma}{\gamma -1}}
$$


This gives the following formula for the isentropic Mach number: 
$$
\operatorname{Ma_s}=\sqrt{\frac{2}{\gamma -1}\left[\left(\frac{p_0}{p}\right)^{\frac{\gamma -1 }{\gamma}}-1\right]}
$$
 where $p_\infty$ is the total pressure in the freestream outside of the boundary layers. Usually the inlet total pressure is used, $p$ is the local static pressure, and $\gamma$ is the ratio of specific heat for air.

#### Step2：Make three meshes by numeca autogrid (80000/130000/300000)



