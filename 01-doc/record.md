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

- Exact the point from the picture, by Enqauqe Diaitizer

![image-20220225125819629](/Users/wjq/Library/Application Support/typora-user-images/image-20220225125819629.png)

| name       | upper                                                        | lower                                                        |
| ---------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| blade      | 0.00122,-3e-06<br/>0.01243,0.00973<br/>0.03328,0.01455<br/>0.0712,0.01934<br/>0.15433,0.02767<br/>0.30953,0.04071<br/>0.42558,0.04896<br/>0.54286,0.05721<br/>0.62212,0.05702<br/>0.71834,0.05191<br/>0.81923,0.03825<br/>0.9152,0.02338<br/>0.98557,0.00857<br/>1.00012,0.00488<br/>1.00119,-0.00122 | 0.00122,-3e-06<br/>0.03165,-0.00252<br/>0.10234,-0.00391<br/>0.22902,-0.01031<br/>0.3509,-0.01305<br/>0.49607,-0.01097<br/>0.65109,-0.00525<br/>0.79504,-0.00316<br/>0.90842,-0.00465<br/>0.99137,-0.00364<br/>1.00122,-3e-06 |
| Mas-Ma1.05 | 0.05898,1.4614<br/>0.10861,1.39079<br/>0.15728,1.36316<br/>0.20829,1.31711<br/>0.25799,1.25877<br/>0.30786,1.22807<br/>0.34243,1.08991<br/>0.37595,0.77982<br/>0.41156,0.61404<br/>0.44637,0.51579<br/>0.55459,0.44211<br/>0.60341,0.43904<br/>0.68139,0.4114<br/>0.76152,0.33772<br/>0.84283,0.25482<br/>0.94992,0.19649 | 0.04936,0.08596<br/>0.15027,0.21491<br/>0.25207,0.2886<br/>0.31195,0.29781<br/>0.35097,0.2886<br/>0.45192,0.22412<br/>0.55022,0.12588<br/>0.65008,0.08289<br/>0.75137,0.07368<br/>0.85027,0.07368<br/>0.95039,0.07368 |
| Mas-Ma1.12 | 0.05883,1.43684<br/>0.10735,1.38465<br/>0.1574,1.38158<br/>0.20859,1.36623<br/>0.25823,1.29868<br/>0.30809,1.26491<br/>0.34339,1.24649<br/>0.37932,1.13289<br/>0.41284,0.82281<br/>0.44832,0.63553<br/>0.55247,0.4943<br/>0.60253,0.4943<br/>0.68293,0.4636<br/>0.7631,0.39605<br/>0.84194,0.31009<br/>0.95262,0.23947 | 0.04963,0.12895<br/>0.15056,0.26096<br/>0.25237,0.33772<br/>0.31225,0.34693<br/>0.35123,0.33158<br/>0.45216,0.26404<br/>0.5517,0.16886<br/>0.65032,0.12281<br/>0.75039,0.1136<br/>0.84927,0.11053<br/>0.94941,0.1136 |


