# GraAda3D

Developer: Yiyuan Zhong 钟乙源

zhongyiyuan123@gmail.com

## 1 Introduction

It is a 3D gravity inversion program implemented with C++. The inversion domain is discretized as rectangular prismatic meshes in the Cartesian coordinate system. The inversion mesh can be adaptively refined to boost computational performance and avoid over-parameterization. In addition, users can choose to use L0-norm (minimum support) and L1-norm regularization to obtain a focused image, or use L2-norm regularization to get a smooth result. Furthermore, a-priori information can be incorporated in inversion through cross-gradient coupling or direct parameter relationship. 

Any combination of gravity field components or gravity gradient components (gz, gx, gy, Tzz, Txz, Tyz, Txx, Txy, Tyy) can be used as input data.  Exact analytical solutions of gravity field and gravity gradient tensor are used to ensure accuracy of the forward modeling. 

**Hightlights**

- Adaptively refined mesh for inversion parameters
- Incorporation of  a-priori constraints through cross-gradient coupling or direct parameter relation
- Lp-norm regularization (p=0,1,2,...)

## 2 How to use it

I am reworking on the documentation. It will be available soon.



