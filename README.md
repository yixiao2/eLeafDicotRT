# Three Dimensional Modelling of Dicot Leaf Photosynthesis
eLeafDicotRT is a end-to-end 3D modelling tool for dicot leaf photosynthesis. Users just need to input a list of anatomical and biochemical parameters, run several lines of code in MATLAB, then eLeafDicotRT will 1) reconstruct a representative 3D leaf anatomy, 2) simulate the light absorption profile with a raytracing algorithm. With the simulated light absorption profiles, light response curves of leaf photosynthesis can be calculated. 

**[Reference]**
Xiao Y, Salesse-Smith C, Driever S, Kromdijk J, Long S. (2023) Understanding the effects of variation in dicot leaf structure on internal light environment and light use efficiency from 3D modeling. []

## Getting Started
In eLeafDicotRT model, 3D reconstruction was achieved through operations in COMSOL Multiphysics driven by algorithms coded in MATLAB. Monte-Carlo ray tracing was implemented in C program to simulate light propagation inside the 3D leaf. All these modules were mastered by several MATLAB scripts.  

### - Prerequisites
Here we used a desktop (Intel i7-9700; 8 Cores; 64GB RAM) with Ubuntu 20.04 operating system.

1. COMSOL Multiphyscis:  
   - COMSOL version >= 5.6
   - Required modules of COMSOL:   
     - CFD module  
     - CAD import module  
     - LiveLink for MATLAB module  

2. MATLAB:  
   - Any compatible MATLAB, we used MATLAB 2020a.  

3. gcc:  
   - We used 9.4.0.  

## Versioning

This is eLeaf code version 0.1b, released on May. 30th 2023. 

## Authors

* code by **Yi Xiao** - *yixiao20@outlook.com* - [https://github.com/xiaoyizz78](https://github.com/xiaoyizz78)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

