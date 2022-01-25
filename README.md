# VAWTAero.jl

This repository contains a set of aerodynamic tools for VAWTs both steady and unsteady operation, 2D and
3D (stacked 2D) convenience functions along with coupling to NREL's InflowWind for turbulent inflow. You
will need to provide your own .bts files from turbsim (compile or download OpenFAST and use the turbsim
binary with an .inp file.  There is an example in the test/data/ifw folder)

Double Multiple Streamtube implementation per https://doi.org/10.5194/wes-2019-44

Actuator Cylinder Implementation from https://github.com/byuflowlab/vawt-ac (updated and modified)

CACTUS Dynamic Stall models https://github.com/SNL-WaterPower/CACTUS

3D VAWT error resolution and unsteady method numerical acceleration (RPI) per https://arc.aiaa.org/doi/abs/10.2514/1.J060476
