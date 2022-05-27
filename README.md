# Snow Simulation with the Material Point Method

The project implements the material point method for the simulation of snow dynamics, as introduced in the following [paper](https://disney-animation.s3.amazonaws.com/uploads/production/publication_asset/94/asset/SSCTS13_2.pdf). Pioneered as a joint effort between engineers at Disney, and Researchers at UCLA, the paper documents an approach to mimic snow dynamics from a dual grid and particle system, which enables implicit and powerful snow phenomena. This project attempts to replicate the process and explore the material point method as a means of extending computer graphics to handle snow phenomena in excess of simple accumulation.

## Problems and Solutions:

### The Problem with Snow:

Snow is a fundamentally *weird* particle effect. As detailed in the paper, researchers have struggled to centralize "snow solvers", instead building dedicated snow accumulation, or frost generation algorithms, without a catch all general snow modelling algorithm. The authors of the paper recognize that often times, this leads to inefficient, and unrealistic snow simulation, and suggest that primarily, snow is inherently annoying to work with due to its multi-phasing properties. Describing their material point method, the authors suggest, that this MPM has the capacity to deal with inherent tension, collision, and phase adjustment which is endemic to snow simulation. The choice to adopt MPM over pure snow particle simulation was due to the relative inefficiency, and exponentially increasing compute power that explicit particle simulation begets.

As with many things in graphics, this simulation methodology is so distinctly important because snow is simply pretty to look at. Developed as the basis of tools like Matterhorn, the MPM simulator was necessary for Disney's 2015 feature film, *Frozen*. Disney successfully wrote a centralized solver that was able to deal with the intricate minutiae of snow

### The Snow-lution:

This project is a direct implementation (or as close as it can be) to Disney's MPM snow model. The goal is to replicate the entire algorithm (detailed step by step in section 4) almost exactly as Disney does, however, the visualization with be with OpenGL rather than any proprietary technology. The main difficulty lies in digesting the mind-numbing amount of physics in the paper, but beyond that, there will be a necessary research into particle simulation with OpenGL, and the actual pipeline surrounding particle system rendering. As a *very* simplified description of the general solver, the code will generate a series of particles, project them to a grid, solve a system of equations to determine the deformation gradients, new velocities, and new positions, then reload these parameters from the grid to the particles.
