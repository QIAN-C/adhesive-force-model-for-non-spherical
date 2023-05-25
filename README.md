# a-adhesive-force-model-for-non-spherical 
# General information
adhesive-force-model-for-non-spherical is a computational program for modeling the adhesive contact forces between non-spherical particles. This program operates based on the open-source software Mechsys, so you need to install Mechsys completely before running it.
# Installation
The library is developed based on Mechsys. Therefore, a complete installation of Mechsys is required. You can install it using the following link: https://mechsys.nongnu.org/index.html; or https://github.com/topics/mechsys;
# Usage
Once you have downloaded Mechsys software completely according to the official website instructions, you need to replace the "distance", "interacton" and "particle" file in the folder with the one from "mechsys\lib\dem". 
After that, you can compile it using "ccmake ." and then run it with the command "make" by using ".singlefiber".
# Attention

The program in the "smooth surface" folder corresponds to the computational program for the bouncing process of cubic particles on a smooth surface. This program includes a contact force model for non-spherical particles. Any adhesive force calculation program can be extended based on this.
