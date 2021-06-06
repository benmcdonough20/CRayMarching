# CRayMarching
A Ray marching render engine written in c++

The Ray.cpp file is the first iteration of the engine that I wrote in the summer of 2020 (quarantine project!). It implements a recursive ray marching render engine to produce recursive shadows, phong shading with specular, diffuse, and ambient lighting, as well as point and sun lights, a ground plane, a sphere, and a box (shadows still don't work well on the box). The code is purely my own (no copypaste), and I wrote it after reading about rendering fractals with ray marching render engines. The render of the repeated spheres is meant to show the power of the ray marching engine to easily distort space and still produce correct lighting.

The Mandelbulb.cpp file is a tweak I made to submit this project for CPSP 223 at Yale. It was my original goal to render fractals. I added these quickly to showcase the versatility of the renderer. Unfortunately, my original code asked each object to report the surface normal, which was problematic for rendering a fractal. To fix this, I set the ambient color coefficient using an orbital trapping algorithm (which is common for coloring fractals) and I ignored the shadows and other lighting parameters. The mandelbulb distance function is not my own, I borrowed it from this website: http://blog.hvidtfeldts.net/index.php/2011/09/distance-estimated-3d-fractals-v-the-mandelbulb-different-de-approximations/

This project is meant to showcase recursion. The ray marching method calls itself until it registers a hit defined by the minimum distance parameter. In addition, shadows are produced by ray marching in the direction of a bounce determined by the surface normal, which is another example of recursion.

Enjoy!

