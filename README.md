EdgyEggs is a collection of Houdini's operators providing an early access to some of the popular open source libraries for geometry processing [1],[2],[3],[4],[5],[6]. As name hopefully implies, these nodes meant to be experimental, not polished enough nor crash free enough to become self contained projects yet. Basically the idea is so that a single github account collects some of publicaly available algorthms. Iterested audence can compare them or use in case they can't wait or develop solid, well rounded plugin. 

### Some general considerations
Most of the dependancies shell be provided as git submoudles, unless they are too big to manage them this way. In such case they shell be optional (don't break a build). Intention is to make building and installing process as easy as possible (also in exachange or robustness or flexibility). 

Some of available operators:
### (onging process on wrapping most of libigl functionality)
* IGL Deform - harmonic and As Rigid As Possible mesh deformation
* IGL DiscriteGeometry - compute curvature, laplacian (smooth and sharpen), gradient, eigenvectors on mesh
* IGL UVProject - ARAP and LSCM projections
### (ShapeOp based)
* ShapeOp - optimization of mesh shape under constraints
### (various libraries: SICP/IFGR/CPD)
* Point Cloud Align - implemenetation of SparseICP, Intel Fast Global Registration and Coherent Point Drift algorithms for pointcloud aligment.



## References:
libigl - A simple C++ geometry processing library
[1]: http://libigl.github.io/libigl

		@misc{libigl,
		  title = {{libigl}: A simple {C++} geometry processing library},
		  author = {Alec Jacobson and Daniele Panozzo and others},
		  note = {http://libigl.github.io/libigl/},
		  year = {2016},
		}


ShapeOp - library for static and dynamic geometry processing, using a unified framework for optimization under constraints.
[2]: http://shapeop.org 

		Mario Deuss, Anders Holden Deleuran, Sofien Bouaziz, Bailin Deng, Daniel Piker, Mark Pauly
		ShapeOp - A Robust and Extensible Geometric Modelling Paradigm, Design Modelling Symposium, 2015


El Topo - Robust Topological Operations for Dynamic Explicit Surfaces
[3]: https://www.cs.ubc.ca/labs/imager/tr/2009/eltopo/eltopo.html 

		@article{brochu09,
		author = {Tyson Brochu and Robert Bridson},
		title = {Robust Topological Operations for Dynamic Explicit Surfaces},
		publisher = {SIAM},
		year = {2009},
		journal = {SIAM Journal on Scientific Computing},
		volume = {31},
		number = {4},
		pages = {2472-2493},
		keywords = {interface tracking; dynamic surfaces; triangle meshes; geometric flows},
		url = {http://link.aip.org/link/?SCE/31/2472/1},
		doi = {10.1137/080737617}
		}
CPD - C++ implementation of the Coherent Point Drift point set registration algorithm. http://www.gadom.ski/cpd/
[4]: Myronenko A., Song X. (2010): "Point-Set Registration: Coherent Point Drift", IEEE Trans. on Pattern Analysis and Machine Intelligence, vol. 32, issue 12, pp. 2262-2275, 

Fast Global Registration
[5]:  Qian-Yi Zhou, Jaesik Park, and Vladlen Koltun, ECCV 2016

Sparse Iterative Closest Point (SparseICP)
[6]: C++ implementation for the paper: 

    "Sparse Iterative Closest Point"
    Sofien Bouaziz, Andrea Tagliasacchi, Mark Pauly
    Symposium on Geometry Processing 2013
    Journal: Computer Graphics Forum.
    

