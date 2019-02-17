EdgyEggs is a collection of Houdini's operators providing an early access to some of the popular open source libraries for geometry processing. As name hopefully implies, these nodes are meant to be experimental. Basically the idea is so that a single github account collects some of publicly available algorithms. Interested audience can compare them or use in case they can't wait for or develop own solid well rounded plugin. 

Most of the dependancies we provid as git submodules, unless they are too big to manage them this way. In such case they shell be optional.

### IGL
* IGL Deform - harmonic and As Rigid As Possible mesh deformation (SOP)
* IGL Discrete Geometry - compute curvature, laplacian (to smooth and sharpen mesh), gradient, eigenvectors (SOP)
* IGL UV Project - ARAP and LSCM projections (SOP)
### ShapeOp 
* ShapeOp - optimization of mesh under constraints (SOP)
### Various libraries for point-cloud alignment (SparseICP/IntelFGR/CPD)
* Point Cloud Align - implementation of SparseICP, Intel Fast Global Registration and Coherent Point Drift algorithms for point-cloud alignment (SOP)



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
		ShapeOp - A Robust and Extensible Geometric Modeling Paradigm, Design Modeling Symposium, 2015


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
    
