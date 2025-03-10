// geo file for meshing with Gmsh meshing software created by FreeCAD

// enable multi-core processing
General.NumThreads = 12;

// open brep geometry
Merge "beso_031_solid001_Geometry.brep";

// group data
Physical Surface("ConstraintFixed") = {155, 1703, 1783, 1851, 1931, 227, 306, 385, 458, 541, 613, 685, 762, 77};
Physical Surface("ConstraintForce") = {859};

// Characteristic Length
// no boundary layer settings for this mesh
// min, max Characteristic Length
Mesh.CharacteristicLengthMax = 1e+22;
Mesh.CharacteristicLengthMin = 0.0;
Mesh.MeshSizeFromCurvature = 12; // number of elements per 2*pi radians, 0 to deactivate

// optimize the mesh
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 0;
// High-order meshes optimization (0=none, 1=optimization, 2=elastic+optimization, 3=elastic, 4=fast curving)
Mesh.HighOrderOptimize = 0;

// mesh order
Mesh.ElementOrder = 1;

// mesh algorithm, only a few algorithms are usable with 3D boundary layer generation
// 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad, 9=Packing Parallelograms, 11=Quasi-structured Quad)
Mesh.Algorithm = 2;
// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 7=MMG3D, 9=R-tree, 10=HTX)
Mesh.Algorithm3D = 1;

// subdivision algorithm
Mesh.SubdivisionAlgorithm = 0;

// incomplete second order elements
Mesh.SecondOrderIncomplete = 0;

// meshing
Geometry.Tolerance = 1e-06; // set geometrical tolerance (also used for merging nodes)
Mesh  3;
Coherence Mesh; // Remove duplicate vertices

// save
// For each group save not only the elements but the nodes too.;
Mesh.SaveGroupsOfNodes = 1;
// Needed for Group meshing too, because for one material there is no group defined;
// Ignore Physical definitions and save all elements;
Mesh.SaveAll = 1;
Save "beso_031_solid001_Mesh.unv";


// **********************************************************************
// Gmsh documentation:
// https://gmsh.info/doc/texinfo/gmsh.html#Mesh
//
// We do not check if something went wrong, like negative jacobians etc. You can run Gmsh manually yourself: 
//
// to see full Gmsh log, run in bash:
// /Applications/Gmsh.app/Contents/MacOS/gmsh - /Users/yimingliu/Desktop/CellPackage/Result/beso_30/FEMMeshGmsh/shape2mesh.geo
//
// to run Gmsh and keep file in Gmsh GUI (with log), run in bash:
// /Applications/Gmsh.app/Contents/MacOS/gmsh /Users/yimingliu/Desktop/CellPackage/Result/beso_30/FEMMeshGmsh/shape2mesh.geo
