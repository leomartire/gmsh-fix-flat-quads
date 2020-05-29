# gmsh-fix-flat-quads
Fixes the flat quadrilaterals in a GMSH .msh file, i.e. the quadrilaterals having jacobian=0 somewhere.

# usage
python GmshFixFlatQuads.py -i extMesh.msh -o extMesh_out.msh

python GmshFixFlatQuads.py -i extMesh.msh -o extMesh_out.msh -v 1
