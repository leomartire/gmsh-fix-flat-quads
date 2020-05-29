# gmsh-fix-flat-quads
Fixes the flat quadrilaterals in a GMSH .msh file, i.e. the quadrilaterals having jacobian=0 somewhere.

# Usage
```
python GmshFixFlatQuads.py -i extMesh.msh -o extMesh_out.msh
```
```
python GmshFixFlatQuads.py -i extMesh.msh -o extMesh_out.msh -v 1
```

# Illustration
mesh with element 1 having Jacobian = 0 (node 3)             |  modified mesh
:-------------------------:|:-------------------------:
![](https://raw.githubusercontent.com/leomartire/gmsh-fix-flat-quads/master/extMesh.png)  |  ![](https://raw.githubusercontent.com/leomartire/gmsh-fix-flat-quads/master/extMesh_out.png)
