/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "occ.h"
#include "OCCSurface.h"
#include "TopoDS.hxx"
#include "TopoDS_Face.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "BRepMesh_IncrementalMesh.hxx"
#include "Poly_Triangulation.hxx"
#include "BRep_Tool.hxx"
#include "Standard_Version.hxx"

//===========================================================
// occmesh: mesh with TRI on a CAD shape
// use occ mesher based only on hausd
// exports faces and edges
//===========================================================
PyObject* K_OCC::occmesh(PyObject* self, PyObject* args)
{
  PyObject* hook;
  E_Float hausd;
  if (!PYPARSETUPLE_(args, O_ R_, &hook, &hausd)) return NULL;  
  
  GETSHAPE;
  
  // Triangulate shape
  E_Float angularDeflection = 0.5; // en degres
  Standard_Boolean relative = Standard_False;
  if (hausd < 0) { relative = Standard_True; hausd = -hausd; }
  BRepMesh_IncrementalMesh Mesh(*shape, hausd, relative, angularDeflection, Standard_True);
  Mesh.Perform(); 

  // recupere les edges
  PyObject* edges = PyList_New(0);

  // recupere les faces
  PyObject* faces = PyList_New(0);

  E_Int nbNodes = 0;
  E_Int nbTris = 0;
  
  // Dimensionnement
  for (TopExp_Explorer exp(*shape, TopAbs_FACE); exp.More(); exp.Next())
  {
    TopoDS_Face face = TopoDS::Face(exp.Current());

    // get the triangulation
    TopLoc_Location loc;
    Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
    if (tri.IsNull()) continue;
    
#if OCC_VERSION_MAJOR >= 7 && OCC_VERSION_MINOR >= 6
    Handle(TColgp_HArray1OfPnt) nodes = tri->MapNodeArray();
#else
    const TColgp_Array1OfPnt* nodes = &(tri->Nodes());
#endif
    for (Standard_Integer iCount = nodes->Lower(); iCount <= nodes->Upper(); iCount++)
    {
      nbNodes++;  
    }
    nbTris += tri->NbTriangles();
  }
  
  printf("INFO: total number of nodes: " SF_D_ "\n", nbNodes);
  printf("INFO: total number of triangles: " SF_D_ "\n", nbTris);
  
  // buildArray
  PyObject* o = K_ARRAY::buildArray3(3, "x,y,z", nbNodes, nbTris, "TRI", false, 1);
  FldArrayF* f; FldArrayI* c;
  K_ARRAY::getFromArray3(o, f, c);
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);

  E_Int cNodes = 0; E_Int cTris = 0;
  for (TopExp_Explorer exp(*shape, TopAbs_FACE); exp.More(); exp.Next())
  {
    TopoDS_Face face = TopoDS::Face(exp.Current());

    // get the triangulation
    TopLoc_Location loc;
    Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
    if (tri.IsNull()) continue;

    // create a transformation from the location
    gp_Trsf xloc = loc;

#if OCC_VERSION_MAJOR >= 7 && OCC_VERSION_MINOR >= 6
    Handle(TColgp_HArray1OfPnt) nodes = tri->MapNodeArray();
#else
    const TColgp_Array1OfPnt* nodes = &(tri->Nodes());
#endif    
    // get the nodes
    for (Standard_Integer iCount = nodes->Lower(); iCount <= nodes->Upper(); iCount++)
    {
      Standard_Real x,y,z;
      (*nodes)(iCount).Coord(x,y,z);
      //tri->Node(iCount);
      xloc.Transforms(x,y,z);
      fx[cNodes+iCount-1] = x;
      fy[cNodes+iCount-1] = y;
      fz[cNodes+iCount-1] = z;
      //if (iCount > nbNodes) printf("danger iCount>nbNodes %d %d\n", iCount, nbNodes);
      //printf("point %d: %f %f %f\n",iCount,x,y,z);
    }
    
    // copy the polygons
    Standard_Integer i1, i2, i3;
#if OCC_VERSION_MAJOR < 7
    const Poly_Array1OfTriangle& tris = tri->Triangles();
#endif
    for (Standard_Integer iCount = 1; iCount <= tri->NbTriangles(); iCount++) 
    {
      // get the node indexes for this triangle
#if OCC_VERSION_MAJOR < 7
      Poly_Triangle tril = tris(iCount);
#else
      Poly_Triangle tril = tri->Triangle(iCount);
#endif
      tril.Get(i1, i2, i3);
      //c1[stride*(cTris+iCount-1)] = i1+cNodes;
      //c2[stride*(cTris+iCount-1)] = i2+cNodes;
      //c3[stride*(cTris+iCount-1)] = i3+cNodes;
      //if (iCount > nbTris) printf("danger nbTris %d %d\n", iCount, nbNodes);
      //if (i1 > nbNodes) printf("danger1 %d %d\n",i1,nbNodes);
      //if (i2 > nbNodes) printf("danger2 %d %d\n",i2,nbNodes);
      //if (i3 > nbNodes) printf("danger3 %d %d\n",i3,nbNodes);
      //printf("TRI %d: %d %d %d\n", iCount, i1,i2,i3);
    }
    cNodes += nodes->Upper()-nodes->Lower()+1;
    cTris += tri->NbTriangles();
  }

  // export both
  PyObject* tpl = PyList_New(0);
  PyList_Append(tpl, edges); Py_DECREF(edges);
  PyList_Append(tpl, faces); Py_DECREF(faces);
  
  //tpl = K_ARRAY::buildArray3(*coords, "x,y,z,u,v", *cn, -1, "TRI");
  //delete coords; delete cn;

  return tpl;
}
