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
// revolve a list of edges

#include "occ.h"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRep_Builder.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "TopoDS.hxx"
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax1.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>

//=====================================================================
// Revolve a list of edges
//=====================================================================
PyObject* K_OCC::revolve(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listEdges;
  E_Float cx, cy, cz;
  E_Float ax, ay, az; 
  E_Float angle; 
  char* name;
  if (!PYPARSETUPLE_(args, OO_ TRRR_ TRRR_ R_ S_, &hook, &listEdges, 
    &cx, &cy, &cz, &ax, &ay, &az, &angle, &name)) return NULL;

  GETSHAPE;
  GETMAPEDGES;

  gp_Pnt center(cx, cy, cz);
  gp_Ax1 axis(center, gp_Dir(ax, ay, az));
  
  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);

  E_Int nedges = PyList_Size(listEdges);
  for (E_Int i = 0; i < nedges; i++)
  {
    PyObject* noO = PyList_GetItem(listEdges, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(no));
    builder.Add(compound, E);
  }
  TopoDS_Shape revolvedSurface = BRepPrimAPI_MakeRevol(compound, axis);

  // rebuild
  TopoDS_Shape* newshp = new TopoDS_Shape(revolvedSurface);
  
#ifdef USEXCAF

  TDocStd_Document* doc = (TDocStd_Document*)packet[5];
  addShape2OCAF(*newshp, name, *doc);
  newshp = copyOCAF2TopShape(*doc);
  delete shape;
  SETSHAPE(newshp);
  Py_INCREF(Py_None);
  return Py_None;

#else

  GETMAPSURFACES;
  
  // Rebuild a single compound
  BRep_Builder builder2;
  TopoDS_Compound compound2;
  builder2.MakeCompound(compound2);
    
  for (E_Int i = 1; i <= surfaces.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(surfaces(i));
    builder2.Add(compound2, F);
  }
  for (E_Int i = 1; i <= edges.Extent(); i++)
  {
    TopoDS_Edge E = TopoDS::Edge(edges(i));
    builder2.Add(compound2, E);
  }

  // Add the revolved faces
  TopTools_IndexedMapOfShape* sfs = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_FACE, *sfs);
  TopTools_IndexedMapOfShape& surfaces2 = *(TopTools_IndexedMapOfShape*)sfs;
  for (E_Int i = 1; i <= surfaces2.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(surfaces2(i));
    builder2.Add(compound2, F);
  }
  delete sfs;

  newshp = new TopoDS_Shape(compound2);

  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after revolve: Nb edges=%d\n", se->Extent());
  printf("INFO: after revolve: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
#endif
}
