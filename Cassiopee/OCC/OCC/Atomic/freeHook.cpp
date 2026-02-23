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
// free CAD hook
#include "occ.h"

#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "TopoDS_Wire.hxx"
#include "BRepCheck_Wire.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "BRep_Builder.hxx"
#include "TDocStd_Document.hxx"

//=====================================================================
// Split shape by max area
//=====================================================================
PyObject* K_OCC::freeHook(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PYPARSETUPLE_(args, O_, &hook)) return NULL;

  GETPACKET;

  // free
  TopTools_IndexedMapOfShape* edges = (TopTools_IndexedMapOfShape*)packet[2];
  TopTools_IndexedMapOfShape* surfaces = (TopTools_IndexedMapOfShape*)packet[1];
  TopoDS_Shape* shape = (TopoDS_Shape*)packet[0];
  char* fileName = (char*)packet[3];
  char* fileFmt = (char*)packet[4];
  
  delete shape;
  delete edges;
  delete surfaces;
  delete [] fileName;
  delete [] fileFmt;
  //delete doc; // can i free it? not sure

  packet[0] = NULL;
  packet[1] = NULL;
  packet[2] = NULL;
  packet[3] = NULL;
  packet[4] = NULL;
  packet[5] = NULL;

  Py_INCREF(Py_None);
  return Py_None;
}
