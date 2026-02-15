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

#include "TDF_Label.hxx"
#include "TDF_LabelSequence.hxx"
#include "TDF_Tool.hxx"

#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFDoc_DocumentTool.hxx"
#include "XCAFDoc_ShapeMapTool.hxx"

#include "TDocStd_Document.hxx"
#include "TDataStd_Name.hxx"

#include "TDF_ChildIterator.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Face.hxx"
#include "Geom_Surface.hxx"
#include "BRep_Tool.hxx"

//=====================================================================
// set face name label to faces
//=====================================================================
PyObject* K_OCC::setFaceNameInOCAF(PyObject* self, PyObject* args)
{
  PyObject* hook;
  char* name;
  PyObject* listFaces;
  if (!PYPARSETUPLE_(args, OO_ S_, &hook, &listFaces, &name)) return NULL;

  GETPACKET;
  GETMAPSURFACES;

  TDocStd_Document* doc = (TDocStd_Document*)packet[5];
  if (doc == NULL) 
  {
    PyErr_SetString(PyExc_TypeError, "setFaceNameInOCAF: no OCAF document.");
    return NULL;
  }

  Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());

  // create label
  //TDF_Label label = doc->Main().NewChild();

  // set label to faces
  E_Int nfaces = PyList_Size(listFaces);
  for (E_Int no = 0; no < nfaces; no++)
  {
    PyObject* noFaceO = PyList_GetItem(listFaces, no);
    E_Int noFace = PyInt_AsLong(noFaceO);
    const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
    TDF_Label faceLabel = shapeTool->AddShape(F);
    TDataStd_Name::Set(faceLabel, name);
  }

  Py_INCREF(Py_None);
  return Py_None;
}