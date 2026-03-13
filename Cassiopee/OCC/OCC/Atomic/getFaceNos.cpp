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

#include "Standard_Version.hxx"

//=====================================================================
// Get face numbers from label name
// Return [face no]
//=====================================================================
PyObject* K_OCC::getFaceNos(PyObject* self, PyObject* args)
{
  PyObject* hook;
  char* labelName;
  if (!PYPARSETUPLE_(args, O_ S_, &hook, &labelName)) return NULL;

  GETPACKET;
  GETDOC;  
  GETSHAPETOOL;
  TDF_LabelSequence labels;
  shapeTool->GetShapes(labels);
  
  Handle(TDataStd_Name) name = new TDataStd_Name();

  E_Int istart = 1; E_Int iend = 1;
  E_Bool found = false;
  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    TDF_Label label = labels.Value(i);

    if (shapeTool->IsAssembly(label) == true || shapeTool->IsCompound(label) == true) 
      continue;

    if (label.FindAttribute(TDataStd_Name::GetID(), name)) // retourne tous les attributs de type string
    { 
      TCollection_ExtendedString labelName2 = name->Get();
      if (labelName2 == TCollection_ExtendedString(labelName)) 
      { found = true; }
    }

    TopoDS_Shape shape = shapeTool->GetShape(label);
    TopTools_IndexedMapOfShape sf = TopTools_IndexedMapOfShape();
    TopExp::MapShapes(shape, TopAbs_FACE, sf);
    
    iend = istart + sf.Extent();
    if (found) break;
    istart = iend;
  }

  PyObject* out = PyList_New(0);
  if (found)
  {
    for (E_Int i = istart; i < iend; i++)
    {
      PyObject* o = PyLong_FromLong(i);
      PyList_Append(out, o); Py_DECREF(o);
    }
  }
  return out;
}

//=====================================================================
// Get edge numbers from label name
// Return [edge no]
//=====================================================================
PyObject* K_OCC::getEdgeNos(PyObject* self, PyObject* args)
{
  PyObject* hook;
  char* labelName;
  if (!PYPARSETUPLE_(args, O_ S_, &hook, &labelName)) return NULL;

  GETPACKET;
  GETDOC;
  GETSHAPETOOL;
  TDF_LabelSequence labels;
  shapeTool->GetShapes(labels);
  
  Handle(TDataStd_Name) name = new TDataStd_Name();

  E_Int istart = 1; E_Int iend = 1;
  E_Bool found = false;
  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    TDF_Label label = labels.Value(i);

    if (shapeTool->IsAssembly(label) == true || shapeTool->IsCompound(label) == true) 
      continue;

    if (label.FindAttribute(TDataStd_Name::GetID(), name)) // retourne tous les attributs de type string
    { 
      TCollection_ExtendedString labelName2 = name->Get();
      if (labelName2 == TCollection_ExtendedString(labelName)) 
      { found = true; }
    }

    TopoDS_Shape shape = shapeTool->GetShape(label);
    TopTools_IndexedMapOfShape se = TopTools_IndexedMapOfShape();
    TopExp::MapShapes(shape, TopAbs_EDGE, se);
    
    iend = istart + se.Extent();
    if (found) break;
    istart = iend;
  }

  PyObject* out = PyList_New(0);
  if (found)
  {
    for (E_Int i = istart; i < iend; i++)
    {
      PyObject* o = PyLong_FromLong(i);
      PyList_Append(out, o); Py_DECREF(o);
    }
  }
  return out;
}
