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

//=====================================================================
// Change label name in ocaf
//=====================================================================
PyObject* K_OCC::changeLabelNameInOCAF(PyObject* self, PyObject* args)
{
  PyObject* hook;
  char* oldName; char* newName;
  if (!PYPARSETUPLE_(args, O_ SS_, &hook, &oldName, &newName)) return NULL;

  GETPACKET;
  
  TDocStd_Document* doc = (TDocStd_Document*)packet[5];
  
  // Get labels corresponding to shapes
  TDF_LabelSequence labels;
  Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
  shapeTool->GetShapes(labels);

  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    TDF_Label label = labels.Value(i);
    Handle(TDataStd_Name) name = new TDataStd_Name();

    if (label.FindAttribute(TDataStd_Name::GetID(), name)) // retourne tous les attributs de type string
    { 
      TCollection_ExtendedString labelName = name->Get();
      if (labelName == TCollection_ExtendedString(oldName))
      {
        TDataStd_Name::Set(label, newName);
      }
    }
  }
  Py_INCREF(Py_None);
  return Py_None;
}

