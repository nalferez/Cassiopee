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
// CAD trim a list of face with a list of edges

#include "occ.h"
#include "ShapeFix_Shape.hxx"
#include "ShapeFix_Wireframe.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepFeat_SplitShape.hxx"
#include "BRep_Builder.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "TopoDS.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepAlgoAPI_Splitter.hxx"

#include "XCAFDoc_DocumentTool.hxx"
#include "TDocStd_Document.hxx"
#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFDoc_ShapeMapTool.hxx"
#include "TDataStd_Name.hxx"

//=====================================================================
// Trim faces
// trim face1 set with face2 set
// if face2 is solid -> use cut (algo=0)
// if face2 is a set of faces -> splitter (algo=1)
// if mode=0, face2 cut face1
// if mode=1, face1 cut face2
// if mode=2, face1 cut face2 and face2 cut face1
//=====================================================================
PyObject* K_OCC::trimFaces(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listOfFaceNo1; PyObject* listOfFaceNo2;
  E_Int algo=0; E_Int mode=0;
  if (!PYPARSETUPLE_(args, OOO_ II_, &hook, &listOfFaceNo1, &listOfFaceNo2,
    &mode, &algo)) return NULL;

  GETPACKET;
  GETSHAPE;
  GETMAPSURFACES;
  
  // Build untouched faces list
  std::list<E_Int> pl;
  E_Int nfaces = surfaces.Extent();
  for (E_Int i = 1; i <= nfaces; i++) pl.push_back(i);

  // Build compound from face1
  BRep_Builder builder1;
  TopoDS_Compound compound1;
  builder1.MakeCompound(compound1);
  for (E_Int i = 0; i < PyList_Size(listOfFaceNo1); i++)
  {
    PyObject* faceNoO = PyList_GetItem(listOfFaceNo1, i);
    E_Int faceNo = PyLong_AsLong(faceNoO);
    auto it = std::find(pl.begin(), pl.end(), faceNo);
    if (it != pl.end()) pl.erase(it);
    const TopoDS_Face& F = TopoDS::Face(surfaces(faceNo));
    builder1.Add(compound1, F);
  }

  // Build compound from face2
  BRep_Builder builder2;
  TopoDS_Compound compound2;
  builder1.MakeCompound(compound2);
  for (E_Int i = 0; i < PyList_Size(listOfFaceNo2); i++)
  {
    PyObject* faceNoO = PyList_GetItem(listOfFaceNo2, i);
    E_Int faceNo = PyLong_AsLong(faceNoO);
    auto it = std::find(pl.begin(), pl.end(), faceNo);
    if (it != pl.end()) pl.erase(it);
    const TopoDS_Face& F = TopoDS::Face(surfaces(faceNo));
    builder2.Add(compound2, F);
  }  

  // For rebuild
  BRep_Builder builder3;
  TopoDS_Compound compound3;  
  builder3.MakeCompound(compound3);
#ifndef USEXCAF
  for (auto& i : pl) // untouched faces
  {
    TopoDS_Face F = TopoDS::Face(surfaces(i));
    builder3.Add(compound3, F);
  }
#endif

  if (algo == 0) // cut for closed shapes
  {
    // trim the compound with cut
    TopoDS_Shape trimmedCompound1 = BRepAlgoAPI_Cut(compound1, compound2);
    TopoDS_Shape trimmedCompound2 = BRepAlgoAPI_Cut(compound2, compound1);
    
    if (mode == 0 || mode == 2)
    {
      TopExp_Explorer expl1(trimmedCompound1, TopAbs_FACE);
      while (expl1.More())
      {
        TopoDS_Shape shape = expl1.Current();
        TopoDS_Face face = TopoDS::Face(shape);
        builder3.Add(compound3, face);
        expl1.Next();
      }
    }
    else
    {
      TopExp_Explorer expl1(compound1, TopAbs_FACE);
      while (expl1.More())
      {
        TopoDS_Shape shape = expl1.Current();
        TopoDS_Face face = TopoDS::Face(shape);
        builder3.Add(compound3, face);
        expl1.Next();
      }
    }
    if (mode == 1 || mode == 2)
    {
      TopExp_Explorer expl2(trimmedCompound2, TopAbs_FACE);
      while (expl2.More()) 
      {
        TopoDS_Shape shape = expl2.Current();
        TopoDS_Face face = TopoDS::Face(shape);
        builder3.Add(compound3, face);
        expl2.Next();
      }
    }
    else
    {
      TopExp_Explorer expl2(compound2, TopAbs_FACE);
      while (expl2.More()) 
      {
        TopoDS_Shape shape = expl2.Current();
        TopoDS_Face face = TopoDS::Face(shape);
        builder3.Add(compound3, face);
        expl2.Next();
      }
    }
  }
  else // splitter for face soup
  {
    // trim with splitter
    BRepAlgoAPI_Splitter splitter;
    TopTools_ListOfShape args;
    args.Append(compound1);
    TopTools_ListOfShape tools;
    tools.Append(compound2);
    splitter.SetTools(tools);
    splitter.SetArguments(args);
    splitter.Build();
    TopoDS_Shape trimmedCompound1 = splitter.Shape();
    BRepAlgoAPI_Splitter splitter2;
    TopTools_ListOfShape args2;
    args2.Append(compound2);
    TopTools_ListOfShape tools2;
    tools2.Append(compound1);
    splitter2.SetTools(tools2);
    splitter2.SetArguments(args2);
    splitter2.Build();
    TopoDS_Shape trimmedCompound2 = splitter2.Shape();

    if (mode == 0 || mode == 2)
    {
      TopExp_Explorer expl1(trimmedCompound1, TopAbs_FACE);
      while (expl1.More())
      {
        TopoDS_Shape shape = expl1.Current();
        TopoDS_Face face = TopoDS::Face(shape);
        builder3.Add(compound3, face);
        expl1.Next();
      }
    }
    else
    {
      TopExp_Explorer expl1(compound1, TopAbs_FACE);
      while (expl1.More())
      {
        TopoDS_Shape shape = expl1.Current();
        TopoDS_Face face = TopoDS::Face(shape);
        builder3.Add(compound3, face);
        expl1.Next();
      }
    }
    if (mode == 1 || mode == 2)
    {
      TopExp_Explorer expl2(trimmedCompound2, TopAbs_FACE);
      while (expl2.More()) 
      {
        TopoDS_Shape shape = expl2.Current();
        TopoDS_Face face = TopoDS::Face(shape);
        builder3.Add(compound3, face);
        expl2.Next();
      }
    }
    else
    {
      TopExp_Explorer expl2(compound2, TopAbs_FACE);
      while (expl2.More()) 
      {
        TopoDS_Shape shape = expl2.Current();
        TopoDS_Face face = TopoDS::Face(shape);
        builder3.Add(compound3, face);
        expl2.Next();
      }
    }
  }
  
#ifdef USEXCAF

  GETDOC;
  GETSHAPETOOL;

  // suppress somes faces
  std::map< E_Int, std::vector<E_Int> > label2Faces;
  std::map< E_Int, std::vector<E_Int> > label2Edges;
  getLabel2Edges(*doc, label2Edges);
  getLabel2Faces(*doc, label2Faces);

  for (E_Int no = 0; no < PyList_Size(listOfFaceNo1); no++)
  {
    PyObject* noFaceO = PyList_GetItem(listOfFaceNo1, no);
    E_Int noFace = PyInt_AsLong(noFaceO);
    for (size_t j = 0; j < label2Faces.size(); j++)
    {
      std::vector< E_Int >& f = label2Faces[j];
      f.erase(std::remove(f.begin(), f.end(), noFace), f.end());
    }
  }
  for (E_Int no = 0; no < PyList_Size(listOfFaceNo2); no++)
  {
    PyObject* noFaceO = PyList_GetItem(listOfFaceNo2, no);
    E_Int noFace = PyInt_AsLong(noFaceO);
    for (size_t j = 0; j < label2Faces.size(); j++)
    {
      std::vector< E_Int >& f = label2Faces[j];
      f.erase(std::remove(f.begin(), f.end(), noFace), f.end());
    }
  }
  copyTopShape2OCAF(*shape, label2Edges, label2Faces, *doc);

  // Add compound3 as a new shape
  TDF_Label label = shapeTool->AddShape(*shape);
  TDataStd_Name::Set(label, "trimmed");

  // back copy
  TopoDS_Shape* newshp = copyOCAF2TopShape(*doc);
  delete shape;
  SETSHAPE(newshp);
  Py_INCREF(Py_None);
  return Py_None;

#else
  TopoDS_Shape* newshp = new TopoDS_Shape(compound3);

  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after trim: Nb edges=%d\n", se->Extent());
  printf("INFO: after trim: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
#endif
}