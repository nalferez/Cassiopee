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
#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFDoc_DocumentTool.hxx"
#include "XCAFDoc_ShapeMapTool.hxx"
#include "TDocStd_Document.hxx"
#include "TDataStd_Name.hxx"

//=====================================================================
// return label name if attribute exits, return "unknown" otherwise
//=====================================================================
E_Int K_OCC::getLabelName(TDF_Label& label, TCollection_ExtendedString& labelName)
{
  labelName = "Unnamed";
  Handle(TDataStd_Name) name = new TDataStd_Name();
  if (label.FindAttribute(TDataStd_Name::GetID(), name))
  { 
    labelName = name->Get();
    return 0;
  }
  return 1;
}