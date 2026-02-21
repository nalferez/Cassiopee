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

# include "transform.h"
using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;

// ============================================================================
/* Conservative smoothing */
// ============================================================================
PyObject* K_TRANSFORM::consSmooth(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array))
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =  K_ARRAY::getFromArray3(array, varString,
                                      f, im, jm, km, cn, eltType);

  if (res == 1) 
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: array must be unstructured.");
    return NULL;
  }

  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: array is invalid.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  PyObject* tpl = K_ARRAY::buildArray3(*f, varString, *cn, eltType, f->getApi());
  FldArrayF* fo; FldArrayI* cno;
  K_ARRAY::getFromArray3(array,fo, cno);

  RELEASESHAREDU(array, f, cn); 
  RELEASESHAREDU(tpl, fo, cno); 
  return tpl;
}