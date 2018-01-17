"""
================================================================================

  This file is part of the AMUN source code, a program to perform
  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
  adaptive mesh.

  Copyright (C) 2018 Grzegorz Kowal <grzegorz@amuncode.org>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

================================================================================

 module: AMUN

  Python module with subroutines to read AMUN code HDF5 files.

  The only requirements for this module are:
     - h5py
     - numpy

--------------------------------------------------------------------------------
"""
import h5py    as h5
import numpy   as np
import os.path as op

def amun_attribute(fname, aname):
  '''
      Subroutine to reads global attributes from AMUN HDF5 snapshots.

      Arguments:

        fname - the AMUN HDF5 file name;
        aname - the attribute name;

      Return values:

        ret   - the value read from the attribute;

      Examples:

        time = amun_attribute('p000010_00000.h5', 'time')

  '''
  # open the file
  #
  f = h5.File(fname, 'r')

  # check if the file is written in the AMUN format
  #
  if f.attrs.get('code')[0].astype(str) != "AMUN" or \
        not 'attributes' in f or \
        not 'coordinates' in f or \
        not 'variables' in f:
    print('It seems this HDF5 file is corrupted or not compatible with the AMUN format!')
    f.close()
    return False

  # open the group of attributes
  #
  g  = f['attributes'].attrs

  # get attribute's value
  #
  ret = g.get(aname)
  if len(ret) == 1:
    ret = ret[0]

  # close the file
  #
  f.close()

  # return the value of attribute
  #
  return ret

def amun_dataset(fname, vname):
  '''
      Subroutine to reads dataset from AMUN HDF5 snapshots.

      Arguments:

        fname - the AMUN HDF5 file name;
        vname - the variable name;

      Return values:

        ret   - the array of values for the variable;

      Examples:

        dn = amun_dataset('p000010_00000.h5', 'dens')

  '''
  # open the file
  #

  f = h5.File(fname, 'r')

  # check if the file is written in the AMUN format
  #
  if f.attrs.get('code')[0].astype(str) != "AMUN" or \
        not 'attributes' in f or \
        not 'coordinates' in f or \
        not 'variables' in f:
    print('It seems this HDF5 file is corrupted or not compatible with the AMUN format!')
    return False

  # build the list of supported variables
  #
  variables = []
  for var in f['variables'].keys():
    variables.append(var)

  # check if the requested variable is in the variable list
  #
  if not vname in variables:
    print('The requested variable cannot be extracted from the file datasets!')
    return False

  # get attributes necessary to reconstruct the domain
  #
  g  = f['attributes'].attrs

  # get the snapshot number, the number of domain files, and the number of blocks
  #
  nr = g.get('isnap')[0]
  nc = g.get('nprocs')[0]
  nl = g.get('nleafs')[0]

  # prepare array to hold data
  #
  bm = g.get('dims')
  if bm[2] > 1:
    ndims = 3
  else:
    ndims = 2
  rm = g.get('rdims')
  ng = g.get('nghosts')
  ml = g.get('maxlev')[0]
  dm = rm[0:ndims] * bm[0:ndims] * 2**(ml - 1)
  ret = np.zeros(dm[::-1])

  f.close()

  # iterate over all subdomain files
  #
  dname = op.dirname(fname)
  for n in range(nc):
    lname = op.join(dname, 'p%06d_%05d.h5' % (nr, n))
    f = h5.File(lname, 'r')
    g  = f['attributes'].attrs
    dblocks = g.get('dblocks')[0]
    if dblocks > 0:
      g  = f['coordinates']
      levels  = g['levels'][()]
      coords  = g['coords'][()]
      g  = f['variables']
      if ndims == 3:
        dataset = g[vname][ng:-ng,ng:-ng,ng:-ng,:]
        for l in range(dblocks):
          nn = 2**(ml - levels[l])
          cm = bm * nn
          ibeg = coords[:,l] * cm[:]
          iend = ibeg + cm[:]
          ib, jb, kb = ibeg[0], ibeg[1], ibeg[2]
          ie, je, ke = iend[0], iend[1], iend[2]
          ret[kb:ke,jb:je,ib:ie] = np.repeat(np.repeat(np.repeat(dataset[:,:,:,l], nn, axis = 0), nn, axis = 1), nn, axis = 2)
      else:
        dataset = g[vname][0,ng:-ng,ng:-ng,:]
        for l in range(dblocks):
          nn = 2**(ml - levels[l])
          cm = bm[0:ndims] * nn
          ibeg = coords[:,l] * cm[:]
          iend = ibeg + cm[:]
          ib, jb = ibeg[0], ibeg[1]
          ie, je = iend[0], iend[1]
          ret[jb:je,ib:ie] = np.repeat(np.repeat(dataset[:,:,l], nn, axis = 0), nn, axis = 1)
    f.close()

  return ret

if __name__ == "__main__":
  fname = './p000030_00000.h5'

  ret = amun_attribute(fname, 'time')
  print(ret)
  ret = amun_attribute(fname, 'dims')
  print(ret)

  ret = amun_dataset(fname, 'dens')
  print(ret.shape, ret.min(), ret.max())