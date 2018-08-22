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
import sys

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

def amun_dataset(fname, vname, shrink = 1, progress = False):
  '''
      Subroutine to reads dataset from AMUN HDF5 snapshots.

      Arguments:

        fname    - the AMUN HDF5 file name;
        vname    - the variable name;
        shrink   - the shrink factor (must be the power of 2 and not larger
                   than the block size);
        progress - the progress bar switch;

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

  # get the file path
  #
  dname = op.dirname(fname)

  if progress:
    sys.stdout.write("Data file path:\n  '%s'\n" % (dname))

  # get attributes necessary to reconstruct the domain
  #
  g  = f['attributes'].attrs

  # get the set of equations used to perform the simulation
  # and the equation of state
  #
  eqsys = g.get('eqsys')[0].astype(str)
  eos   = g.get('eos')[0].astype(str)

  # get the snapshot number, the number of domain files, and the number of blocks
  #
  nr = g.get('isnap')[0]
  nc = g.get('nprocs')[0]
  nl = g.get('nleafs')[0]
  if eos == 'adi':
    gm = g.get('gamma')[0]

  # build the list of supported variables
  #
  variables = []
  for var in f['variables'].keys():
    variables.append(var)

  # add derived variables if possible
  #
  if 'velx' in variables and 'vely' in variables and 'velz' in variables:
    variables.append('velo')
    variables.append('divv')
    variables.append('vort')
  if 'magx' in variables and 'magy' in variables and 'magz' in variables:
    variables.append('magn')
    variables.append('divb')
    variables.append('curr')
  if (eqsys == 'hd' or eqsys == 'mhd') and eos == 'adi' \
                    and 'pres' in variables:
    variables.append('eint')
  if 'dens' in variables and 'pres' in variables:
    variables.append('temp')
  if (eqsys == 'hd' or eqsys == 'mhd') \
                    and 'dens' in variables \
                    and 'velx' in variables \
                    and 'vely' in variables \
                    and 'velz' in variables:
    variables.append('ekin')
  if (eqsys == 'mhd' or eqsys == 'srmhd') \
                     and 'magx' in variables \
                     and 'magy' in variables \
                     and 'magz' in variables:
    variables.append('emag')
  if eqsys == 'hd' and 'ekin' in variables and 'eint' in variables:
    variables.append('etot')
  if eqsys == 'mhd' and 'eint' in variables \
                    and 'ekin' in variables \
                    and 'emag' in variables:
    variables.append('etot')
  if (eqsys == 'srhd' or eqsys == 'srmhd') and 'velo' in variables:
    variables.append('lore')

  # check if the requested variable is in the variable list
  #
  if not vname in variables:
    print('The requested variable cannot be extracted from the file datasets!')
    return False

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

  f.close()

  # check if the shrink parameter is correct
  #
  sh = bm[0:ndims].min()
  while(sh > shrink):
    sh /= 2
  shrink = int(sh)

  # prepare dimensions of the output array and allocate it
  #
  dm = np.array(rm[0:ndims] * bm[0:ndims] * 2**(ml - 1) / shrink, \
                                                      dtype = np.int32)
  ret = np.zeros(dm[::-1])

  # iterate over all subdomain files
  #
  nb = 0
  for n in range(nc):
    fname = 'p%06d_%05d.h5' % (nr, n)
    lname = op.join(dname, fname)
    f = h5.File(lname, 'r')
    g  = f['attributes'].attrs
    dblocks = g.get('dblocks')[0]
    if dblocks > 0:
      g  = f['coordinates']
      levels  = g['levels'][()]
      coords  = g['coords'][()]
      dx      = g['dx'][()]
      dy      = g['dy'][()]
      dz      = g['dz'][()]
      g       = f['variables']
      if vname == 'velo':
        dataset = np.sqrt(g['velx'][:,:,:,:]**2 \
                        + g['vely'][:,:,:,:]**2 \
                        + g['velz'][:,:,:,:]**2)
      elif vname == 'magn':
        dataset = np.sqrt(g['magx'][:,:,:,:]**2 \
                        + g['magy'][:,:,:,:]**2 \
                        + g['magz'][:,:,:,:]**2)
      elif vname == 'eint':
        dataset = 1.0 / (gm - 1.0) * g['pres'][:,:,:,:]
      elif vname == 'ekin':
        dataset = 0.5 * g['dens'][:,:,:,:] * (g['velx'][:,:,:,:]**2 \
                                            + g['vely'][:,:,:,:]**2 \
                                            + g['velz'][:,:,:,:]**2)
      elif vname == 'emag':
        dataset = 0.5 * (g['magx'][:,:,:,:]**2 \
                       + g['magy'][:,:,:,:]**2 \
                       + g['magz'][:,:,:,:]**2)
      elif vname == 'etot':
        dataset = 1.0 / (gm - 1.0) * g['pres'][:,:,:,:] \
                + 0.5 * g['dens'][:,:,:,:] * (g['velx'][:,:,:,:]**2 \
                                            + g['vely'][:,:,:,:]**2 \
                                            + g['velz'][:,:,:,:]**2)
        if eqsys == 'mhd':
          dataset += 0.5 * (g['magx'][:,:,:,:]**2 \
                          + g['magy'][:,:,:,:]**2 \
                          + g['magz'][:,:,:,:]**2)
      elif vname == 'temp':
        dataset = g['pres'][:,:,:,:] / g['dens'][:,:,:,:]
      elif vname == 'lore':
        dataset = 1.0 / np.sqrt(1.0 - (g['velx'][:,:,:,:]**2 \
                                     + g['vely'][:,:,:,:]**2 \
                                     + g['velz'][:,:,:,:]**2))
      elif vname == 'divv':
        dataset = np.zeros(g['velx'].shape)
        fields  = [ 'velx', 'vely', 'velz' ]
        h       = (dx, dy, dz)
        for i in range(ndims):
          v = fields[i]
          dataset += 0.5 * (np.roll(g[v][:,:,:,:], -1, axis = 2)  \
                          - np.roll(g[v][:,:,:,:],  1, axis = 2)) \
                                                    / h[i][levels[:] - 1]
      elif vname == 'divb':
        dataset = np.zeros(g['magx'].shape)
        fields  = [ 'magx', 'magy', 'magz' ]
        h       = (dx, dy, dz)
        for i in range(ndims):
          v = fields[i]
          dataset += 0.5 * (np.roll(g[v][:,:,:,:], -1, axis = 2)  \
                          - np.roll(g[v][:,:,:,:],  1, axis = 2)) \
                                                    / h[i][levels[:] - 1]
      elif vname == 'vort':
        if ndims == 3:
          wx = 0.5 * (np.roll(g['velz'][:,:,:,:], -1, axis = 1)  \
                    - np.roll(g['velz'][:,:,:,:],  1, axis = 1)) \
                                                / dy[levels[:]-1] \
             - 0.5 * (np.roll(g['vely'][:,:,:,:], -1, axis = 0)  \
                    - np.roll(g['vely'][:,:,:,:],  1, axis = 0)) \
                                                / dz[levels[:]-1]
          wy = 0.5 * (np.roll(g['velx'][:,:,:,:], -1, axis = 0)  \
                    - np.roll(g['velx'][:,:,:,:],  1, axis = 0)) \
                                                / dz[levels[:]-1] \
             - 0.5 * (np.roll(g['velz'][:,:,:,:], -1, axis = 2)  \
                    - np.roll(g['velz'][:,:,:,:],  1, axis = 2)) \
                                                / dx[levels[:]-1]
          wz = 0.5 * (np.roll(g['vely'][:,:,:,:], -1, axis = 2)  \
                    - np.roll(g['vely'][:,:,:,:],  1, axis = 2)) \
                                                / dx[levels[:]-1] \
             - 0.5 * (np.roll(g['velx'][:,:,:,:], -1, axis = 1)  \
                    - np.roll(g['velx'][:,:,:,:],  1, axis = 1)) \
                                                / dy[levels[:]-1]
        else:
          wx =   0.5 * (np.roll(g['velz'][:,:,:,:], -1, axis = 1)  \
                      - np.roll(g['velz'][:,:,:,:],  1, axis = 1)) \
                                                 / dy[levels[:]-1]
          wy = - 0.5 * (np.roll(g['velz'][:,:,:,:], -1, axis = 2)  \
                      - np.roll(g['velz'][:,:,:,:],  1, axis = 2)) \
                                                 / dx[levels[:]-1]
          wz =   0.5 * (np.roll(g['vely'][:,:,:,:], -1, axis = 2)  \
                      - np.roll(g['vely'][:,:,:,:],  1, axis = 2)) \
                                                 / dx[levels[:]-1] \
               - 0.5 * (np.roll(g['velx'][:,:,:,:], -1, axis = 1)  \
                      - np.roll(g['velx'][:,:,:,:],  1, axis = 1)) \
                                                 / dy[levels[:]-1]
        dataset = np.sqrt(wx * wx + wy * wy + wz * wz)
      elif vname == 'curr':
        if ndims == 3:
          wx = 0.5 * (np.roll(g['magz'][:,:,:,:], -1, axis = 1)  \
                    - np.roll(g['magz'][:,:,:,:],  1, axis = 1)) \
                                                / dy[levels[:]-1] \
             - 0.5 * (np.roll(g['magy'][:,:,:,:], -1, axis = 0)  \
                    - np.roll(g['magy'][:,:,:,:],  1, axis = 0)) \
                                                / dz[levels[:]-1]
          wy = 0.5 * (np.roll(g['magx'][:,:,:,:], -1, axis = 0)  \
                    - np.roll(g['magx'][:,:,:,:],  1, axis = 0)) \
                                                / dz[levels[:]-1] \
             - 0.5 * (np.roll(g['magz'][:,:,:,:], -1, axis = 2)  \
                    - np.roll(g['magz'][:,:,:,:],  1, axis = 2)) \
                                                / dx[levels[:]-1]
          wz = 0.5 * (np.roll(g['magy'][:,:,:,:], -1, axis = 2)  \
                    - np.roll(g['magy'][:,:,:,:],  1, axis = 2)) \
                                                / dx[levels[:]-1] \
             - 0.5 * (np.roll(g['magx'][:,:,:,:], -1, axis = 1)  \
                    - np.roll(g['magx'][:,:,:,:],  1, axis = 1)) \
                                                / dy[levels[:]-1]
        else:
          wx =   0.5 * (np.roll(g['magz'][:,:,:,:], -1, axis = 1)  \
                      - np.roll(g['magz'][:,:,:,:],  1, axis = 1)) \
                                                 / dy[levels[:]-1]
          wy = - 0.5 * (np.roll(g['magz'][:,:,:,:], -1, axis = 2)  \
                      - np.roll(g['magz'][:,:,:,:],  1, axis = 2)) \
                                                 / dx[levels[:]-1]
          wz =   0.5 * (np.roll(g['magy'][:,:,:,:], -1, axis = 2)  \
                      - np.roll(g['magy'][:,:,:,:],  1, axis = 2)) \
                                                 / dx[levels[:]-1] \
               - 0.5 * (np.roll(g['magx'][:,:,:,:], -1, axis = 1)  \
                      - np.roll(g['magx'][:,:,:,:],  1, axis = 1)) \
                                                 / dy[levels[:]-1]
        dataset = np.sqrt(wx * wx + wy * wy + wz * wz)
      else:
        dataset = g[vname][:,:,:,:]

      # rescale all blocks to the effective resolution
      #
      for l in range(dblocks):
        nn   = 2**(ml - levels[l])
        cm   = np.array(bm[0:ndims] * nn / shrink, dtype = np.int32)
        ibeg = coords[0:ndims,l] * cm[0:ndims]
        iend = ibeg + cm[0:ndims]
        if ndims == 3:
          ib, jb, kb = ibeg[0], ibeg[1], ibeg[2]
          ie, je, ke = iend[0], iend[1], iend[2]
          ret[kb:ke,jb:je,ib:ie] = rebin(dataset[ng:-ng,ng:-ng,ng:-ng,l], cm)
        else:
          ib, jb = ibeg[0], ibeg[1]
          ie, je = iend[0], iend[1]
          ret[jb:je,ib:ie]       = rebin(dataset[     0,ng:-ng,ng:-ng,l], cm)

        nb += 1

        # print progress bar if desired
        #
        if progress:
          sys.stdout.write('\r')
          sys.stdout.write("Reading '%s' from '%s': block %d from %d" \
                        % (vname, fname, nb, nl))
          sys.stdout.flush()

    f.close()

  if (progress):
    sys.stdout.write('\n')
    sys.stdout.flush()

  # return the dataset
  #
  return ret

def rebin(a, newshape):
  '''
      Subroutine changes the size of the input array to to new shape,
      by copying cells or averaging them.
  '''
  assert len(a.shape) == len(newshape)

  m = a.ndim - 1
  if (a.shape[m] > newshape[m]):
    if a.ndim == 3:
      nn = [newshape[0], a.shape[0] / newshape[0],
            newshape[1], a.shape[1] / newshape[1],
            newshape[2], a.shape[2] / newshape[2]]
      return a.reshape(nn).mean(5).mean(3).mean(1)
    else:
      nn = [newshape[0], a.shape[0] / newshape[0],
            newshape[1], a.shape[1] / newshape[1]]
      return a.reshape(nn).mean(3).mean(1)
  else:
    for n in range(a.ndim):
      a = np.repeat(a, newshape[n] / a.shape[n], axis = n)
    return(a)


def amun_integrals(field, filename, pathlist):
    '''
        get_integral: iterate over pathlist and read and merge field values from filename files in the provided paths
    '''
    # Initiate the return values with empty array and file number.
    #
    vals = array([])
    num  = 1

    # Iterate over all paths provided in the list 'pathlist'.
    #
    for path in pathlist:

      # Iterate over all files in the current path.
      #
      while True:

        # Generate file name.
        #
        dfile = path + '/' + filename + '_' + str(num).zfill(2) + '.dat'

        # Check if the file exists.
        #
        if isfile(dfile):

          # Read values from the current integrals file.
          #
          lvals = read_integrals(dfile, field)

          # Append to the return array.
          #
          vals = append(vals, lvals)

          # Increase the number file.
          #
          num = num + 1

        else:

          # File does not exists, so go to the next path.
          #
          break

    # Return appended values.
    #
    return vals

def read_integrals(filename, column):
    '''
        read_integrals: reads a given column from an integral file.
    '''
    # Open the given file and check if it is text file.
    #
    f = open(filename, 'r')

    # Read fist line and store it in h, since it will be used to obtain the
    # column headers.
    #
    l = f.readline()
    h = l

    # Read first line which is not comment in order to determine the number of
    # columns and store the number of columns in nc.  Calculate the column width
    # and store it in wc.
    #
    while l.startswith('#'):
        l = f.readline()
    nc = len(l.rsplit())
    wc = int((len(l) - 9) / (nc - 1))

    # Split header line into a list.
    #
    lh = [h[1:9].strip()]
    for i in range(nc - 1):
      ib = i * wc + 10
      ie = ib + wc - 1
      lh.append(h[ib:ie].strip())


    ic = lh.index(column)

    # Read given column.
    #
    if (ic > -1):
      lc = [float(l.split()[ic])]
      for l in f:
        lc.append(float(l.split()[ic]))

    # Close the file.
    #
    f.close()

    # Return values.
    #
    return(array(lc))

if __name__ == "__main__":
  fname = './p000030_00000.h5'

  ret = amun_attribute(fname, 'time')
  print(ret)
  ret = amun_attribute(fname, 'dims')
  print(ret)

  ret = amun_dataset(fname, 'dens')
  print(ret.shape, ret.min(), ret.max())
