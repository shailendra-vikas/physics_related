#!/usr/bin/env python
import sys
import time
import numpy
import scipy
from scipy.spatial import KDTree
from cosmology import Cosmology

"""
    Collection of function for calculating the pair counting of the two different catalog.
"""


def get_cart_coord(ra,dec,r):
    '''Convert spherical to cartisian for given ra, dec and r
    '''
    phi_rad=numpy.deg2rad(ra)
    theta_rad=numpy.deg2rad(90-dec)
    return numpy.transpose(numpy.asarray((r*numpy.cos(phi_rad)*numpy.sin(theta_rad),r*numpy.sin(phi_rad)*numpy.sin(theta_rad),r*numpy.cos(theta_rad)))) #Transposed for the KDTree input

def get_cord_for_cat(data):
    '''
    Convert the spherical cordinate to cartisian system for the catalog
    '''
    phi_rad=numpy.deg2rad(data['RA'])
    theta_rad=numpy.deg2rad(90-data['DEC'])
    cosmo=Cosmology(Om=0.26,Ol=0.74,w=-1.0,h=0.72)
    r=cosmo.Dc(data['Z'])
    return_cat=numpy.empty(len(data),dtype={'names':('X','Y','Z'),'formats':('f','f','f')})
    return_cat['X']=r*numpy.cos(phi_rad)*numpy.sin(theta_rad)
    return_cat['Y']=r*numpy.sin(phi_rad)*numpy.sin(theta_rad)
    return_cat['Z']=r*numpy.cos(theta_rad)
    return return_cat

def count_pairs_parallel(coord1,coord2,dist,cosmo=None,disttype='ang',prefix='distances'):
    ''' Pair counting fucntion between two catalogs. Supports parallel processing.
        coord1,coord2 : The catalog files
        dist : Maximum distance of th calculation
        cosmo: Instance of the cosmology, If None it is made using default parameters
        disttype : 'ang' default, 'proj' or  '3d'
        prefix : prefix of the result file
    '''
    try:
        from mpi4py import MPI
        comm=MPI.COMM_WORLD
        rank,size=comm.Get_rank(),comm.Get_size()
    except ImportError:
        rank,size=0,1

    #Check the file name are passed
    if not ( isinstance(coord1,str) and isinstance(coord2,str) ):
        print 'Expect numpy file for input'
        sys.exit(-1)
    cosmo= cosmo if cosmo is not None else Cosmology(Om=0.26,Ol=0.74,w=-1.0,h=0.72)
    coord1_data,coord2_data=map(numpy.load,(coord1,coord2)) #load the catalog files

    #For multi processing devide one of the catalog in smaller catalog for parellel pair counting
    chunk=len(coord2_data)/size
    if rank==size-1:
        coord2_data=coord2_data[rank*chunk:]
    else:
        coord2_data=coord2_data[rank*chunk:(rank+1)*chunk]

    # For pair counting as angular ditance
    if disttype=='ang':
        coord1_t=KDTree(get_cart_coord(coord1_data['RA'],coord1_data['DEC'],1.0))
        coord2_t=KDTree(get_cart_coord(coord2_data['RA'],coord2_data['DEC'],1.0))
        dist_cart=2.0*numpy.sin(numpy.deg2rad(dist)/2.0) # convert deg to cord distance
        coord1_t_i,coord2_t_i,distance= scipy.sparse.find(coord1_t.sparse_distance_matrix(coord2_t,dist_cart))
        distance=numpy.rad2deg(2.0*numpy.arcsin(distance/2.0))

    # For pair counting as projected distance
    if disttype=='proj':
        coord1_z=cosmo.Dc(coord1_data['Z'])
        coord2_z=cosmo.Dc(coord2_data['Z'])
        coord1_t=KDTree(get_cart_coord(coord1_data['RA'],coord1_data['DEC'],coord1_z))
        coord2_t=KDTree(get_cart_coord(coord2_data['RA'],coord2_data['DEC'],coord2_z))
        dist_cart=numpy.sqrt(2)*dist
        coord1_t_i,coord2_t_i,distance= scipy.sparse.find(coord1_t.sparse_distance_matrix(coord2_t,dist_cart))
        los_distance=numpy.abs(coord1_z[coord1_t_i]-coord2_z[coord2_t_i])
        proj_distance=numpy.sqrt(distance**2-los_distance**2)
        distance=numpy.transpose(numpy.array((proj_distance[proj_distance<=dist],los_distance[proj_distance<=dist])))
    # For pair counting as 3d distance
    else:
        coord1_z=cosmo.Dc(coord1_data['Z'])
        coord2_z=cosmo.Dc(coord2_data['Z'])
        coord1_t=KDTree(get_cart_coord(coord1_data['RA'],coord1_data['DEC'],coord1_z))
        coord2_t=KDTree(get_cart_coord(coord2_data['RA'],coord2_data['DEC'],coord2_z))
        coord1_t_i,coord2_t_i,distance= scipy.sparse.find(coord1_t.sparse_distance_matrix(coord2_t,dist))
    # Save the individual files for multiprocess.
    numpy.save('%s%i'%(prefix,rank),distance)

def count_pairs(cat1,cat2,dist,prefix='distances'):
    """ This method calculate the line of sight  and projected distance for the distance pair.
        cat1,cat2: are catalogs
        dist: maximum distance of the pair
        prefix: The prefix for the name of result file.
    """
    cat1,cat2=map(numpy.load,(cat1,cat2))
    # if not cartisian already than convert it to cartisian
    if 'X' not in cat1.dtype.names:
        cat1=get_cord_for_cat(cat1)
    if 'X' not in cat2.dtype.names:
        cat2=get_cord_for_cat(cat2)
    r1=numpy.sqrt(cat1['X']**2 + cat1['Y']**2 + cat1['Z']**2)
    r2=numpy.sqrt(cat2['X']**2 + cat2['Y']**2 + cat2['Z']**2)

    cat1_tree=KDTree(numpy.transpose(numpy.asarray((cat1['X'],cat1['Y'],cat1['Z'])) ) )
    del(cat1)
    cat2_tree=KDTree(numpy.transpose(numpy.asarray((cat2['X'],cat2['Y'],cat2['Z'])) ) )
    del(cat2)
    
    index1,index2,distances=scipy.sparse.find(cat1_tree.sparse_distance_matrix(cat2_tree,numpy.sqrt(2)*dist))
    los=(distances**2+numpy.abs(r1[index1]**2-r2[index2]**2))/(2*numpy.where(r1[index1]>r2[index2],r1[index1],r2[index2]))
    proj=numpy.sqrt(distances**2-los**2)
    

    dtype={'names':('INDEX1','INDEX2','LOS','PROJ'),'formats':(numpy.uint32,numpy.uint32,numpy.float64,numpy.float64)}
    result=numpy.empty(len(distances),dtype=dtype)
    result['INDEX1']=index1
    result['INDEX2']=index2
    result['LOS']=los
    result['PROJ']=proj

    numpy.save('%s'%prefix,result)

def get_options():
    from optparse import OptionParser,OptionGroup
    parser=OptionParser(usage='%prog [options] arg')
    parser.add_option('--prefix',default='distance',help='Prefix for the file created.')
    parser.add_option('--type',default='ang',help='Distance type "ang" or "comov" or "proj"')
    parser.add_option('--maxdist',default=2,type='float',help='Max distance')
    parser.add_option('--parallel',action='store_true',default=False,help='Use parallel function')
    options,args=parser.parse_args()
    if not args or len(args)<2:
       parser.print_help()
       sys.exit(-1) 

    if options.parallel:
        count_pairs_parallel(args[0],args[1],options.maxdist,prefix=options.prefix,disttype=options.type)
    else:
        count_pairs(args[0],args[1],options.maxdist,prefix=options.prefix)
    
if __name__=='__main__':
    get_options()
