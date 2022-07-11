/*
 * Filename: dip_distance.h
 *
 * (C) Copyright 1995-1997               Pattern Recognition Group
 *     All rights reserved               Faculty of Applied Physics
 *                                       Delft University of Technology
 *                                       Lorentzweg 1
 *                                       2628 CJ Delft
 *                                       The Netherlands
 *
 * Contact: Dr. ir. Lucas J. van Vliet
 * email : lucas@ph.tn.tudelft.nl
 *
 * Author: Geert M.P. van Kempen
 *
 */

#ifndef DIP_DISTANCE_H
#define DIP_DISTANCE_H
#ifdef __cplusplus
extern "C" {
#endif

typedef enum
{
	DIP_DISTANCE_EDT_FAST = 0,
	DIP_DISTANCE_EDT_TIES = 1,
	DIP_DISTANCE_EDT_TRUE = 2,
	DIP_DISTANCE_EDT_BRUTE_FORCE = 3
} dipf_DistanceTransform;

typedef struct
{
  dip_int x;
  dip_int y;
} dip_XYPosition;

typedef struct
{
  dip_int x;
  dip_int y;
  dip_int z;
} dip_XYZPosition;

#define DIP_DST_GDT_CHUNK_SIZE 100
#define DIP_DST_GDT_MAXDIST 50000

typedef enum {
  DIP_DBPS_MAX = 0,
  DIP_DBPS_INDEX_MAX,
  DIP_DBPS_COUNT_MAX,
  DIP_DBPS_MIN,
  DIP_DBPS_INDEX_MIN,
  DIP_DBPS_COUNT_MIN,
  DIP_DBPS_MEAN,
  DIP_DBPS0_DISTANCE,
  DIP_DBPSS_DIM_OFFSET,
  DIP_DBPSS_NO_FIELDS=DIP_DBPSS_DIM_OFFSET,
  DIP_DBPSS_NO_ELEMENTS=DIP_DBPSS_NO_FIELDS*2-1,
} dip_DistanceBetweenPointSetsIDs;

DIP_ERROR dip_EuclideanDistanceTransform    ( dip_Image, dip_Image,
	                                           dip_FloatArray, dip_Boolean, dipf_DistanceTransform );
DIP_ERROR dip_VectorDistanceTransform       ( dip_Image, dip_ImageArray,
	                                           dip_FloatArray, dip_Boolean, dipf_DistanceTransform );
DIP_ERROR dip_GreyWeightedDistanceTransform ( dip_Image, dip_Image, dip_Image,
	                                           dip_Image, dip_int, dip_IntegerArray, dip_FloatArray );
DIP_ERROR dip_FastMarching_PlaneWave (dip_Image,dip_Image,dip_Image,dip_Image,dip_Image,dip_Image,dip_Image,dip_Image,dip_Image);
DIP_ERROR dip_FastMarching_SphericalWave (dip_Image,dip_Image,dip_Image,dip_Image,dip_Image,dip_Image,dip_Image,dip_Image,dip_Image);
/* DIP_ERROR dip_FastMarching_PlanarCost (); */

DIP_ERROR dip_fm                            ( dip_Image, dip_Image,
                                              dip_IntegerArray, dip_IntegerArray*, dip_IntegerArray,
                                              dip_ImageArray*, dip_ImageArray*, dip_ImageArray* );
DIP_ERROR dip_fmgradientdecent              ( dip_Image, dip_int, dip_Image,
                                              dip_int*, dip_int*, dip_Resources );
DIP_ERROR dip_DistanceBetweenPointSets( dip_Image p1, dip_Image p2, dip_Image pdims, dip_ImageArray out, dip_float tolerance, dip_Boolean fsqrt );

#ifdef __cplusplus
}
#endif
#endif
