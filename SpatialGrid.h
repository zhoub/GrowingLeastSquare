#pragma once

#include <vector>

#include <OpenEXR/ImathVec.h>
#include <OpenEXR/ImathBox.h>

class SpatialGrid
{
public :

    SpatialGrid( float CellSize ) ;

    ~SpatialGrid() ;

    void Initialize( const Imath::V3f * PArray , size_t NumP ) ;

    void Search( std::vector< int > & NeighborIdxVec ,  const Imath::V3f * PArray , const Imath::V3f & X , float R ) const ;

private :

    float mCellSize ;

    Imath::Box3f mBound ;

    Imath::V3i mNumCells ;
    std::vector< Imath::V2i > mIdxVec ;
} ;
