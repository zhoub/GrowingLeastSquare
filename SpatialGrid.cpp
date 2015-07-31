#include <algorithm>

#include "SpatialGrid.h"

template < int N = 1 >
struct RadixSort
{
    bool operator()( const Imath::V2i & A , const Imath::V2i & B )
    {
        return A[N] < B[N] ;
    }
} ;

SpatialGrid::SpatialGrid( float CellSize )
:
mCellSize( CellSize ) ,
mNumCells( 0 )
{
    mBound.makeEmpty() ;
}

SpatialGrid::~SpatialGrid()
{
}

void SpatialGrid::Initialize( const Imath::V3f * PArray , size_t NumP )
{
    mIdxVec.reserve( NumP ) ;

    //
    for ( size_t i = 0 ; i < NumP ; ++ i )
    {
        mBound.extendBy( PArray[i] ) ;
    }

    mBound.min -= Imath::V3f( FLT_EPSILON ) ;
    mBound.max += Imath::V3f( FLT_EPSILON ) ;

    //
    const Imath::V3f & BoundSize = mBound.size() ;

    Imath::V3f NumCells = BoundSize / mCellSize ;
    mNumCells.x = static_cast< int >( ceilf( NumCells.x ) ) ;
    mNumCells.y = static_cast< int >( ceilf( NumCells.y ) ) ;
    mNumCells.z = static_cast< int >( ceilf( NumCells.z ) ) ;

    for ( size_t i = 0 ; i < NumP ; ++ i )
    {
        Imath::V3f FltIdx = ( PArray[i] - mBound.min ) / mCellSize ;
        Imath::V3i IntIdx( static_cast< int >( floorf( FltIdx.x ) + FLT_EPSILON ) ,
                           static_cast< int >( floorf( FltIdx.y ) + FLT_EPSILON ) ,
                           static_cast< int >( floorf( FltIdx.z ) + FLT_EPSILON ) ) ;
        int LinearIdx = IntIdx.x + IntIdx.y * mNumCells.x + IntIdx.z * mNumCells.x * mNumCells.y ;
        mIdxVec.push_back( Imath::V2i( i , LinearIdx ) ) ;
    }

    std::sort( mIdxVec.begin() , mIdxVec.end() , RadixSort< 1 >() ) ;
}

void SpatialGrid::Search( std::vector< int > & NeighborIdxVec ,  const Imath::V3f * PArray , const Imath::V3f & X , float R ) const
{
    Imath::V3f FltIdx = ( X - mBound.min ) / mCellSize ;
    Imath::V3i CentricCellIdx( static_cast< int >( floorf( FltIdx.x ) + FLT_EPSILON ) ,
                               static_cast< int >( floorf( FltIdx.y ) + FLT_EPSILON ) ,
                               static_cast< int >( floorf( FltIdx.z ) + FLT_EPSILON ) ) ;

    int CellRange = static_cast< int >( ceilf( R / mCellSize ) ) ;

    //
    for ( int k = CentricCellIdx.z - CellRange ;  k <= CentricCellIdx.z + CellRange ; ++ k )
    {
        if ( k < 0 || k > mNumCells.z - 1 )
        {
            continue ;
        }
        for ( int j = CentricCellIdx.y - CellRange ; j <= CentricCellIdx.y + CellRange ; ++ j )
        {
            if ( j < 0 || j > mNumCells.y - 1 )
            {
                continue ;
            }
            for ( int i = CentricCellIdx.x - CellRange ; i <= CentricCellIdx.x + CellRange ; ++ i )
            {
                if ( i < 0 || i > mNumCells.x - 1 )
                {
                    continue ;
                }
                int LinearIdx = i + j * mNumCells.x + k * mNumCells.x * mNumCells.y ;
                Imath::V2i CompoundIdx( - 1 , LinearIdx ) ;

                std::vector< Imath::V2i >::const_iterator LowerBoundIter = std::lower_bound( mIdxVec.begin() , mIdxVec.end() , CompoundIdx , RadixSort< 1 >() ) ;
                if ( LowerBoundIter->y != LinearIdx )
                {
                    continue ;
                }

                int Idx = - 1 ;
                while ( LowerBoundIter->y == LinearIdx )
                {
                    int Idx = LowerBoundIter->x ;
                    if ( ( PArray[Idx] - X ).length() < R )
                    {
                        NeighborIdxVec.push_back( Idx ) ;
                    }

                     ++ LowerBoundIter ;
                }
            }
        }
    }
}
