/**
 * Growing Least Square
 * Copyright (C) 2017 Bo Zhou<bo.schwarzstein@gmai.com>

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/


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
