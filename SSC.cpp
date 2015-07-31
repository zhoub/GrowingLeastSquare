#include <Windows.h>

#include <iostream>
#include <iterator>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include <OpenEXR/ImathVec.h>

#include <sfcnn.hpp>

#include "SpatialGrid.h"

void GenerateSamples( std::vector< Imath::V3f > & PVec , std::vector< Imath::V3f > & NVec , int NumSamples , const Imath::V3f & Center , float R )
{
    for ( int i = 0 ; i < NumSamples ; ++ i )
    {
        const float Theta = static_cast< float >( rand() ) / RAND_MAX * static_cast< float >( M_PI ) ;
        const float Phi   = static_cast< float >( rand() ) / RAND_MAX * static_cast< float >( M_PI ) * 2.0f ;

        Imath::V3f XYZ( R * cosf( Theta ) * sinf( Phi ) ,
                        R * sinf( Theta ) * sinf( Phi ) ,
                        R * cosf( Phi ) ) ;
        Imath::V3f Offset( 0.1f * R * ( static_cast< float >( rand() ) / RAND_MAX * 2.0f - 1.0f ) ,
                           0.1f * R * ( static_cast< float >( rand() ) / RAND_MAX * 2.0f - 1.0f ) ,
                           0.1f * R * ( static_cast< float >( rand() ) / RAND_MAX * 2.0f - 1.0f ) ) ;

        const Imath::V3f & P = XYZ + Center + Offset ;
        const Imath::V3f & N = XYZ.normalize() ;
        PVec.push_back( P ) ;
        NVec.push_back( N ) ;
    }
}

float EvaluateWeight( const Imath::V3f X , const Imath::V3f & P , float R , float H )
{
    float x = ( P - X ).length() / ( R * H ) ;
    if ( x < 1.0f )
    {
        return powf( 1 - x * x , 4.0f ) ;
    }
    else
    {
        return 0.0f ;
    }
}

Imath::V3f Summate( const std::vector< float > & WVec , const std::vector< Imath::V3f > & XVec )
{
    Imath::V3f Result( 0.0f ) ;
    for ( std::vector< float >::size_type i = 0 ; i < WVec.size() ; ++ i )
    {
        Result += WVec[i] * XVec[i] ;
    }
    return Result ;
}

float Summate( const std::vector< float > & WVec , const std::vector< Imath::V3f > & XVec , const std::vector< Imath::V3f > & YVec )
{
    float Result = 0.0f ;
    for ( std::vector< float >::size_type i = 0 ; i < WVec.size() ; ++ i )
    {
        Result += WVec[i] * XVec[i].dot( YVec[i] ) ;
    }
    return Result ;
}

void SolveU( float * U , float * UHat , const Imath::V3f & X , float RSearch , float RPoint , std::vector< Imath::V3f > & PVec , const std::vector< Imath::V3f > & NVec )
{
    //
    std::vector< Imath::V3f > PNearVec ;
    std::vector< Imath::V3f > NNearVec ;

    float TotalW = 0.0f ;
    std::vector< float > WVec ;

    for ( std::vector< Imath::V3f >::size_type i = 0 ; i < PVec.size() ; ++ i )
    {
        float Distance = ( PVec[i] - X ).length() ;
        if ( Distance < RSearch )
        {
            PNearVec.push_back( PVec[i] ) ;
            NNearVec.push_back( NVec[i] ) ;

            float W = EvaluateWeight( X , PVec[i] , RPoint , 1.0f ) ;
            WVec.push_back( W ) ;
            TotalW += W ;
        }
    }

    std::vector< float > WNVec( WVec ) ;
    for ( std::vector< float >::size_type i = 0 ; i < WNVec.size() ; ++ i )
    {
        WNVec[i] /= TotalW ;
    }

    if ( ! PNearVec.size() )
    {
        exit( EXIT_FAILURE ) ;
    }

    //
    float UQuadric = 0.5f
                   * ( Summate( WVec , PNearVec , NNearVec ) - Summate( WNVec , PNearVec ) .dot( Summate( WVec , NNearVec ) ) )
                   / ( Summate( WVec , PNearVec , PNearVec ) - Summate( WNVec , PNearVec ) .dot( Summate( WVec , PNearVec ) ) ) ;

    Imath::V3f ULinear = Summate( WNVec , NVec ) - 2.0f * UQuadric * Summate( WNVec , PNearVec ) ;

    float UConstant = - ULinear.dot( Summate( WNVec , PNearVec ) ) - UQuadric * Summate( WNVec , PNearVec , PNearVec ) ;

    //
    U[0] = UConstant ;
    U[1] = ULinear.x ;
    U[2] = ULinear.y ;
    U[3] = ULinear.z ;
    U[4] = UQuadric ;

    //
    float ULength = sqrtf( ULinear.length2() - 4.0f * UConstant * UQuadric ) ;
    for ( int i = 0 ; i < 5 ; ++ i )
    {
        UHat[i] = U[i] / ULength ;
    }
}

void FitSphere()
{
    const int NumSamples = 10000 ;

    const Imath::V3f Center( 0 ) ;
    const float R = 100.0f ;

    //
    std::vector< Imath::V3f > PVec ;
    std::vector< Imath::V3f > NVec ;

    GenerateSamples( PVec , NVec , NumSamples , Center , R ) ;

    //
    const Imath::V3f & X = PVec.front() ;

    float U[5] = { 0.0f , 0.0f , 0.0f , 0.0f , 0.0f } ;
    float UHat[5] = { 0.0f , 0.0f , 0.0f , 0.0f , 0.0f } ;
    SolveU( U , UHat , X , 5.0f , 10.0f , PVec , NVec ) ;

    std::cout << "U : " << std::endl ;
    std::copy(     U ,   U + 5 , std::ostream_iterator< float >( std::cout , " " ) ) ;
    std::cout << std::endl ;

    std::cout << "UHat : " << std::endl ;
    std::copy( UHat , UHat + 5 , std::ostream_iterator< float >( std::cout , " " ) ) ;
    std::cout << std::endl ;

    float S = U[0] * 1.0f +
              U[1] * X.x +
              U[2] * X.y +
              U[3] * X.z +
              U[4] * X.length2() ;
    std::cout << "S(x) = " << S << std::endl ;

    float Ka = 2.0f * UHat[4] ;
    std::cout << "Ka = " << Ka << std::endl ;
}

void SearchPoint()
{
    const float CellSize = 2.5f ;

    //
    const int NumSamples = 50000 ;

    std::vector< Imath::V3f > PVec ; PVec.reserve( NumSamples ) ;
    for ( int i = 0 ; i < NumSamples ; ++ i )
    {
        Imath::V3f P( static_cast< float >( rand() ) / RAND_MAX * 10.0f ,
                      static_cast< float >( rand() ) / RAND_MAX * 10.0f ,
                      static_cast< float >( rand() ) / RAND_MAX * 10.0f ) ;

        PVec.push_back( P ) ;
    }

    //
    Imath::V3f X( 0.5f , 1.5f , 2.5f ) ;
    float R = 1.0f ;

    //
    LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
    LARGE_INTEGER Frequency;



    //
    std::vector< unsigned long > KNNIdxVec ;

    sfcnn< Imath::V3f , 3 , float > ANN( & PVec[0] , NumSamples ) ;
    std::vector< double > DistVec ;

    QueryPerformanceFrequency(&Frequency); 
    QueryPerformanceCounter(&StartingTime);
    ANN.ksearch( X , 200 , KNNIdxVec , DistVec ) ;

    //
    QueryPerformanceCounter(&EndingTime);
    ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    std::cout << ElapsedMicroseconds.QuadPart << std::endl ;

    //
    std::vector< double >::const_iterator LastPos = std::upper_bound( DistVec.begin() , DistVec.end() , static_cast< double >( R ) , std::less< double >() ) ;
    std::sort( KNNIdxVec.begin() , KNNIdxVec.begin() + ( LastPos - DistVec.begin() ) , std::less< unsigned long >() ) ;
    std::copy( KNNIdxVec.begin() , KNNIdxVec.begin() + ( LastPos - DistVec.begin() ) , std::ostream_iterator< int >( std::cout , "," ) ) ;
    std::cout << std::endl << std::endl ;

    //

    std::auto_ptr< SpatialGrid > SG( new SpatialGrid( CellSize ) ) ;
    SG->Initialize( & PVec[0] , NumSamples ) ;

    QueryPerformanceFrequency(&Frequency); 
    QueryPerformanceCounter(&StartingTime);
    std::vector< int > NeighborIdxVec ;
    SG->Search( NeighborIdxVec , &PVec[0] , X , R ) ;

    QueryPerformanceCounter(&EndingTime);
    ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    std::cout << ElapsedMicroseconds.QuadPart << std::endl ;

    std::sort( NeighborIdxVec.begin() , NeighborIdxVec.end() , std::less< int >() ) ;
    std::copy( NeighborIdxVec.begin() , NeighborIdxVec.end() , std::ostream_iterator< int >( std::cout , "," ) ) ;
}

int main( int Argc , char * Argv[] )
{
    SearchPoint() ;

    return EXIT_SUCCESS ;
}
