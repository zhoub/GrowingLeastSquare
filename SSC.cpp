#include <iostream>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include <OpenEXR/ImathVec.h>

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

void SolveU( float * U , const Imath::V3f & X , float RSearch , float RPoint , std::vector< Imath::V3f > & PVec , const std::vector< Imath::V3f > & NVec , const std::vector< float > & RVec )
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

    //

}

int main( int Argc , char * Argv[] )
{
    const int NumSamples = 16 ;

    const Imath::V3f Center( 0.5f , 1.5f , 2.5f ) ;
    const float R = 100.0f ;

    //
    std::vector< Imath::V3f > PVec ;
    std::vector< Imath::V3f > NVec ;

    GenerateSamples( PVec , NVec , NumSamples , Center , R ) ;

    //
    float U[5] = { 0.0f , 0.0f , 0.0f , 0.0f , 0.0f } ;
    for ( int i = 0 ; i < NumSamples ; ++ i )
    {
    }

    return EXIT_SUCCESS ;
}
