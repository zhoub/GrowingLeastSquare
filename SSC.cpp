#include <Windows.h>

#include <iostream>
#include <iterator>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include <OpenEXR/ImathVec.h>

#include <sfcnn.hpp>

#include <pointcloud.h>

#include <omp.h>

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

void SolveU( float * U , float * UHat , const Imath::V3f & X , float RSearch , float RPoint , float H , std::vector< Imath::V3f > & PNearVec , const std::vector< Imath::V3f > & NNearVec )
{
    float TotalW = 0.0f ;
    std::vector< float > WVec ;

    for ( std::vector< Imath::V3f >::size_type i = 0 ; i < PNearVec.size() ; ++ i )
    {
        float W = EvaluateWeight( X , PNearVec[i] , RPoint , H ) ;
        WVec.push_back( W ) ;
        TotalW += W ;
    }

    std::vector< float > WNVec( WVec ) ;
    for ( std::vector< float >::size_type i = 0 ; i < WNVec.size() ; ++ i )
    {
        WNVec[i] /= TotalW ;
    }

    //
    float UQuadric = 0.5f
                   * ( Summate( WVec , PNearVec , NNearVec ) - Summate( WNVec , PNearVec ) .dot( Summate( WVec , NNearVec ) ) )
                   / ( Summate( WVec , PNearVec , PNearVec ) - Summate( WNVec , PNearVec ) .dot( Summate( WVec , PNearVec ) ) ) ;

    Imath::V3f ULinear = Summate( WNVec , NNearVec ) - 2.0f * UQuadric * Summate( WNVec , PNearVec ) ;

    float UConstant = - ULinear.dot( Summate( WNVec , PNearVec ) ) - UQuadric * Summate( WNVec , PNearVec , PNearVec ) ;

    //
    U[0] = UConstant ;
    U[1] = ULinear.x ;
    U[2] = ULinear.y ;
    U[3] = ULinear.z ;
    U[4] = UQuadric ;

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
    const float H = 1.0f ;

    //
    std::vector< Imath::V3f > PVec ;
    std::vector< Imath::V3f > NVec ;

    GenerateSamples( PVec , NVec , NumSamples , Center , R ) ;

    //
    const Imath::V3f & X = PVec.front() ;

    float U[5] = { 0.0f , 0.0f , 0.0f , 0.0f , 0.0f } ;
    float UHat[5] = { 0.0f , 0.0f , 0.0f , 0.0f , 0.0f } ;
    SolveU( U , UHat , X , 5.0f , 10.0f , H , PVec , NVec ) ;

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
    const float CellSize = 0.5f ;

    //
    const int NumSamples = 500 ;

    std::vector< Imath::V3f > PVec ; PVec.reserve( NumSamples ) ;
    for ( int i = 0 ; i < NumSamples ; ++ i )
    {
        Imath::V3f P( static_cast< float >( rand() ) / RAND_MAX *  1.0f ,
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
    std::vector< unsigned long > KNNIdxVec ; KNNIdxVec.reserve( 1000 ) ;
    std::vector< double > DistVec ; DistVec.reserve( 1000 ) ;

    sfcnn< Imath::V3f , 3 , float > ANN( & PVec[0] , NumSamples , 4 ) ;
    

    QueryPerformanceFrequency(&Frequency); 
    QueryPerformanceCounter(&StartingTime);
    ANN.ksearch( X , 100 , KNNIdxVec , DistVec ) ;

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

void ReadPointCloud( Imath::Box3f & Bound , std::vector< Imath::V3f > & PVec , std::vector< Imath::V3f > & NVec , float & RAverage , const char * FilePath )
{
    PtcPointCloud PTCHandle = PtcSafeOpenPointCloudFile( FilePath ) ;
    if ( ! PTCHandle )
    {
        exit( EXIT_FAILURE ) ;
    }

    int NumPoints = 0 ;
    PtcGetPointCloudInfo( PTCHandle , "npoints" , & NumPoints ) ;

    PtcGetPointCloudInfo( PTCHandle , "bbox" , & Bound.min.x ) ;

    PVec.reserve( NumPoints ) ;
    NVec.reserve( NumPoints ) ;
    RAverage = 0.0f ;
    float PerPointData[32];

    for ( int i = 0 ; i < NumPoints ; ++ i )
    {
        Imath::V3f P( 0.0f ) , N( 0.0f ) ;
        float R = - 1.0f ;
        PtcReadDataPoint( PTCHandle , & P.x , & N.x , & R , PerPointData ) ;

        PVec.push_back( P ) ;
        NVec.push_back( N ) ;

        RAverage += R ;
    }
    RAverage /= NumPoints ;

    PtcClosePointCloudFile( PTCHandle ) , PTCHandle = NULL ;
}

int main( int Argc , char * Argv[] )
{
    -- Argc , ++ Argv ;
    if ( Argc != 2 )
    {
        return EXIT_FAILURE ;
    }

    //
    const char * PTCFilePath = Argv[0] ;
    const float H = static_cast< float >( atof( Argv[1] ) ) ;

    //
    Imath::Box3f Bound ;
    std::vector< Imath::V3f > PVec ;
    std::vector< Imath::V3f > NVec ;
    float RAverage = - 1.0f ;
    ReadPointCloud( Bound , PVec , NVec , RAverage , PTCFilePath ) ;

    const float R = RAverage * 10 ;

    //
    sfcnn< Imath::V3f , 3 , float > ANN( & PVec[0] , PVec.size() ) ;

    const int NumThreads = 1 ;
    omp_set_num_threads( NumThreads ) ;

    std::vector< unsigned long > LocalIVecPool[NumThreads] ;
    std::vector< double > LocalDVecPool[NumThreads] ;
    std::vector< Imath::V3f > LocalPVecPool[NumThreads] ;
    std::vector< Imath::V3f > LocalNVecPool[NumThreads];

    for ( int i = 0 ; i < NumThreads ; ++ i )
    {
        LocalIVecPool[i].reserve( 1024 ) ;
        LocalDVecPool[i].reserve( 1024 ) ;
        LocalPVecPool[i].reserve( 1024 ) ;
        LocalNVecPool[i].reserve( 1024 ) ;
    }

    #pragma omp parallel
    {
        #pragma omp for
        for ( int i = 0 ; i < static_cast< int >( PVec.size() ) ; ++ i )
        {
            const Imath::V3f & X = PVec[i];

            //
            const int ThreadId = omp_get_thread_num() ;

            std::vector< unsigned long > & LocalIVec = LocalIVecPool[ThreadId] ;
            LocalIVec.clear() ;

            std::vector< double > & LocalDVec = LocalDVecPool[ThreadId] ;
            LocalDVec.clear() ;

            //
            int NumLocalSamples = 100 ;
            ANN.ksearch( X , NumLocalSamples , LocalIVec , LocalDVec ) ;

            while ( sqrt( LocalDVec.back() ) < R )
            {
                NumLocalSamples *= 1.5 ;

                LocalIVec.clear() ;
                LocalDVec.clear() ;
                ANN.ksearch( X , NumLocalSamples , LocalIVec , LocalDVec ) ;
            }

            //
            std::vector< Imath::V3f > & LocalPVec = LocalPVecPool[ThreadId] ;
            LocalPVec.clear() ;

            std::vector< Imath::V3f > & LocalNVec = LocalNVecPool[ThreadId] ;
            LocalNVec.clear() ;

            size_t NumNeighbors = std::upper_bound( LocalDVec.begin() , LocalDVec.end() , static_cast< double >( R * R ) , std::less< double >() ) - LocalDVec.begin() ;
            for ( std::vector< unsigned long >::const_iterator itr = LocalIVec.begin() ; itr != ( LocalIVec.begin() + NumNeighbors ) ; ++ itr )
            {
                LocalPVec.push_back( PVec[* itr] ) ;
                LocalNVec.push_back( NVec[* itr] ) ;
            }

            float U[5] = { 0.0f , 0.0f , 0.0f , 0.0f , 0.0f } ;
            float UHat[5] = { 0.0f , 0.0f , 0.0f , 0.0f , 0.0f } ;
            SolveU( U , UHat , X , R , RAverage , H , LocalPVec , LocalNVec ) ;

            //
            float S = U[0] * 1.0f +
                      U[1] * X.x +
                      U[2] * X.y +
                      U[3] * X.z +
                      U[4] * X.length2() ;
            float Ka = 2.0f * UHat[4] ;
            printf( "S(x) = %f, Ka = %f\n" , S , Ka ) ;
        }
    }

    return EXIT_SUCCESS ;
}
