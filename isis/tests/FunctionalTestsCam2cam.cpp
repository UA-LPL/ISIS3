#include <iostream>
#include <QTemporaryFile>

#include "cam2cam.h"

#include "Cube.h"
#include "CubeAttribute.h"
#include "IException.h"
#include "PixelType.h"
#include "Pvl.h"
#include "PvlGroup.h"
#include "PvlKeyword.h"
#include "TestUtilities.h"
#include "FileName.h"
#include "ProjectionFactory.h"
#include "CameraFixtures.h"
#include "Mocks.h"

using namespace Isis;
using ::testing::Return;
using ::testing::AtLeast;

static QString APP_XML = FileName("$ISISROOT/bin/xml/cam2cam.xml").expanded();

TEST_F(DefaultCube, FunctionalTestCam2CamNoChange) {

  QVector<QString> args = {"to="+tempDir.path()+"/Cam2CamNoChange.cub", "INTERP=BILINEAR"};
  UserInterface ui(APP_XML, args);

  testCube->reopen("r");
  QString inFile = testCube->fileName();
  Cube mcube(inFile,"r");

  cam2cam(testCube, &mcube, ui);

  Cube icube(inFile);
  PvlGroup icubeInstrumentGroup = icube.label()->findGroup("Instrument", Pvl::Traverse);

  Cube ocube(tempDir.path()+"/Cam2CamNoChange.cub");
  PvlGroup ocubeInstrumentGroup = ocube.label()->findGroup("Instrument", Pvl::Traverse);

  ASSERT_EQ(icubeInstrumentGroup.findKeyword("SpacecraftName"), ocubeInstrumentGroup.findKeyword("SpacecraftName"));
  ASSERT_EQ(icubeInstrumentGroup.findKeyword("InstrumentId"), ocubeInstrumentGroup.findKeyword("InstrumentID"));
  ASSERT_EQ(icubeInstrumentGroup.findKeyword("TargetName"), ocubeInstrumentGroup.findKeyword("TargetName"));
  ASSERT_EQ(icubeInstrumentGroup.findKeyword("StartTime"), ocubeInstrumentGroup.findKeyword("StartTime"));
  ASSERT_EQ(icubeInstrumentGroup.findKeyword("ExposureDuration"), ocubeInstrumentGroup.findKeyword("ExposureDuration"));
  ASSERT_EQ(icubeInstrumentGroup.findKeyword("SpacecraftClockCount"), ocubeInstrumentGroup.findKeyword("SpacecraftClockCount"));
  ASSERT_EQ(icubeInstrumentGroup.findKeyword("FloodModeId"), ocubeInstrumentGroup.findKeyword("FloodModeId"));
  ASSERT_EQ(icubeInstrumentGroup.findKeyword("GainModeId"), ocubeInstrumentGroup.findKeyword("GainModeId"));
  ASSERT_EQ(icubeInstrumentGroup.findKeyword("OffsetModeId"), ocubeInstrumentGroup.findKeyword("OffsetModeId"));
}

// This test evaluates the behavior of example 1 in the cam2cam app. The pixel
// locations are arbitrary but in the regions described in the example
// explanation. Unfortunately, the OREX OCAMS ALE implementation is broken
// in that the lat/lons are way off. However, the mapping behavior seems 
// consistent as described...maybe. ALE is preferred since it does not require
// repo storage of a complete or subsampled binary ISIS image.
TEST_F( OrexManyIsdCameraCubes, FunctionalTestCam2CamOffbody ) {
  typedef std::unique_ptr<Cube>   LocalCubePtr;

  const double tolerance_d = 0.05;
  const double tolerance_p = 0.00001;

  const bool OffBodyTrue      = true;
  const bool OffBodyFalse     = false;
  const bool OffBodyTrimTrue  = true;
  const bool OffBodyTrimFalse = false;

  QString fromcube( "cam2cam_bennu_from.cub" );
  QString matchcube( "cam2cam_bennu_match.cub" );

  LocalCubePtr from( make_cube( fromcube ) );
  LocalCubePtr match( make_cube( matchcube ) );

  Camera *from_c = from->camera();
  Camera *match_c = match->camera();

  double v_match_sample;
  double v_match_line;
  double v_from_sample;
  double v_from_line;  

  double v_latitude;
  double v_longitude;
  
  //----------------------------------------------------------------
  // See how well pixel locations map to lat/lon and back where
  // both cameras use the DSK (by default).
  // Check few mapping points
  v_match_sample = 512.0;
  v_match_line   = 512.0;
  EXPECT_TRUE( match_c->SetImage( v_match_sample, v_match_line ) );
  v_latitude  = match_c->UniversalLatitude( );
  v_longitude = match_c->UniversalLongitude( );

  // ALE OREX/OCAMS geometry is severely inaccurate, but behavior seems consistent.
  // CORRECT VALUES
  // EXPECT_NEAR( v_latitude,  -31.317275598917,  tolerance_d );  
  // EXPECT_NEAR( v_longitude, 329.17715075029,   tolerance_d );

  // ALE INCORRECT VALUES
  EXPECT_NEAR( v_latitude,  -27.297110494335211,  tolerance_d );
  EXPECT_NEAR( v_longitude, 332.68030130833739,   tolerance_d );  

  EXPECT_TRUE( from_c->SetUniversalGround( v_latitude, v_longitude) );
  double v_insample = from_c->Sample();
  double v_inline   = from_c->Line();
  EXPECT_NEAR( v_insample, v_match_sample, tolerance_p );
  EXPECT_NEAR( v_inline,   v_match_line,   tolerance_p );

  // ALE OREX/OCAMS geometry is severely inaccurate, but behavior seems consistent.
  // EXPECT_NEAR( match_c->RightAscension( ),     179.20061977189,   tolerance_d );
  // EXPECT_NEAR( match_c->Declination( ),         -1.9093150733987, tolerance_d );

  EXPECT_NEAR( match_c->RightAscension( ),     179.38309125410601,   tolerance_d );
  EXPECT_NEAR( match_c->Declination( ),         -1.9093150733987, tolerance_d );

  // Check ellipsoid vs DSK mapping...
  from_c->IgnoreElevationModel( false );  // Turns on shape model
  match_c->IgnoreElevationModel( true );  // Uses PCK ellipsoid

  v_match_sample = 365.0;
  v_match_line   =  50.0;
  EXPECT_TRUE( from_c->SetImage( v_match_sample, v_match_line ) );
  EXPECT_TRUE( match_c->SetImage( v_match_sample, v_match_line ) );

  EXPECT_TRUE( match_c->UniversalLatitude( )  != from_c->UniversalLatitude() );
  EXPECT_TRUE( match_c->UniversalLongitude( ) != from_c->UniversalLongitude() );

  // Reverse ordering and check results...
  from_c->IgnoreElevationModel( true );    // Turns off shape model
  match_c->IgnoreElevationModel( false );  // Uses DSK ellipsoid

  EXPECT_TRUE( from_c->SetImage( v_match_sample, v_match_line ) );
  EXPECT_TRUE( match_c->SetImage( v_match_sample, v_match_line ) );

  EXPECT_TRUE( match_c->UniversalLatitude( )  != from_c->UniversalLatitude() );
  EXPECT_TRUE( match_c->UniversalLongitude( ) != from_c->UniversalLongitude() );

  // Off lower left limb of DSK...
  v_match_sample = 300.0;
  v_match_line   = 970.0;
     
  EXPECT_TRUE( from_c->SetImage( v_match_sample, v_match_line ) );
  EXPECT_FALSE( match_c->SetImage( v_match_sample, v_match_line ) );

  //----------------------------------------------------------------
  // Create indivudual mapper functors
  Cam2CamMapper map_from_default(    from.get(), match.get(), OffBodyFalse, OffBodyTrimFalse  ); // 2nd image column
  Cam2CamMapper map_from_off_notrim( from.get(), match.get(), OffBodyTrue,  OffBodyTrimFalse );  // 3th image column
  Cam2CamMapper map_from_off_trim(   from.get(), match.get(), OffBodyTrue,  OffBodyTrimTrue );   // 4rd image column

  // Validate functor configurations
  EXPECT_TRUE( map_from_default.IsValid() );
  EXPECT_TRUE( map_from_off_trim.IsValid() );
  EXPECT_TRUE( map_from_off_notrim.IsValid() );

  // Tests clear off body in upper left corner
  v_match_sample = 50.0;
  v_match_line   = 50.0;

  EXPECT_FALSE( from_c->SetImage( v_match_sample, v_match_line) );
  EXPECT_FALSE( match_c->SetImage( v_match_sample, v_match_line) );

  EXPECT_FALSE( map_from_default.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) ); 
  EXPECT_TRUE( map_from_off_notrim.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) );
  EXPECT_TRUE( map_from_off_trim.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) ); 

  //----------------------------------------------------------------
  // Row A->D uses FROM w/DSK and MATCH w/ellipsoid
  // Yellow (NULL) section in lower left in panels 
  v_match_sample = 300.0;
  v_match_line   = 950.0;

  // Config for top panel A->D
  from_c->IgnoreElevationModel( false );  // Turns on shape model
  match_c->IgnoreElevationModel( true );  // Uses PCK ellipsoid

  // Intercepts Bennu with ellipsoid, not DSK
  EXPECT_FALSE( from_c->SetImage( v_match_sample, v_match_line) );
  EXPECT_TRUE( match_c->SetImage( v_match_sample, v_match_line) );

  EXPECT_FALSE( map_from_default.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) );   // A->B
  EXPECT_FALSE( map_from_off_notrim.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) ); // A->C
  EXPECT_FALSE( map_from_off_trim.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) );  // A->D

  // NULL portion in images B,C and D doesn't register here because we are using the 
  // NAIF DSK toolkit for these tests to maintain efficiency (load/init times for Bullet
  // are inefficent). ISIS NAIF DSK currently does not test for occlusion as Bullet does.
  v_match_sample = 440.0;
  v_match_line   = 80.0;

  // Config for top panel
  from_c->IgnoreElevationModel( false );
  match_c->IgnoreElevationModel( true );

  EXPECT_TRUE( from_c->SetImage( v_match_sample, v_match_line) );
  EXPECT_TRUE( match_c->SetImage( v_match_sample, v_match_line) );

  EXPECT_TRUE( map_from_default.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) );    // A->B
  EXPECT_TRUE( map_from_off_notrim.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) ); // A->C 
  EXPECT_TRUE( map_from_off_trim.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) );   // A->D


  //----------------------------------------------------------------
  // Row E->H uses FROM w/ellipsoid and MATCH w/DSK
  v_match_sample = 300.0;
  v_match_line   = 950.0;  

  from_c->IgnoreElevationModel( true );
  match_c->IgnoreElevationModel( false );

  EXPECT_TRUE( from_c->SetImage( v_match_sample, v_match_line) );
  EXPECT_FALSE( match_c->SetImage( v_match_sample, v_match_line) );

  EXPECT_TRUE( map_from_default.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) );     // E->F
  EXPECT_FALSE( map_from_off_notrim.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) ); // E->G
  EXPECT_FALSE( map_from_off_trim.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) );   // E->H


  // Test center of image. Should be good for all instances...
  v_match_sample = 512.0;
  v_match_line   = 512.0;

  from_c->IgnoreElevationModel( true );
  match_c->IgnoreElevationModel( false );

  EXPECT_TRUE( from_c->SetImage( v_match_sample, v_match_line) );
  EXPECT_TRUE( match_c->SetImage( v_match_sample, v_match_line) );

  EXPECT_NEAR( match_c->UniversalLatitude( ),  from_c->UniversalLatitude(),  tolerance_d );
  EXPECT_NEAR( match_c->UniversalLongitude( ), from_c->UniversalLongitude(),  tolerance_d );

  EXPECT_TRUE( map_from_default.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) );     // E->F
  EXPECT_TRUE( map_from_off_notrim.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) );  // E->G
  EXPECT_TRUE( map_from_off_trim.mapit( v_from_sample, v_from_line, v_match_sample, v_match_line) );    // E->H

}

