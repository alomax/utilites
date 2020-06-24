// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file tests/TestGeometry.h
 *
 * @brief C++ TestGeometry object
 *
 * C++ unit testing for TestGeometry.
 */

#if !defined(cencalvm_storage_testgeometry_h)
#define cencalvm_storage_testgeometry_h

#include <cppunit/extensions/HelperMacros.h>

namespace cencalvm {
  namespace storage {
    class TestGeometry;
  } // storage
} // cencalvm

/// C++ unit testing for Geometry
class cencalvm::storage::TestGeometry : public CppUnit::TestFixture
{ // class TestGeometry

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestGeometry );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testAddress );
  CPPUNIT_TEST( testInvAddress );
  CPPUNIT_TEST( testEdgeLen );
  CPPUNIT_TEST( testLevel );
  CPPUNIT_TEST( testFindAncestor );
  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor
  void testConstructor(void);

  /// Test lonLatElevToAddr()
  void testAddress(void);

  /// Test addrToLonLatElev()
  void testInvAddress(void);

  /// Test edgeLen()
  void testEdgeLen(void);

  /// Test level()
  void testLevel(void);

  /// Test findAncestor()
  void testFindAncestor(void);

}; // class TestGeometry

#endif // cencalvm_storage_testgeometry

// version
// $Id: TestGeometry.h,v 1.6 2005/10/11 18:28:21 brad Exp $

// End of file 
