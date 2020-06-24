// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// <LicenseText>
//
// ----------------------------------------------------------------------
//

/** @file extensions/cencalvmdb/lib/CenCalVMDB.h
 *
 * @brief C++ CenCalVMDB.
 *
 * Wrap VMQuery to create a SpatialDB object.
 */

#if !defined(cencalvm_cencalvmdb_cencalvmdb_h)
#define cencalvm_cencalvmdb_cencalvmdb_h

#include "cencalvm/query/VMQuery.h" // USES VMQuery

namespace cencalvm {
  namespace extensions {
    namespace cencalvmdb {
      class CenCalVMDB;
    } // cencalvmdb
  } // extensions
} // cencalvm

namespace spatialdata {
  namespace geocoords {
    class CSGeo; // forward declaration
  } // geocoords
} // spatialdata

/// C++ object for CenCalVMDB.
class cencalvm::extensions::cencalvmdb::CenCalVMDB : 
  public spatialdata::spatialdb::SpatialDB
{ // class CenCalVMDB

 public :
  // PUBLIC METHODS /////////////////////////////////////////////////////
  
  /// Default constructor.
  CenCalVMDB(void);
  
  /** Constructor with label.
   *
   * @param label Label for database
   */
  CenCalVMDB(const char* label);
  
  /// Default destructor.
  ~CenCalVMDB(void);

  /// Open the database and prepare for querying.
  void open(void);

  /// Close the database.
  void close(void);

  /** Set query type.
   *
   * @param queryType Set type of query
   */
  void queryType(const query::VMQuery::QueryEnum queryType);

  /** Set query resolution. Resolution is associated with vertical
   * resolution, which is 4 times greater than the horizontal
   * resolution.
   *
   * Meaning depends on type of query:
   *   @li MAXRES Resolution is not applicable
   *   @li FIXEDRES Query etree at level corresponding to resolution
   *   @li WAVERES Resolution corresponds to minimum period of waves;
   *     etree is queried at level corresponding to minimum wavelength
   *     for shear waves.
   *
   * @param res Resolution of query.
   */
  void queryRes(const double res);

  /** Set minimum shear-wave speed (clip values lower than this).
   *
   * @param vs Shear wave speed
   */
  void minVs(const double vs);
  
  /** Set values to be returned by queries.
   *
   * @pre Must call open() before queryVals()
   *
   * @param names Names of values to be returned in queries
   * @param numVals Number of values to be returned in queries
   */
  void queryVals(const char** names,
		 const int numVals);

  /** Set the database filename.
   *
   * @param filename Name of database file
   */
  void filename(const char* filename);

  /** Set size of cache during queries.
   *
   * @param size Size of cache in MB
   */
  void cacheSize(const int size);

  /** Set the extended database filename.
   *
   * @param filename Name of extended database file
   */
  void filenameExt(const char* filename);

  /** Set size of cache during queries of extended database
   *
   * @param size Size of cache in MB
   */
  void cacheSizeExt(const int size);

  /** Query the database.
   *
   * @pre Must call open() before query()
   *
   * @param pVals Pointer to computed values (output from query)
   * @param numVals Number of values expected (size of pVals array)
   * @param x X coordinate of location for query
   * @param y Y coordinate of location for query
   * @param z Z coordinate of location for query
   * @param pCSQuery Coordinate system of coordinates
   */
  void query(double** pVals,
	     const int numVals,
	     const double x,
	     const double y,
	     const double z,
	     const spatialdata::geocoords::CoordSys* pCSQuery);

 private :
  // PRIVATE METHODS ////////////////////////////////////////////////////
  
  CenCalVMDB(const CenCalVMDB& data); ///< Not implemented
  const CenCalVMDB& operator=(const CenCalVMDB& data); ///< Not implemented
  
private :
 // PRIVATE MEMBERS ////////////////////////////////////////////////////
  
  double _minVs; ///< Minimum shear wave speed
  int _vsVal; ///< Index of query value corresponding to Vs

  cencalvm::query::VMQuery* _pQuery; ///< Pointer to velocity model query
  spatialdata::geocoords::CSGeo* _pCS; ///< Pointer to coord system of VMQuery

}; // class CenCalVMDB

#include "CenCalVMDB.icc" // inline methods

#endif // cencalvm_cencalvmdb_cencalvmdb_h

// version
// $Id: CenCalVMDB.h,v 1.3 2005/11/01 02:36:22 brad Exp $

// End of file 
