// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/** @file lib/Averager.h
 *
 * @brief C++ manager for spatial averaging of an etree database
 * containing the USGS central CA velocity model.
 *
 * This C++ code is based on the C convertdb application written by
 * Julio Lopez (Carnegie Mellon University).
 */

#if !defined(cencalvm_average_averager_h)
#define cencalvm_average_averager_h

#include <string> // HASA std::string
#include "cencalvm/storage/etreefwd.h" // HOLDSA etree_t

namespace cencalvm {
  namespace average {
    class Averager;
  } // namespace average
  namespace storage {
    class ErrorHandler; // HOLDSA ErrorHandler
  } // namespace storage
} // namespace cencalvm

/// C++ manager for spatial averaging of an etree database containing
/// the USGS central CA velocity model
class cencalvm::average::Averager
{ // Projector

public :
  // PUBLIC METHODS /////////////////////////////////////////////////////

  /// Constructor.
  Averager(void);

  /// Destructor
  ~Averager(void);

  /** Set filename of input database.
   *
   * @param filename Name of file
   */
  void filenameIn(const char* filename);

  /** Set filename of output database.
   *
   * @param filename Name of file
   */
  void filenameOut(const char* filename);

  /** Spatially average etree database by filling in etree octants
   * with average of their children.
   */
  void average(void);

  /** Get handle to error handler.
   *
   * @returns Pointer to Error handler
   */
  cencalvm::storage::ErrorHandler* errorHandler(void);

  /** Set flag indicating creation should be quiet (no progress reports).
   *
   * Default behavior is for to give progress reports.
   *
   * @param flag True for quiet operation, false to give progress reports
   */
  void quiet(const bool flag);

private :
  // PRIVATE METHODS ////////////////////////////////////////////////////

  Averager(const Averager& p); ///< Not implemented
  const Averager& operator=(const Averager& p); ///< Not implemented
  
private :
  // PRIVATE METHODS ////////////////////////////////////////////////////

private :
  // PRIVATE MEMBERS ////////////////////////////////////////////////////

  etree_t* _dbIn; ///< Input database
  etree_t* _dbAvg; ///< Database with averaging

  std::string _filenameIn; ///< Filename of input database
  std::string _filenameOut; ///< Filename of output database
  
  cencalvm::storage::ErrorHandler* _pErrHandler; ///< Error handler

  bool _quiet; ///< Flag to eliminate progress reports

}; // Averager

#include "Averager.icc" // inline methods

#endif // cencalvm_average_averager_h

// version
// $Id: Averager.h,v 1.4 2005/10/11 18:24:22 brad Exp $

// End of file 
