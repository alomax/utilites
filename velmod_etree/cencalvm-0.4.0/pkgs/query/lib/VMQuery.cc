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

#include "VMQuery.h" // implementation of class methods

#include "cencalvm/storage/Payload.h" // USES PayloadStruct
#include "cencalvm/storage/Geometry.h" // USES Geometry
#include "cencalvm/storage/ErrorHandler.h" // USES ErrorHandler

extern "C" {
#include "etree.h"
}

#include <sstream> // USES std::ostringstream
#include <iomanip> // USES setw(), setiosflags(), resetiosflags()
#include <strings.h> // USES strcasecmp()
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
/// Default constructor
cencalvm::query::VMQuery::VMQuery(void) :
  _queryRes(0),
  _db(0),
  _dbExt(0),
  _filename(""),
  _filenameExt(""),
  _pQueryVals(0),
  _querySize(0),
  _cacheSize(128),
  _cacheSizeExt(128),
  _queryFn(&cencalvm::query::VMQuery::_queryMax),
  _pErrHandler(new cencalvm::storage::ErrorHandler),
  _pGeom(new cencalvm::storage::Geometry(*_pErrHandler))
{ // constructor
  const int querySize = 8;
  _pQueryVals = (querySize > 0) ? new int[querySize] : 0;
  for (int i=0; i < querySize; ++i)
    _pQueryVals[i] = i;
  _querySize = querySize;
} // constructor
  
// ----------------------------------------------------------------------
// Default destructor.
cencalvm::query::VMQuery::~VMQuery(void)
{ // destructor
  delete[] _pQueryVals; _pQueryVals = 0;
  delete _pGeom; _pGeom = 0;
  delete _pErrHandler; _pErrHandler = 0;
  close();
} // destructor

// ----------------------------------------------------------------------
// Open the database and prepare for querying.
void
cencalvm::query::VMQuery::open(void)
{ // open
  if (0 == _db) { // database is not already open
    assert(_cacheSize > 0);
    _db = etree_open(_filename.c_str(), O_RDONLY, _cacheSize, 0, 0);
    if (0 == _db) {
      std::ostringstream msg;
      msg << "Could not open the etree database '" << _filename
	  << "' for querying.";
      _pErrHandler->error(msg.str().c_str());
    } // if
  } // if

  if (0 != strcmp(_filenameExt.c_str(), "") && 0 == _dbExt) {
    assert(_cacheSizeExt > 0);
    _dbExt = etree_open(_filenameExt.c_str(), O_RDONLY, _cacheSizeExt, 0, 0);
    if (0 == _dbExt) {
      std::ostringstream msg;
      msg << "Could not open the etree (regional) database '"
	  << _filenameExt << "' for querying.";
      _pErrHandler->error(msg.str().c_str());
    } // if
  } // if
} // open
  
// ----------------------------------------------------------------------
// Close the database.
void
cencalvm::query::VMQuery::close(void)
{ // close
  if (0 != _db && 0 != etree_close(_db)) {
    std::ostringstream msg;
    msg << "Could not close the etree database '" << _filename << "'.";
    _pErrHandler->error(msg.str().c_str());
  } // if
  _db = 0;

  if (0 != _dbExt && 0 != etree_close(_dbExt)) {
    std::ostringstream msg;
    msg << "Could not close the etree database '" << _filename << "'.";
    _pErrHandler->error(msg.str().c_str());
  } // if
  _dbExt = 0;
} // close
  
// ----------------------------------------------------------------------
// Set query type.
void
cencalvm::query::VMQuery::queryType(const QueryEnum queryType)
{ // queryType
  switch (queryType)
    { // switch
    case MAXRES :
      _queryFn = &cencalvm::query::VMQuery::_queryMax;
      break;
    case FIXEDRES :
      _queryFn = &cencalvm::query::VMQuery::_queryFixed;
      break;
    case WAVERES :
      _queryFn = &cencalvm::query::VMQuery::_queryWave;
      break;
    default :
      _pErrHandler->error("Could not find query function for requested "
			 "query type.");
    } // switch
} // queryType

// ----------------------------------------------------------------------
// Set values to be returned by queries.
void
cencalvm::query::VMQuery::queryVals(const char** names,
				    const int numVals)
{ // queryVals
  assert(0 != names);

  if (numVals > 0) {
    delete[] _pQueryVals; 
    _pQueryVals = (numVals > 0) ? new int[numVals] : 0;
    _querySize = numVals;
    for (int iVal=0; iVal < numVals; ++iVal) {
      if (0 == strcasecmp("Vp", names[iVal]))
	_pQueryVals[iVal] = 0;
      else if (0 == strcasecmp("Vs", names[iVal]))
	_pQueryVals[iVal] = 1;
      else if (0 == strcasecmp("Density", names[iVal]))
	_pQueryVals[iVal] = 2;
      else if (0 == strcasecmp("Qp", names[iVal]))
	_pQueryVals[iVal] = 3;
      else if (0 == strcasecmp("Qs", names[iVal]))
	_pQueryVals[iVal] = 4;
      else if (0 == strcasecmp("DepthFreeSurf", names[iVal]))
	_pQueryVals[iVal] = 5;
      else if (0 == strcasecmp("FaultBlock", names[iVal]))
	_pQueryVals[iVal] = 6;
      else if (0 == strcasecmp("Zone", names[iVal]))
	_pQueryVals[iVal] = 7;
      else {
	std::ostringstream msg;
	msg << "Value name '" << names[iVal] << "' does not match any "
	    << "of the values in the velocity database.\n"
	    << "Values in database are:\n"
	    << "  Vp, Vs, Density, Qp, Qs, DepthFreeSurf, FaultBlock, Zone.";
	_pErrHandler->error(msg.str().c_str());
      } // else
    } // for
  } // if
} // queryVals

// ----------------------------------------------------------------------
// Query the database.
void
cencalvm::query::VMQuery::query(double** ppVals,
				const int numVals,
				const double lon,
				const double lat,
				const double elev)
{ // query
  assert(0 != _queryFn);
  assert(0 != ppVals);
  assert(numVals == _querySize);
  assert(0 != _pGeom);

  bool useAddr = false;
  etree_addr_t addr;
  cencalvm::storage::PayloadStruct payload;
  (this->*_queryFn)(&payload, &addr, _db, lon, lat, elev, useAddr);

  // If not found in detailed model, query the regional model
  char buf[ETREE_MAXBUF];
  if (cencalvm::storage::Payload::NODATABLOCK == payload.FaultBlock &&
      0 != _dbExt) {
    useAddr = true;
    (this->*_queryFn)(&payload, &addr, _dbExt, lon, lat, elev, useAddr);
  } // if

  // If not found in any model, set data to NODATA values
  if (cencalvm::storage::Payload::NODATABLOCK == payload.FaultBlock) {
    std::ostringstream msg;
    msg
      << std::resetiosflags(std::ios::fixed)
      << std::setiosflags(std::ios::scientific)
      << std::setprecision(6)
      << lon << ", " << lat << ", " << elev << ", No data\n";
    _pErrHandler->log(msg.str().c_str());
    std::ostringstream warning;
    warning
      << "WARNING: No data for "
      << std::resetiosflags(std::ios::fixed)
      << std::setiosflags(std::ios::scientific)
      << std::setprecision(6)
      << lon << ", " << lat << ", " << elev << ".\n";
    _pErrHandler->warning(warning.str().c_str());
  } // if

  // Copy values from payload into array
  for (int i=0; i < _querySize; ++i) {
    switch (_pQueryVals[i])
      { // switch
      case 0 :
	(*ppVals)[i] = payload.Vp;
	break;
      case 1 :
	(*ppVals)[i] = payload.Vs;
	break;
      case 2 :
	(*ppVals)[i] = payload.Density;
	break;
      case 3 :
	(*ppVals)[i] = payload.Qp;
	break;
      case 4 :
	(*ppVals)[i] = payload.Qs;
	break;
      case 5 :
	(*ppVals)[i] = payload.DepthFreeSurf;
	break;
      case 6 :
	(*ppVals)[i] = payload.FaultBlock;
	break;
      case 7 :
	(*ppVals)[i] = payload.Zone;
	break;
      default :
	_pErrHandler->error("Could not parse requested query value.");
      } // switch
  } // for
} // query

// ----------------------------------------------------------------------
// Query database at maximum resolution possible.
void
cencalvm::query::VMQuery::_queryMax(cencalvm::storage::PayloadStruct* pPayload,
				    etree_addr_t* pAddr,
				    etree_t* pDB,
				    const double lon,
				    const double lat,
				    const double elev,
				    const bool useAddr)
{ // _queryMax
  assert(0 != pPayload);
  assert(0 != pDB);
  assert(0 != pAddr);

  if (!useAddr) {
    pAddr->level = ETREE_MAXLEVEL;
    pAddr->type = ETREE_LEAF;
    _pGeom->lonLatElevToAddr(pAddr, lon, lat, elev);
    if (_pErrHandler->status() != cencalvm::storage::ErrorHandler::OK)
      return;
  } // if

  etree_addr_t resAddr;
  int err = etree_search(pDB, *pAddr, &resAddr, "*", pPayload);
  // If search returned interior octant (averaged), return no data
  // instead of averaged values since query request is for maximum
  // resolution and we don't have a leaf octant (data) at that
  // location.
  if (0 != err ||
      ETREE_INTERIOR == resAddr.type)
    _setNoData(pPayload);
} // _queryMax

// ----------------------------------------------------------------------
// Query database at fixed resolution. Resolution is specified by queryRes().
void
cencalvm::query::VMQuery::_queryFixed(cencalvm::storage::PayloadStruct* pPayload,
				    etree_addr_t* pAddr,
				    etree_t* pDB,
				    const double lon,
				    const double lat,
				    const double elev,
				    const bool useAddr)
{ // _queryFixed
  assert(0 != pPayload);
  assert(0 != pDB);
  assert(0 != pAddr);
  assert(0 != _pGeom);

  if (!useAddr) {
    const double vertExag = cencalvm::storage::Geometry::vertExag();
    pAddr->level = _pGeom->level(vertExag * _queryRes);
    pAddr->type = ETREE_LEAF;
    _pGeom->lonLatElevToAddr(pAddr, lon, lat, elev);
    if (_pErrHandler->status() != cencalvm::storage::ErrorHandler::OK)
      return;
  } // if

  etree_addr_t resAddr;
  const int err = etree_search(pDB, *pAddr, &resAddr, "*", pPayload);
  // if search returned interior octant at coarser resolution than
  // what we want, return no data instead of averaged octant since
  // query request was for a given resolution and we don't have a leaf
  // octant at that resolution or one higher.
  if (0 != err ||
      (ETREE_INTERIOR == resAddr.type && pAddr->level > resAddr.level))
    _setNoData(pPayload);
} // _queryFixed

// ----------------------------------------------------------------------
// Query database at resolution specified by wavelength. Resolution
// is specified by queryRes().
void
cencalvm::query::VMQuery::_queryWave(cencalvm::storage::PayloadStruct* pPayload,
				    etree_addr_t* pAddr,
				    etree_t* pDB,
				    const double lon,
				    const double lat,
				    const double elev,
				    const bool useAddr)
{ // _queryWave
  assert(0 != pPayload);
  assert(0 != pDB);
  assert(0 != pAddr);

  if (!useAddr) {
    pAddr->level = ETREE_MAXLEVEL;
    pAddr->type = ETREE_LEAF;
    _pGeom->lonLatElevToAddr(pAddr, lon, lat, elev);
    if (_pErrHandler->status() != cencalvm::storage::ErrorHandler::OK)
      return;
  } // if

  etree_addr_t resAddr;
  const int err = etree_search(pDB, *pAddr, &resAddr, "*", pPayload);

  const double vertExag = cencalvm::storage::Geometry::vertExag();
  const double minPeriod = vertExag * _queryRes;

  // If search returned interior octant (averaged), return no data
  // instead of averaged values if interior octant is coarser than we want
  if (0 != err ||
      (ETREE_INTERIOR == resAddr.type &&
       _pGeom->edgeLen(resAddr.level) / pPayload->Vs > minPeriod)) {
    _setNoData(pPayload);
    return;
  } // if
  
  cencalvm::storage::PayloadStruct childPayload = *pPayload;
  while (pPayload->Vs > 0.0 && 
	 _pGeom->edgeLen(resAddr.level) / pPayload->Vs < minPeriod &&
	 resAddr.level > 0) {
    childPayload = *pPayload;
    etree_addr_t parentAddr;
    _pGeom->findAncestor(&parentAddr, resAddr, resAddr.level-1);
    if (0 != etree_search(pDB, parentAddr, &resAddr, "*", pPayload)) {
      char buf[ETREE_MAXBUF];
      std::ostringstream msg;
      msg
	<< "Could not find parent octant " << 
	etree_straddr(pDB, buf, parentAddr)
	<< "\nof child octant " << etree_straddr(pDB, buf, resAddr) << "\n"
	<< "for location " << lon << ", " << lat << ", " << elev
	<< ".\nUsing values from child octant.";
      _pErrHandler->warning(msg.str().c_str());
      etree_search(pDB, *pAddr, &resAddr, "*", pPayload);
      return;
    } // if
  } // while
  
  // Values returned are averages between values in payload that first
  // meets the wavelength criteria and the child that didn't. This
  // helps reduce the bias towards faster wave speeds. Fault block
  // values are taken from parent which is already in pPayload.
  pPayload->Vp = 0.5 * (childPayload.Vp + pPayload->Vp);
  pPayload->Vs = 0.5 * (childPayload.Vs + pPayload->Vs);
  pPayload->Density = 0.5 * (childPayload.Density + pPayload->Density);
  pPayload->Qp = 0.5 * (childPayload.Qp + pPayload->Qp);
  pPayload->Qs = 0.5 * (childPayload.Qs + pPayload->Qs);
  pPayload->DepthFreeSurf = 
    0.5 * (childPayload.DepthFreeSurf + pPayload->DepthFreeSurf);
} // _queryWave

// ----------------------------------------------------------------------
// Set payload to NODATA values.
void
cencalvm::query::VMQuery::_setNoData(cencalvm::storage::PayloadStruct* 
				     pPayload)
{ // _setNoData
  pPayload->Vp = cencalvm::storage::Payload::NODATAVAL;
  pPayload->Vs = cencalvm::storage::Payload::NODATAVAL;
  pPayload->Density = cencalvm::storage::Payload::NODATAVAL;
  pPayload->Qp = cencalvm::storage::Payload::NODATAVAL;
  pPayload->Qs = cencalvm::storage::Payload::NODATAVAL;
  pPayload->DepthFreeSurf = cencalvm::storage::Payload::NODATAVAL;
  pPayload->FaultBlock = cencalvm::storage::Payload::NODATABLOCK;
  pPayload->Zone = cencalvm::storage::Payload::NODATAZONE;
} // _setNoData

// version
// $Id: VMQuery.cc,v 1.16 2005/10/29 05:54:02 brad Exp $

// End of file 
