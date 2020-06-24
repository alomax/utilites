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

#include "GridIngester.h" // implementation of class methods

#include "cencalvm/storage/Payload.h" // USES PayloadStruct
#include "cencalvm/storage/Geometry.h" // USES Geometry
#include "cencalvm/storage/ErrorHandler.h" // USES ErrorHandler

extern "C" {
#include "etree.h"
}

#include <fstream> // USES std::ifstream
#include <iostream> // USES std::cout
#include <sstream> // USES std::ostringstream
#include <iomanip> // USES std::resetiosflags(), std::setprecision()
#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Add gridded data to database.
void
cencalvm::create::GridIngester::addGrid(etree_t** pDB,
					const char* filename,
				  cencalvm::storage::ErrorHandler& errHandler,
					const bool quiet)
{ // addGrid
  assert(0 != pDB);

  if (cencalvm::storage::ErrorHandler::ERROR != errHandler.status())
    errHandler.resetStatus();
  else
    return;

  std::ifstream fin(filename);
  if (!fin.is_open()) {
    std::ostringstream msg;
    msg << "Could not open grid file '" << filename
	<< "' for reading.";
    errHandler.error(msg.str().c_str());
    return;
  } // if

  cencalvm::storage::Geometry vmgeom(errHandler);

  if (!quiet)
    std::cout
      << "Beginning processing of '" << filename << "'..." << std::endl;

  double resHoriz = 0.0;
  double resVert = 0.0;
  int numX = 0;
  int numY = 0;
  int numZ = 0;
  int numTotal = 0;
  fin >> resHoriz >> resVert >> numX >> numY >> numZ >> numTotal;
  if (0.0 == resHoriz ||
      0.0 == resVert) {
    std::ostringstream msg;
    msg << "Could not read horizontal and vertical resolution in '"
	<< filename << "'.";
    errHandler.error(msg.str().c_str());
    return;
  } // if
  // convert resHoriz and resVert from km to m
  resHoriz *= 1.0e+3;
  resVert *= 1.0e+3;
  
  const double tolerance = 1.0e-6;
  const double vertExag = resHoriz / resVert;
  const double vertExagE = vmgeom.vertExag();
  if (fabs(1.0 - vertExag/vertExagE) > tolerance) {
    std::ostringstream msg;
    msg << "Vertical exaggeration of " << vertExag << " in '" << filename
	<< "' does not match velocity model vertical exaggeration of "
	<< vertExagE << ".";
    errHandler.error(msg.str().c_str());
    return;
  } // if

  const etree_tick_t level = vmgeom.level(resHoriz);
  const double edgeLen = vmgeom.edgeLen(level);
  if ( fabs(1.0 - resHoriz/edgeLen) > tolerance) {
    std::ostringstream msg;
    msg << "Horizontal resolution of " << resHoriz << " in '" << filename
	<< "' does not fit resolution of nearest level in database of "
	<< edgeLen << ".";
    errHandler.error(msg.str().c_str());
    return;
  } // if    

  int numAdded = 0;
  int numIgnored = 0;
  for (int i=0; i < numTotal; ++i) {
    double lon = 0.0;
    double lat = 0.0;
    double elev = 0.0;
    int volID = 0;
    cencalvm::storage::PayloadStruct payload;
    fin
      >> lon
      >> lat
      >> elev
      >> payload.Vp
      >> payload.Vs
      >> payload.Density
      >> payload.Qp
      >> payload.Qs
      >> payload.DepthFreeSurf
      >> payload.FaultBlock
      >> payload.Zone
      >> volID;
    if (!fin.good()) {
      errHandler.error("Couldn't parse line.");
      break;
    } // if
    if (payload.FaultBlock != cencalvm::storage::Payload::NODATABLOCK &&
	payload.Zone != cencalvm::storage::Payload::NODATAZONE) {
      // convert elev and depth from km to m
      if (elev != cencalvm::storage::Payload::NODATAVAL)
	elev *= 1.0e+3;
      if (payload.DepthFreeSurf != cencalvm::storage::Payload::NODATAVAL)
	payload.DepthFreeSurf *= 1.0e+3;
      
      // convert Vp & Vs from km/s to m/s
      if (payload.Vp != cencalvm::storage::Payload::NODATAVAL)
	payload.Vp *= 1.0e+3;
      if (payload.Vs != cencalvm::storage::Payload::NODATAVAL)
	payload.Vs *= 1.0e+3;
      
      // convert Density from g/cm^3 to kg/m^3
      if (payload.Density != cencalvm::storage::Payload::NODATAVAL)
	payload.Density *= 1.0e+3;
      
      // add data to etree
      etree_addr_t addr;
      addr.level = level;
      addr.type = ETREE_LEAF;
      vmgeom.lonLatElevToAddr(&addr, lon, lat, elev);
      
      if (0 != etree_insert(*pDB, addr, &payload)) {
	errHandler.error(etree_strerror(etree_errno(*pDB)));
	break;
      } // if
      numAdded++;
    } else {
      // convert elev and depth from km to m
      if (elev != cencalvm::storage::Payload::NODATAVAL)
	elev *= 1.0e+3;
      if (payload.DepthFreeSurf != cencalvm::storage::Payload::NODATAVAL)
	payload.DepthFreeSurf *= 1.0e+3;

      if (payload.FaultBlock != cencalvm::storage::Payload::NODATABLOCK) {
	std::ostringstream msg;
	msg
	  << std::resetiosflags(std::ios::fixed)
	  << std::setiosflags(std::ios::scientific)
	  << std::setprecision(6)
	  << lon << ", " << lat << ", " << elev << ", No fault block\n";
	errHandler.log(msg.str().c_str());
	numIgnored++;
      } else {
	std::ostringstream msg;
	msg
	  << std::resetiosflags(std::ios::fixed)
	  << std::setiosflags(std::ios::scientific)
	  << std::setprecision(6)
	  << lon << ", " << lat << ", " << elev << ", Ignoring\n";
	errHandler.log(msg.str().c_str());
	numIgnored++;
      } // if/else
      
    } // if/else
  } // for

  fin.close();

  if (cencalvm::storage::ErrorHandler::OK != errHandler.status()) {
    std::ostringstream msg;
    msg << "Caught error while reading grid from '" << filename << "'.\n"
	<< "Successfully added " << numAdded << " points and ignored "
	<< numIgnored << " others.\n"
	<< "Error message: " << errHandler.message();
    errHandler.error(msg.str().c_str());
    return;
  } // if

  if (!quiet)
    std::cout << "Done procesing '" << filename << "'."
	      << "  # points added: " << numAdded
	      << ",  # points ignored: " << numIgnored
	      << std::endl;
} // addGrid

// version
// $Id: GridIngester.cc,v 1.14 2005/10/11 18:25:27 brad Exp $

// End of file
