#include "bindings.hpp"
#include "Real3D.hpp"
#include "Property.hpp"
#include "esutil/Collectives.hpp"
#include <hello/bindings.hpp>
#include <bc/bindings.hpp>
#include <particles/bindings.hpp>
#include <pairs/bindings.hpp>
#include <interaction/bindings.hpp>
#include <integrator/bindings.hpp>
#include <thermostat/bindings.hpp>
#include <esutil/PyLogger.hpp>

void espresso::registerPython() {

  espresso::registerPythonReal3D();
  espresso::registerPythonParticle();
  espresso::registerPythonProperties();
  espresso::esutil::Collectives::registerPython();
  espresso::hello::registerPython();
  espresso::bc::registerPython();
  espresso::pairs::registerPython();
  espresso::interaction::registerPython();
  espresso::particles::registerPython();
  espresso::integrator::registerPython();
  espresso::thermostat::registerPython();

  log4espp::PyLogger::registerPython();
}
