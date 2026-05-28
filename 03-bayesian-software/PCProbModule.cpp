
/*
 * PCProbModule.cpp
 * ----------------
 * Registers the pcprob module with JAGS.
 *
 * The global variable pattern (extern "C" jags::Module *pcprob_module)
 * is the mechanism used by this JAGS version: the dynamic loader finds
 * the symbol at load time and JAGS picks it up automatically, without
 * needing an initModule() entry point.
 */

#include "PCProb.h"
#include <module/Module.h>

namespace jags {
namespace custom {

class PCProbModule : public Module {
public:
  PCProbModule() : Module("pcprob") {
    insert(new PCProb());
  }
};

} // namespace custom
} // namespace jags

extern "C" {
  jags::Module *pcprob_module = new jags::custom::PCProbModule();
}
