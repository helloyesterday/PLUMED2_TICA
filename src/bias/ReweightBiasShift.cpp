/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ReweightBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_BIAS
/*
Calculate weights for ensemble averages that negate the effect the bias has on the region of phase space explored

If a static or pseudo-static bias \f$V(x,t')\f$ is acting on 
the system we can remove the bias and get the unbiased probability distribution using:

\f[
\langle P(s',t) \rangle = \frac{ \sum_{t'}^t \delta( s(x) - s' ) \exp\left( +\frac{V(x,t')}{k_B T} \right) }{ \sum_t'^t \exp\left( +\frac{V(x,t')}{k_B T} \right) } 
\f]

The weights calculated by this action are equal to \f$\exp\left( +\frac{V(x,t')}{k_B T} \right)\f$ these weights can then be used in any action 
that computes ensemble averages.  For example this action can be used in tandem with \ref HISTOGRAM or \ref AVERAGE. 

\par Examples

In the following example there is a fixed restraint on the distance between atoms 1 and 2.  Clearly, this 
restraint will have an effect on the region of phase space that will be sampled when an MD simulation is
run using this variable.  Consequently, when the histogram as a function of the distance, \f$x\f$, is accumulated,
we use reweighting into order to discount the effect of the bias from our final histogram.

\verbatim 
x: DISTANCE ATOMS=1,2
RESTRAINT ARG=x SLOPE=1.0 AT=0.0
as: REWEIGHT_BIAS_SHIFT TEMP=300

HISTOGRAM ...
  ARG=x
  GRID_MIN=0.0
  GRID_MAX=3.0
  GRID_BIN=100
  BANDWIDTH=0.1
  LOGWEIGHTS=bias
  LABEL=hB
... HISTOGRAM

DUMPGRID GRID=hB FILE=histoB STRIDE=1 FMT=%8.4f
\endverbatim

Sometimes, when the value of bias potential are too large or too small, the calculation of reweighting factor at each frame may be beyond the data range.
In this situation, a constant value is needed to shift the value of bias potential. For example, the values of bias potential during the biased simulation
are almost larger than 2000, we can add a shift value -2000 (SHIFT=-2000) to the bias at reweighting:

\verbatim 
x: DISTANCE ATOMS=1,2
RESTRAINT ARG=x SLOPE=1.0 AT=0.0
as: REWEIGHT_BIAS_SHIFT TEMP=300 SHIFT=-2000

HISTOGRAM ...
  ARG=x
  GRID_MIN=0.0
  GRID_MAX=3.0
  GRID_BIN=100
  BANDWIDTH=0.1
  LOGWEIGHTS=bias
  LABEL=hB
... HISTOGRAM

DUMPGRID GRID=hB FILE=histoB STRIDE=1 FMT=%8.4f
\endverbatim


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightBiasShift : public ReweightBase {
private:
  double bshift;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightBiasShift(const ActionOptions&ao);
  double getLogWeight() const ;
};

PLUMED_REGISTER_ACTION(ReweightBiasShift,"REWEIGHT_BIAS_SHIFT")

void ReweightBiasShift::registerKeywords(Keywords& keys ){
  ReweightBase::registerKeywords( keys ); keys.remove("ARG");
  keys.add("compulsory","ARG","*.bias","the biases that must be taken into account when reweighting");
  keys.add("optional","SHIFT","a constant value to shift the bias potential (when they are too large or too small)");
}

ReweightBiasShift::ReweightBiasShift(const ActionOptions&ao):
Action(ao),
ReweightBase(ao),
bshift(0.0)
{
  parse("SHIFT",bshift);
}

double ReweightBiasShift::getLogWeight() const {
   // Retrieve the bias
   double bias=bshift; for(unsigned i=0;i<getNumberOfArguments();++i) bias+=getArgument(i);
   return bias / simtemp;
}

}
}
