/*
 * ClassTemplate.h
 *
 *  Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner

		This program is free software: you can redistribute it and/or modify
		it under the terms of the GNU General Public License as published by
		the Free Software Foundation, either version 3 of the License, or
		(at your option) any later version.

		This program is distributed in the hope that it will be useful,
		but WITHOUT ANY WARRANTY; without even the implied warranty of
		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
		GNU General Public License for more details.

		You should have received a copy of the GNU General Public License
		along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_CLASSPOTENTIALSMRSS_H_
#define SRC_CLASSPOTENTIALSMRSS_H_

#include "ClassPotentialOrigin.h"

class Class_Potential_SMRSS : public Class_Potential_Origin
{
public:
    Class_Potential_SMRSS ();
  virtual
  ~Class_Potential_SMRSS ();


  // Add here your parameters for the Lagrangian as well as for the counterterm potential
  // Add here your variables in which you will save the Debye correction factors

    double mu_hs, lambda_h, mu_s, lambda_s, lambda_m, Dmu_hs, Dlambda_h, Dmu_s, Dlambda_s, Dlambda_m, DT3, DT5;




    void ReadAndSet(const std::string& linestr, std::vector<double>& par);
    std::string addLegendCT();
    std::string addLegendTemp();
    std::string addLegendTripleCouplings();
    std::string addLegendVEV();

    void set_gen(const std::vector<double>& par);
    void set_CT_Pot_Par(const std::vector<double>& par);
    void write();

    void TripleHiggsCouplings();
    void calc_CT(std::vector<double>& par);


    void MinimizeOrderVEV(const std::vector<double>& vevminimizer, std::vector<double>& vevFunction);

    void SetCurvatureArrays();
    bool CalculateDebyeSimplified();
    bool CalculateDebyeGaugeSimplified();
    double VTreeSimplified(const std::vector<double>& v);
    double VCounterSimplified(const std::vector<double>& v);
    void Debugging(const std::vector<double>& input, std::vector<double>& output);
};

#endif /* SRC_CLASSPOTENTIALSMRSS_H_ */
