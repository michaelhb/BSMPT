/*
 * ClassTemplate.cpp
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

#include "ClassPotentialSMRSS.h"
#include "IncludeAllModels.h"
using namespace Eigen;

/**
 * @file
 * Template for adding a new model class
 */

/**
 * Here you have to adjust NNeutralHiggs, NChargedHiggs, nPar (number of Lagrangian parameters AFTER
 *  using the tadpole conditions),
 * nParCT (number of counterterms) as well as nVEV (number of VEVs for minimization)
 */
Class_Potential_SMRSS::Class_Potential_SMRSS ()
{
  Model = C_ModelSMRSS; // global int constant which will be used to tell the program which model is called
  NNeutralHiggs = 5; // number of neutral Higgs bosons at T = 0
  NChargedHiggs=0; // number of charged Higgs bosons  at T = 0 (all d.o.f.)

  nPar = 5; // number of parameters in the tree-Level Lagrangian
  nParCT = 9; // number of parameters in the counterterm potential

  nVEV=2; // number of VEVs to minimize the potential

  NHiggs = NNeutralHiggs+NChargedHiggs;

}

Class_Potential_SMRSS::~Class_Potential_SMRSS ()
{
}

/**
 * returns a string which tells the user the chronological order of the counterterms. Use this to
 * complement the legend of the given input file
 */
std::string Class_Potential_SMRSS::addLegendCT(){
  std::string out = "DmuSq\tDlambda\tDMssq\tDkappa1\tDkappa2\tDkappa\tDlambdas\tDT3\tDT5";
  return out;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and the critical temperature. Use this to
 * complement the legend of the given input file
 */
std::string Class_Potential_SMRSS::addLegendTemp(){
  std::string out = "T_c\tv_c\tv\tvs";
  return out;
}

/**
 * returns a string which tells the user the chronological order of the Triple Higgs couplings. Use this to
 * complement the legend of the given input file
 *
 */
std::string Class_Potential_SMRSS::addLegendTripleCouplings(){
    // MB: leaving this one alone for now since we're not using TripleHiggsNLO

	std::vector<std::string> particles;
	particles.resize(NHiggs);
	//here you have to define the particle names in the vector particles

	particles[0]="H";

	     std::string out = "Tree_";
	     for(int i=0;i<NHiggs;i++)
	       {
	 	for(int j=i;j<NHiggs;j++)
	 	  {
	 	    for(int k=j;k<NHiggs;k++)
	 	      {
	 	    	if(!(i==0 and j==0 and k==0)) out += "\tTree_";
	 		out+=particles.at(i);
	 		out+=particles.at(j);
	 		out+=particles.at(k);
	 		out+="\tCT_";
	 		out+=particles.at(i);
	 		out+=particles.at(j);
	 		out+=particles.at(k);
	 		out+="\tCW_";
	 		out+=particles.at(i);
	 		out+=particles.at(j);
	 		out+=particles.at(k);
	 	      }
	 	  }
	       }

	     return out;
}

/**
 * returns a string which tells the user the chronological order of the VEVs. Use this to
 * complement the legend of the given input file
 */
std::string Class_Potential_SMRSS::addLegendVEV(){
	std::string out;
	//out = "Your VEV order";
	out = "v\tvs";
  return out;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Potential_SMRSS::ReadAndSet(const std::string& linestr, std::vector<double>& par )
{
	std::stringstream ss(linestr);
	double tmp;

	double lcos_theta, lvs, lMS, llambdas, lkappa;

	for(int k=1;k<=5;k++)
	{
	      ss>>tmp;
	      if(k==1) lcos_theta = tmp;
	      else if (k==2) lvs = tmp;
	      else if (k==3) lMS = tmp;
	      else if (k==4) llambdas = tmp;
	      else if (k==5) lkappa = tmp;
	}
	par[0] = lcos_theta;
	par[1] = lvs;
	par[2] = lMS;
	par[3] = llambdas;
	par[4] = lkappa;


	set_gen(par); // This you have to call so that everything will be set
	return ;
}


/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Potential_SMRSS::set_gen(const std::vector<double>& par) {

    // direct input params
    cos_theta = par[0];
    vs = par[1];
    MS = par[2];
    lambdas = par[3];
    kappa = par[4];

    // Input parameters -> gauge basis Lagrangian parameters
    double sin_theta = std::sin(std::acos(cos_theta));
    double cos_2_theta = std::pow(cos_theta, 2) - std::pow(sin_theta,2);
    double sin_2_theta = 2.0*sin_theta*cos_theta;

    double MS2 = std::pow(MS,2);
    double MH2 = std::pow(C_vev0,2);

    lambda = ((2.0)/std::pow(C_vev0,2))*(
      MH2*std::pow(cos_theta,2) +
      MS2*std::pow(sin_theta,2));

    kappa1 = ((2.0*vs)/std::pow(C_vev0,2))*(
      (-1.0)*(MS2 + MH2) +
      vs*(kappa + 4.0*vs*lambdas) +
      (MH2 - MS2)*cos_2_theta);

    kappa2 = ((1.0)/(2*std::pow(C_vev0,2)*vs))*(
      4*vs*(MH2 + MS2 - vs*(kappa + 4*vs*lambdas) + (MS2 - MH2)*cos_2_theta) +
      (MH2 - MS2)*C_vev0*sin_2_theta);

    muSq = (-0.5)*(std::pow(C_vev0,2)*lambda + 2.0*vs*kappa1 + std::pow(vs,2)*kappa2);

    Mssq = (-1.0/(2.0*vs))*(
      2.0*std::pow(vs,2)*kappa + 4*std::pow(vs,3)*lambdas + 
      std::pow(C_vev0, 2)*kappa1 + std::pow(C_vev0,2)*vs*kappa2);

    // This sets the renormalization scale
    scale = C_vev0;

    vevTreeMin.resize(nVEV);
    vevTree.resize(NHiggs);

    // Here you have to set the vector vevTreeMin. The vector vevTree will then be set by the function MinimizeOrderVEV
    vevTreeMin[0] = C_vev0;
    vevTreeMin[1] = vs;

    MinimizeOrderVEV(vevTreeMin,vevTree);
    if(!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the entries of Curvature_Higgs_CT_L1 to
 * Curvature_Higgs_CT_L4.
 */
void Class_Potential_SMRSS::set_CT_Pot_Par(const std::vector<double>& par){

  DmuSq = par[0];
  Dlambda = par[1];
  DMssq = par[2];
  Dkappa1 = par[3];
  Dkappa2 = par[4];
  Dkappa = par[5];
  Dlambdas = par[6];
  DT3 = par[7];
  DT5 = par[8];

  Curvature_Higgs_CT_L1[2] = DT3;
  Curvature_Higgs_CT_L1[4] = DT5;

  Curvature_Higgs_CT_L2[0][0] = DmuSq;
  Curvature_Higgs_CT_L2[1][1] = DmuSq;
  Curvature_Higgs_CT_L2[2][2] = DmuSq;
  Curvature_Higgs_CT_L2[3][3] = DmuSq;
  Curvature_Higgs_CT_L2[4][4] = DMssq;

  Curvature_Higgs_CT_L3[0][0][4] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[0][4][0] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[1][1][4] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[1][4][1] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[2][2][4] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[2][4][2] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[3][3][4] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[3][4][3] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[4][0][0] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[4][1][1] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[4][2][2] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[4][3][3] = 3*Dkappa1;
  Curvature_Higgs_CT_L3[4][4][4] = 2*Dkappa;

  Curvature_Higgs_CT_L4[0][0][0][0] = 3*Dlambda;
  Curvature_Higgs_CT_L4[0][0][1][1] = 6*Dlambda;
  Curvature_Higgs_CT_L4[0][0][2][2] = 6*Dlambda;
  Curvature_Higgs_CT_L4[0][0][3][3] = 6*Dlambda;
  Curvature_Higgs_CT_L4[0][0][4][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[0][1][0][1] = 6*Dlambda;
  Curvature_Higgs_CT_L4[0][1][1][0] = 6*Dlambda;
  Curvature_Higgs_CT_L4[0][2][0][2] = 6*Dlambda;
  Curvature_Higgs_CT_L4[0][2][2][0] = 6*Dlambda;
  Curvature_Higgs_CT_L4[0][3][0][3] = 6*Dlambda;
  Curvature_Higgs_CT_L4[0][3][3][0] = 6*Dlambda;
  Curvature_Higgs_CT_L4[0][4][0][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[0][4][4][0] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[1][0][0][1] = 6*Dlambda;
  Curvature_Higgs_CT_L4[1][0][1][0] = 6*Dlambda;
  Curvature_Higgs_CT_L4[1][1][0][0] = 6*Dlambda;
  Curvature_Higgs_CT_L4[1][1][1][1] = 3*Dlambda;
  Curvature_Higgs_CT_L4[1][1][2][2] = 6*Dlambda;
  Curvature_Higgs_CT_L4[1][1][3][3] = 6*Dlambda;
  Curvature_Higgs_CT_L4[1][1][4][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[1][2][1][2] = 6*Dlambda;
  Curvature_Higgs_CT_L4[1][2][2][1] = 6*Dlambda;
  Curvature_Higgs_CT_L4[1][3][1][3] = 6*Dlambda;
  Curvature_Higgs_CT_L4[1][3][3][1] = 6*Dlambda;
  Curvature_Higgs_CT_L4[1][4][1][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[1][4][4][1] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[2][0][0][2] = 6*Dlambda;
  Curvature_Higgs_CT_L4[2][0][2][0] = 6*Dlambda;
  Curvature_Higgs_CT_L4[2][1][1][2] = 6*Dlambda;
  Curvature_Higgs_CT_L4[2][1][2][1] = 6*Dlambda;
  Curvature_Higgs_CT_L4[2][2][0][0] = 6*Dlambda;
  Curvature_Higgs_CT_L4[2][2][1][1] = 6*Dlambda;
  Curvature_Higgs_CT_L4[2][2][2][2] = 3*Dlambda;
  Curvature_Higgs_CT_L4[2][2][3][3] = 6*Dlambda;
  Curvature_Higgs_CT_L4[2][2][4][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[2][3][2][3] = 6*Dlambda;
  Curvature_Higgs_CT_L4[2][3][3][2] = 6*Dlambda;
  Curvature_Higgs_CT_L4[2][4][2][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[2][4][4][2] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[3][0][0][3] = 6*Dlambda;
  Curvature_Higgs_CT_L4[3][0][3][0] = 6*Dlambda;
  Curvature_Higgs_CT_L4[3][1][1][3] = 6*Dlambda;
  Curvature_Higgs_CT_L4[3][1][3][1] = 6*Dlambda;
  Curvature_Higgs_CT_L4[3][2][2][3] = 6*Dlambda;
  Curvature_Higgs_CT_L4[3][2][3][2] = 6*Dlambda;
  Curvature_Higgs_CT_L4[3][3][0][0] = 6*Dlambda;
  Curvature_Higgs_CT_L4[3][3][1][1] = 6*Dlambda;
  Curvature_Higgs_CT_L4[3][3][2][2] = 6*Dlambda;
  Curvature_Higgs_CT_L4[3][3][3][3] = 3*Dlambda;
  Curvature_Higgs_CT_L4[3][3][4][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[3][4][3][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[3][4][4][3] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][0][0][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][0][4][0] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][1][1][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][1][4][1] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][2][2][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][2][4][2] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][3][3][4] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][3][4][3] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][4][0][0] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][4][1][1] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][4][2][2] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][4][3][3] = 6*Dkappa2;
  Curvature_Higgs_CT_L4[4][4][4][4] = 12*Dlambdas;

}


/**
 * console output of all Parameters
 */
void Class_Potential_SMRSS::write() {

	std::cout << "The input parameters are : " << std::endl;
	std::cout << "cos_theta = " << cos_theta << std::endl
            << "vs = " << vs << std::endl
            << "MS = " << MS << std::endl
            << "lambdas = " << lambdas << std::endl
            << "kappa = " << kappa << std::endl;

	std::cout << "Other lagrangian parameters calculated from these : " << std::endl;
	std::cout << "muSq = " << muSq << std::endl
            << "lambda = " << lambda << std::endl
            << "Mssq = " << Mssq << std::endl
            << "kappa1 = " << kappa1 << std::endl
            << "kappa2 = " << kappa2 << std::endl;

	std::cout << "The counterterm parameters are : " << std::endl;
	std::cout << "DmuSq = " << DmuSq << std::endl
            << "Dlambda = " << Dlambda << std::endl
            << "DMssq = " << DMssq << std::endl
            << "Dkappa1 = " << Dkappa1 << std::endl
            << "Dkappa2 = " << Dkappa2 << std::endl
            << "Dkappa = " << Dkappa << std::endl
            << "Dlambdas = " << Dlambdas << std::endl
            << "DT3 = " << DT3 << std::endl
            << "DT5 = " << DT5 << std::endl;

	std::cout << "The scale is given by mu = " << scale << " GeV " << std::endl;

}


/**
 * Calculates the counterterms. Here you need to work out the scheme and implement the formulas.
 */
void Class_Potential_SMRSS::calc_CT(std::vector<double>& par){
	bool Debug=false;
	if(Debug) std::cout << "Start" << __func__ << std::endl;

	if(!SetCurvatureDone)SetCurvatureArrays();
	if(!CalcCouplingsdone)CalculatePhysicalCouplings();
	if(Debug) {
	std::cout << "Couplings done " << std::endl;
	}
	std::vector<double> WeinbergNabla,WeinbergHesse;
	WeinbergFirstDerivative(WeinbergNabla);
	WeinbergSecondDerivative(WeinbergHesse);

	if(Debug) std::cout << "Finished Derivatives " << std::endl;

	VectorXd Deriv1(NHiggs);
	MatrixXd Deriv2(NHiggs,NHiggs),HiggsRot(NHiggs,NHiggs);
	for(int i=0;i<NHiggs;i++)
	{
		Deriv1[i] = WeinbergNabla[i];
		for(int j=0;j<NHiggs;j++) Deriv2(i,j) = WeinbergHesse.at(j*NHiggs+i);
	}

  DmuSq = (-1.0/(2*C_vev0))*(3.0*Deriv1(2) - vs*Deriv2(2,4) - C_vev0*Deriv2(2,2));
  Dlambda = (-1.0/std::pow(C_vev0,3))*(-Deriv1(2) + C_vev0*Deriv2(2,2));
  Dkappa = (-1.0/(2*std::pow(vs,2)))*(6.0*Deriv1(4) - 2.0*vs*Deriv2(4,4) - 2.0*C_vev0*Deriv2(2,4));
  Dlambdas = (-1.0/(4.0*std::pow(vs,3)))*(-4.0*Deriv1(4) + 2.0*vs*Deriv2(4,4) + C_vev0*Deriv2(2,4));
  Dkappa2 = (-1.0/(C_vev0*vs))*Deriv2(2,4); 
  DT3 = 0;
  DT5 = 0;
  DMssq = 0;
  Dkappa1 = 0;

  par[0] = DmuSq;
  par[1] = Dlambda;
  par[2] = DMssq;
  par[3] = Dkappa1;
  par[4] = Dkappa2;
  par[5] = Dkappa;
  par[6] = Dlambdas;
  par[7] = DT3;
  par[8] = DT5;

}


void Class_Potential_SMRSS::TripleHiggsCouplings()
{
	bool Debug=false;
	if(Debug) std::cout << "Debug turned on in " << __func__ << std::endl;

	if(!SetCurvatureDone)SetCurvatureArrays();
	if(!CalcCouplingsdone)CalculatePhysicalCouplings();


	std::vector<double> HiggsOrder(NHiggs);
	// Here you have to set the vector HiggsOrder. By telling e.g. HiggsOrder[0] = 5 you always want your 6th lightest
	// particle to be the first particle in the vector (which has the index 5 because they are sorted by mass)

	// example for keeping the mass order
	for(int i=0;i<NHiggs;i++) {
		HiggsOrder[i]=i;
	}

	if(Debug) std::cout << "Calculate Derivative" << std::endl;
	std::vector<double> TripleDeriv;
	WeinbergThirdDerivative(TripleDeriv);
	if(Debug) std::cout << "Finished calculating derivatives " << std::endl;
	double GaugeBasis[NHiggs][NHiggs][NHiggs];
	for(int i=0;i<NHiggs;i++)
	  {
		for(int j=0;j<NHiggs;j++)
		{
		  for(int k=0;k<NHiggs;k++)
			{
			  GaugeBasis[i][j][k] = TripleDeriv.at(i+j*NHiggs+k*NHiggs*NHiggs);
			}
		}
	  }

	MatrixXd HiggsRot(NHiggs,NHiggs);
	for(int i=0;i<NHiggs;i++)
	{
		for(int j=0;j<NHiggs;j++)
		{
			HiggsRot(i,j) = HiggsRotationMatrix[i][j];
		}
	}


	MatrixXd HiggsRotSort(NHiggs,NHiggs);




	for(int i=0;i<NHiggs;i++)
	{
		HiggsRotSort.row(i) = HiggsRot.row(HiggsOrder[i]);
	}

	TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
	TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
	TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
	for(int i=0;i<NHiggs;i++) {
		TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
		TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
		for(int j=0;j<NHiggs;j++) {
			TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
			TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
		}
	}

	if(Debug) std::cout << "Setup done " << std::endl;


	for(int i=0;i<NHiggs;i++)
	  {
		for(int j=0;j<NHiggs;j++)
		{
			for(int k=0;k<NHiggs;k++)
			{
			  TripleHiggsCorrectionsCWPhysical[i][j][k] = 0;
			  TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
			  TripleHiggsCorrectionsCTPhysical[i][j][k] = 0;
			  for(int l=0;l<NHiggs;l++)
			  {
				  for(int m=0;m<NHiggs;m++)
				  {
					  for(int n=0;n<NHiggs;n++)
					  {
						  double RotFac = HiggsRotSort(i,l)*HiggsRotSort(j,m)*HiggsRotSort(k,n);
						  TripleHiggsCorrectionsCWPhysical[i][j][k] += RotFac*GaugeBasis[l][m][n];
						  TripleHiggsCorrectionsTreePhysical[i][j][k] += RotFac*LambdaHiggs_3[l][m][n];
						  TripleHiggsCorrectionsCTPhysical[i][j][k] += RotFac*LambdaHiggs_3_CT[l][m][n];

					  }
				  }
			  }
			}
		}
	  }



}


void Class_Potential_SMRSS::SetCurvatureArrays(){
  /*
   *  Here you have to set the vectors
   *  Curvature_Higgs_L1,Curvature_Higgs_L2,Curvature_Higgs_L3,Curvature_Higgs_L4
   *  Curvature_Gauge_G2H2
   *  Curvature_Quark_F2H1, Curvature_Lepton_F2H1
   *  as described in the potential in the paper.
   */

	initVectors();

	for(int i=0;i<NHiggs;i++) HiggsVev[i] = vevTree[i];

    // Set up Yukawa couplings for weak eigenbasis
    Matrix3cd ckm, Mu, Md, Me, Yu, Yd, Ye;
    ckm << C_Vud, C_Vus, C_Vub,
           C_Vcd, C_Vcs, C_Vcb,
           C_Vtd, C_Vts, C_Vtb;

    Mu << C_MassUp, 0, 0,
          0, C_MassCharm, 0,
          0, 0, C_MassTop;

    Md << C_MassDown, 0, 0,
          0, C_MassStrange, 0,
          0, 0, C_MassBottom;

    Me << C_MassElectron, 0, 0,
          0, C_MassMu, 0,
          0, 0, C_MassTau;

    Yu = (std::sqrt(2) / C_vev0)*(ckm.adjoint())*Mu*ckm;
    Yd = (std::sqrt(2) / C_vev0)*ckm*Md*(ckm.adjoint());
    Ye = (std::sqrt(2) / C_vev0)*Me;

    std::complex<double> II(0,1);

    // NB the curvature arrays are initialised to zero, so
    // only writing the non zero entries here.
    Curvature_Quark_F2H1[0][1][2] = 2*std::sqrt(2)*Yu(1,1);
    Curvature_Quark_F2H1[0][1][3] = (-2.0*II)*std::sqrt(2)*Yu(1,1);
    Curvature_Quark_F2H1[0][3][0] = 2*std::sqrt(2)*Yd(1,1);
    Curvature_Quark_F2H1[0][3][1] = (2.0*II)*std::sqrt(2)*Yd(1,1);
    Curvature_Quark_F2H1[0][5][2] = 2*std::sqrt(2)*Yu(1,2);
    Curvature_Quark_F2H1[0][5][3] = (-2.0*II)*std::sqrt(2)*Yu(1,2);
    Curvature_Quark_F2H1[0][7][0] = 2*std::sqrt(2)*Yd(1,2);
    Curvature_Quark_F2H1[0][7][1] = (2.0*II)*std::sqrt(2)*Yd(1,2);
    Curvature_Quark_F2H1[0][9][2] = 2*std::sqrt(2)*Yu(1,3);
    Curvature_Quark_F2H1[0][9][3] = (-2.0*II)*std::sqrt(2)*Yu(1,3);
    Curvature_Quark_F2H1[0][11][0] = 2*std::sqrt(2)*Yd(1,3);
    Curvature_Quark_F2H1[0][11][1] = (2.0*II)*std::sqrt(2)*Yd(1,3);
    Curvature_Quark_F2H1[1][0][2] = 2*std::sqrt(2)*Yu(1,1);
    Curvature_Quark_F2H1[1][0][3] = (-2.0*II)*std::sqrt(2)*Yu(1,1);
    Curvature_Quark_F2H1[1][2][0] = -2*std::sqrt(2)*Yu(1,1);
    Curvature_Quark_F2H1[1][2][1] = (2.0*II)*std::sqrt(2)*Yu(1,1);
    Curvature_Quark_F2H1[1][4][2] = 2*std::sqrt(2)*Yu(2,1);
    Curvature_Quark_F2H1[1][4][3] = (-2.0*II)*std::sqrt(2)*Yu(2,1);
    Curvature_Quark_F2H1[1][6][0] = -2*std::sqrt(2)*Yu(2,1);
    Curvature_Quark_F2H1[1][6][1] = (2.0*II)*std::sqrt(2)*Yu(2,1);
    Curvature_Quark_F2H1[1][8][2] = 2*std::sqrt(2)*Yu(3,1);
    Curvature_Quark_F2H1[1][8][3] = (-2.0*II)*std::sqrt(2)*Yu(3,1);
    Curvature_Quark_F2H1[1][10][0] = -2*std::sqrt(2)*Yu(3,1);
    Curvature_Quark_F2H1[1][10][1] = (2.0*II)*std::sqrt(2)*Yu(3,1);
    Curvature_Quark_F2H1[2][1][0] = -2*std::sqrt(2)*Yu(1,1);
    Curvature_Quark_F2H1[2][1][1] = (2.0*II)*std::sqrt(2)*Yu(1,1);
    Curvature_Quark_F2H1[2][3][2] = 2*std::sqrt(2)*Yd(1,1);
    Curvature_Quark_F2H1[2][3][3] = (2.0*II)*std::sqrt(2)*Yd(1,1);
    Curvature_Quark_F2H1[2][5][0] = -2*std::sqrt(2)*Yu(1,2);
    Curvature_Quark_F2H1[2][5][1] = (2.0*II)*std::sqrt(2)*Yu(1,2);
    Curvature_Quark_F2H1[2][7][2] = 2*std::sqrt(2)*Yd(1,2);
    Curvature_Quark_F2H1[2][7][3] = (2.0*II)*std::sqrt(2)*Yd(1,2);
    Curvature_Quark_F2H1[2][9][0] = -2*std::sqrt(2)*Yu(1,3);
    Curvature_Quark_F2H1[2][9][1] = (2.0*II)*std::sqrt(2)*Yu(1,3);
    Curvature_Quark_F2H1[2][11][2] = 2*std::sqrt(2)*Yd(1,3);
    Curvature_Quark_F2H1[2][11][3] = (2.0*II)*std::sqrt(2)*Yd(1,3);
    Curvature_Quark_F2H1[3][0][0] = 2*std::sqrt(2)*Yd(1,1);
    Curvature_Quark_F2H1[3][0][1] = (2.0*II)*std::sqrt(2)*Yd(1,1);
    Curvature_Quark_F2H1[3][2][2] = 2*std::sqrt(2)*Yd(1,1);
    Curvature_Quark_F2H1[3][2][3] = (2.0*II)*std::sqrt(2)*Yd(1,1);
    Curvature_Quark_F2H1[3][4][0] = 2*std::sqrt(2)*Yd(2,1);
    Curvature_Quark_F2H1[3][4][1] = (2.0*II)*std::sqrt(2)*Yd(2,1);
    Curvature_Quark_F2H1[3][6][2] = 2*std::sqrt(2)*Yd(2,1);
    Curvature_Quark_F2H1[3][6][3] = (2.0*II)*std::sqrt(2)*Yd(2,1);
    Curvature_Quark_F2H1[3][8][0] = 2*std::sqrt(2)*Yd(3,1);
    Curvature_Quark_F2H1[3][8][1] = (2.0*II)*std::sqrt(2)*Yd(3,1);
    Curvature_Quark_F2H1[3][10][2] = 2*std::sqrt(2)*Yd(3,1);
    Curvature_Quark_F2H1[3][10][3] = (2.0*II)*std::sqrt(2)*Yd(3,1);
    Curvature_Quark_F2H1[4][1][2] = 2*std::sqrt(2)*Yu(2,1);
    Curvature_Quark_F2H1[4][1][3] = (-2.0*II)*std::sqrt(2)*Yu(2,1);
    Curvature_Quark_F2H1[4][3][0] = 2*std::sqrt(2)*Yd(2,1);
    Curvature_Quark_F2H1[4][3][1] = (2.0*II)*std::sqrt(2)*Yd(2,1);
    Curvature_Quark_F2H1[4][5][2] = 2*std::sqrt(2)*Yu(2,2);
    Curvature_Quark_F2H1[4][5][3] = (-2.0*II)*std::sqrt(2)*Yu(2,2);
    Curvature_Quark_F2H1[4][7][0] = 2*std::sqrt(2)*Yd(2,2);
    Curvature_Quark_F2H1[4][7][1] = (2.0*II)*std::sqrt(2)*Yd(2,2);
    Curvature_Quark_F2H1[4][9][2] = 2*std::sqrt(2)*Yu(2,3);
    Curvature_Quark_F2H1[4][9][3] = (-2.0*II)*std::sqrt(2)*Yu(2,3);
    Curvature_Quark_F2H1[4][11][0] = 2*std::sqrt(2)*Yd(2,3);
    Curvature_Quark_F2H1[4][11][1] = (2.0*II)*std::sqrt(2)*Yd(2,3);
    Curvature_Quark_F2H1[5][0][2] = 2*std::sqrt(2)*Yu(1,2);
    Curvature_Quark_F2H1[5][0][3] = (-2.0*II)*std::sqrt(2)*Yu(1,2);
    Curvature_Quark_F2H1[5][2][0] = -2*std::sqrt(2)*Yu(1,2);
    Curvature_Quark_F2H1[5][2][1] = (2.0*II)*std::sqrt(2)*Yu(1,2);
    Curvature_Quark_F2H1[5][4][2] = 2*std::sqrt(2)*Yu(2,2);
    Curvature_Quark_F2H1[5][4][3] = (-2.0*II)*std::sqrt(2)*Yu(2,2);
    Curvature_Quark_F2H1[5][6][0] = -2*std::sqrt(2)*Yu(2,2);
    Curvature_Quark_F2H1[5][6][1] = (2.0*II)*std::sqrt(2)*Yu(2,2);
    Curvature_Quark_F2H1[5][8][2] = 2*std::sqrt(2)*Yu(3,2);
    Curvature_Quark_F2H1[5][8][3] = (-2.0*II)*std::sqrt(2)*Yu(3,2);
    Curvature_Quark_F2H1[5][10][0] = -2*std::sqrt(2)*Yu(3,2);
    Curvature_Quark_F2H1[5][10][1] = (2.0*II)*std::sqrt(2)*Yu(3,2);
    Curvature_Quark_F2H1[6][1][0] = -2*std::sqrt(2)*Yu(2,1);
    Curvature_Quark_F2H1[6][1][1] = (2.0*II)*std::sqrt(2)*Yu(2,1);
    Curvature_Quark_F2H1[6][3][2] = 2*std::sqrt(2)*Yd(2,1);
    Curvature_Quark_F2H1[6][3][3] = (2.0*II)*std::sqrt(2)*Yd(2,1);
    Curvature_Quark_F2H1[6][5][0] = -2*std::sqrt(2)*Yu(2,2);
    Curvature_Quark_F2H1[6][5][1] = (2.0*II)*std::sqrt(2)*Yu(2,2);
    Curvature_Quark_F2H1[6][7][2] = 2*std::sqrt(2)*Yd(2,2);
    Curvature_Quark_F2H1[6][7][3] = (2.0*II)*std::sqrt(2)*Yd(2,2);
    Curvature_Quark_F2H1[6][9][0] = -2*std::sqrt(2)*Yu(2,3);
    Curvature_Quark_F2H1[6][9][1] = (2.0*II)*std::sqrt(2)*Yu(2,3);
    Curvature_Quark_F2H1[6][11][2] = 2*std::sqrt(2)*Yd(2,3);
    Curvature_Quark_F2H1[6][11][3] = (2.0*II)*std::sqrt(2)*Yd(2,3);
    Curvature_Quark_F2H1[7][0][0] = 2*std::sqrt(2)*Yd(1,2);
    Curvature_Quark_F2H1[7][0][1] = (2.0*II)*std::sqrt(2)*Yd(1,2);
    Curvature_Quark_F2H1[7][2][2] = 2*std::sqrt(2)*Yd(1,2);
    Curvature_Quark_F2H1[7][2][3] = (2.0*II)*std::sqrt(2)*Yd(1,2);
    Curvature_Quark_F2H1[7][4][0] = 2*std::sqrt(2)*Yd(2,2);
    Curvature_Quark_F2H1[7][4][1] = (2.0*II)*std::sqrt(2)*Yd(2,2);
    Curvature_Quark_F2H1[7][6][2] = 2*std::sqrt(2)*Yd(2,2);
    Curvature_Quark_F2H1[7][6][3] = (2.0*II)*std::sqrt(2)*Yd(2,2);
    Curvature_Quark_F2H1[7][8][0] = 2*std::sqrt(2)*Yd(3,2);
    Curvature_Quark_F2H1[7][8][1] = (2.0*II)*std::sqrt(2)*Yd(3,2);
    Curvature_Quark_F2H1[7][10][2] = 2*std::sqrt(2)*Yd(3,2);
    Curvature_Quark_F2H1[7][10][3] = (2.0*II)*std::sqrt(2)*Yd(3,2);
    Curvature_Quark_F2H1[8][1][2] = 2*std::sqrt(2)*Yu(3,1);
    Curvature_Quark_F2H1[8][1][3] = (-2.0*II)*std::sqrt(2)*Yu(3,1);
    Curvature_Quark_F2H1[8][3][0] = 2*std::sqrt(2)*Yd(3,1);
    Curvature_Quark_F2H1[8][3][1] = (2.0*II)*std::sqrt(2)*Yd(3,1);
    Curvature_Quark_F2H1[8][5][2] = 2*std::sqrt(2)*Yu(3,2);
    Curvature_Quark_F2H1[8][5][3] = (-2.0*II)*std::sqrt(2)*Yu(3,2);
    Curvature_Quark_F2H1[8][7][0] = 2*std::sqrt(2)*Yd(3,2);
    Curvature_Quark_F2H1[8][7][1] = (2.0*II)*std::sqrt(2)*Yd(3,2);
    Curvature_Quark_F2H1[8][9][2] = 2*std::sqrt(2)*Yu(3,3);
    Curvature_Quark_F2H1[8][9][3] = (-2.0*II)*std::sqrt(2)*Yu(3,3);
    Curvature_Quark_F2H1[8][11][0] = 2*std::sqrt(2)*Yd(3,3);
    Curvature_Quark_F2H1[8][11][1] = (2.0*II)*std::sqrt(2)*Yd(3,3);
    Curvature_Quark_F2H1[9][0][2] = 2*std::sqrt(2)*Yu(1,3);
    Curvature_Quark_F2H1[9][0][3] = (-2.0*II)*std::sqrt(2)*Yu(1,3);
    Curvature_Quark_F2H1[9][2][0] = -2*std::sqrt(2)*Yu(1,3);
    Curvature_Quark_F2H1[9][2][1] = (2.0*II)*std::sqrt(2)*Yu(1,3);
    Curvature_Quark_F2H1[9][4][2] = 2*std::sqrt(2)*Yu(2,3);
    Curvature_Quark_F2H1[9][4][3] = (-2.0*II)*std::sqrt(2)*Yu(2,3);
    Curvature_Quark_F2H1[9][6][0] = -2*std::sqrt(2)*Yu(2,3);
    Curvature_Quark_F2H1[9][6][1] = (2.0*II)*std::sqrt(2)*Yu(2,3);
    Curvature_Quark_F2H1[9][8][2] = 2*std::sqrt(2)*Yu(3,3);
    Curvature_Quark_F2H1[9][8][3] = (-2.0*II)*std::sqrt(2)*Yu(3,3);
    Curvature_Quark_F2H1[9][10][0] = -2*std::sqrt(2)*Yu(3,3);
    Curvature_Quark_F2H1[9][10][1] = (2.0*II)*std::sqrt(2)*Yu(3,3);
    Curvature_Quark_F2H1[10][1][0] = -2*std::sqrt(2)*Yu(3,1);
    Curvature_Quark_F2H1[10][1][1] = (2.0*II)*std::sqrt(2)*Yu(3,1);
    Curvature_Quark_F2H1[10][3][2] = 2*std::sqrt(2)*Yd(3,1);
    Curvature_Quark_F2H1[10][3][3] = (2.0*II)*std::sqrt(2)*Yd(3,1);
    Curvature_Quark_F2H1[10][5][0] = -2*std::sqrt(2)*Yu(3,2);
    Curvature_Quark_F2H1[10][5][1] = (2.0*II)*std::sqrt(2)*Yu(3,2);
    Curvature_Quark_F2H1[10][7][2] = 2*std::sqrt(2)*Yd(3,2);
    Curvature_Quark_F2H1[10][7][3] = (2.0*II)*std::sqrt(2)*Yd(3,2);
    Curvature_Quark_F2H1[10][9][0] = -2*std::sqrt(2)*Yu(3,3);
    Curvature_Quark_F2H1[10][9][1] = (2.0*II)*std::sqrt(2)*Yu(3,3);
    Curvature_Quark_F2H1[10][11][2] = 2*std::sqrt(2)*Yd(3,3);
    Curvature_Quark_F2H1[10][11][3] = (2.0*II)*std::sqrt(2)*Yd(3,3);
    Curvature_Quark_F2H1[11][0][0] = 2*std::sqrt(2)*Yd(1,3);
    Curvature_Quark_F2H1[11][0][1] = (2.0*II)*std::sqrt(2)*Yd(1,3);
    Curvature_Quark_F2H1[11][2][2] = 2*std::sqrt(2)*Yd(1,3);
    Curvature_Quark_F2H1[11][2][3] = (2.0*II)*std::sqrt(2)*Yd(1,3);
    Curvature_Quark_F2H1[11][4][0] = 2*std::sqrt(2)*Yd(2,3);
    Curvature_Quark_F2H1[11][4][1] = (2.0*II)*std::sqrt(2)*Yd(2,3);
    Curvature_Quark_F2H1[11][6][2] = 2*std::sqrt(2)*Yd(2,3);
    Curvature_Quark_F2H1[11][6][3] = (2.0*II)*std::sqrt(2)*Yd(2,3);
    Curvature_Quark_F2H1[11][8][0] = 2*std::sqrt(2)*Yd(3,3);
    Curvature_Quark_F2H1[11][8][1] = (2.0*II)*std::sqrt(2)*Yd(3,3);
    Curvature_Quark_F2H1[11][10][2] = 2*std::sqrt(2)*Yd(3,3);
    Curvature_Quark_F2H1[11][10][3] = (2.0*II)*std::sqrt(2)*Yd(3,3);

    Curvature_Lepton_F2H1[0][1][2] = 2*std::sqrt(2)*Ye(1,1);
    Curvature_Lepton_F2H1[0][1][3] = (2.0*II)*std::sqrt(2)*Ye(1,1);
    Curvature_Lepton_F2H1[0][4][2] = 2*std::sqrt(2)*Ye(1,2);
    Curvature_Lepton_F2H1[0][4][3] = (2.0*II)*std::sqrt(2)*Ye(1,2);
    Curvature_Lepton_F2H1[0][7][2] = 2*std::sqrt(2)*Ye(1,3);
    Curvature_Lepton_F2H1[0][7][3] = (2.0*II)*std::sqrt(2)*Ye(1,3);
    Curvature_Lepton_F2H1[1][0][2] = 2*std::sqrt(2)*Ye(1,1);
    Curvature_Lepton_F2H1[1][0][3] = (2.0*II)*std::sqrt(2)*Ye(1,1);
    Curvature_Lepton_F2H1[1][2][0] = 2*std::sqrt(2)*Ye(1,1);
    Curvature_Lepton_F2H1[1][2][1] = (2.0*II)*std::sqrt(2)*Ye(1,1);
    Curvature_Lepton_F2H1[1][3][2] = 2*std::sqrt(2)*Ye(2,1);
    Curvature_Lepton_F2H1[1][3][3] = (2.0*II)*std::sqrt(2)*Ye(2,1);
    Curvature_Lepton_F2H1[1][5][0] = 2*std::sqrt(2)*Ye(2,1);
    Curvature_Lepton_F2H1[1][5][1] = (2.0*II)*std::sqrt(2)*Ye(2,1);
    Curvature_Lepton_F2H1[1][6][2] = 2*std::sqrt(2)*Ye(3,1);
    Curvature_Lepton_F2H1[1][6][3] = (2.0*II)*std::sqrt(2)*Ye(3,1);
    Curvature_Lepton_F2H1[1][8][0] = 2*std::sqrt(2)*Ye(3,1);
    Curvature_Lepton_F2H1[1][8][1] = (2.0*II)*std::sqrt(2)*Ye(3,1);
    Curvature_Lepton_F2H1[2][1][0] = 2*std::sqrt(2)*Ye(1,1);
    Curvature_Lepton_F2H1[2][1][1] = (2.0*II)*std::sqrt(2)*Ye(1,1);
    Curvature_Lepton_F2H1[2][4][0] = 2*std::sqrt(2)*Ye(1,2);
    Curvature_Lepton_F2H1[2][4][1] = (2.0*II)*std::sqrt(2)*Ye(1,2);
    Curvature_Lepton_F2H1[2][7][0] = 2*std::sqrt(2)*Ye(1,3);
    Curvature_Lepton_F2H1[2][7][1] = (2.0*II)*std::sqrt(2)*Ye(1,3);
    Curvature_Lepton_F2H1[3][1][2] = 2*std::sqrt(2)*Ye(2,1);
    Curvature_Lepton_F2H1[3][1][3] = (2.0*II)*std::sqrt(2)*Ye(2,1);
    Curvature_Lepton_F2H1[3][4][2] = 2*std::sqrt(2)*Ye(2,2);
    Curvature_Lepton_F2H1[3][4][3] = (2.0*II)*std::sqrt(2)*Ye(2,2);
    Curvature_Lepton_F2H1[3][7][2] = 2*std::sqrt(2)*Ye(2,3);
    Curvature_Lepton_F2H1[3][7][3] = (2.0*II)*std::sqrt(2)*Ye(2,3);
    Curvature_Lepton_F2H1[4][0][2] = 2*std::sqrt(2)*Ye(1,2);
    Curvature_Lepton_F2H1[4][0][3] = (2.0*II)*std::sqrt(2)*Ye(1,2);
    Curvature_Lepton_F2H1[4][2][0] = 2*std::sqrt(2)*Ye(1,2);
    Curvature_Lepton_F2H1[4][2][1] = (2.0*II)*std::sqrt(2)*Ye(1,2);
    Curvature_Lepton_F2H1[4][3][2] = 2*std::sqrt(2)*Ye(2,2);
    Curvature_Lepton_F2H1[4][3][3] = (2.0*II)*std::sqrt(2)*Ye(2,2);
    Curvature_Lepton_F2H1[4][5][0] = 2*std::sqrt(2)*Ye(2,2);
    Curvature_Lepton_F2H1[4][5][1] = (2.0*II)*std::sqrt(2)*Ye(2,2);
    Curvature_Lepton_F2H1[4][6][2] = 2*std::sqrt(2)*Ye(3,2);
    Curvature_Lepton_F2H1[4][6][3] = (2.0*II)*std::sqrt(2)*Ye(3,2);
    Curvature_Lepton_F2H1[4][8][0] = 2*std::sqrt(2)*Ye(3,2);
    Curvature_Lepton_F2H1[4][8][1] = (2.0*II)*std::sqrt(2)*Ye(3,2);
    Curvature_Lepton_F2H1[5][1][0] = 2*std::sqrt(2)*Ye(2,1);
    Curvature_Lepton_F2H1[5][1][1] = (2.0*II)*std::sqrt(2)*Ye(2,1);
    Curvature_Lepton_F2H1[5][4][0] = 2*std::sqrt(2)*Ye(2,2);
    Curvature_Lepton_F2H1[5][4][1] = (2.0*II)*std::sqrt(2)*Ye(2,2);
    Curvature_Lepton_F2H1[5][7][0] = 2*std::sqrt(2)*Ye(2,3);
    Curvature_Lepton_F2H1[5][7][1] = (2.0*II)*std::sqrt(2)*Ye(2,3);
    Curvature_Lepton_F2H1[6][1][2] = 2*std::sqrt(2)*Ye(3,1);
    Curvature_Lepton_F2H1[6][1][3] = (2.0*II)*std::sqrt(2)*Ye(3,1);
    Curvature_Lepton_F2H1[6][4][2] = 2*std::sqrt(2)*Ye(3,2);
    Curvature_Lepton_F2H1[6][4][3] = (2.0*II)*std::sqrt(2)*Ye(3,2);
    Curvature_Lepton_F2H1[6][7][2] = 2*std::sqrt(2)*Ye(3,3);
    Curvature_Lepton_F2H1[6][7][3] = (2.0*II)*std::sqrt(2)*Ye(3,3);
    Curvature_Lepton_F2H1[7][0][2] = 2*std::sqrt(2)*Ye(1,3);
    Curvature_Lepton_F2H1[7][0][3] = (2.0*II)*std::sqrt(2)*Ye(1,3);
    Curvature_Lepton_F2H1[7][2][0] = 2*std::sqrt(2)*Ye(1,3);
    Curvature_Lepton_F2H1[7][2][1] = (2.0*II)*std::sqrt(2)*Ye(1,3);
    Curvature_Lepton_F2H1[7][3][2] = 2*std::sqrt(2)*Ye(2,3);
    Curvature_Lepton_F2H1[7][3][3] = (2.0*II)*std::sqrt(2)*Ye(2,3);
    Curvature_Lepton_F2H1[7][5][0] = 2*std::sqrt(2)*Ye(2,3);
    Curvature_Lepton_F2H1[7][5][1] = (2.0*II)*std::sqrt(2)*Ye(2,3);
    Curvature_Lepton_F2H1[7][6][2] = 2*std::sqrt(2)*Ye(3,3);
    Curvature_Lepton_F2H1[7][6][3] = (2.0*II)*std::sqrt(2)*Ye(3,3);
    Curvature_Lepton_F2H1[7][8][0] = 2*std::sqrt(2)*Ye(3,3);
    Curvature_Lepton_F2H1[7][8][1] = (2.0*II)*std::sqrt(2)*Ye(3,3);
    Curvature_Lepton_F2H1[8][1][0] = 2*std::sqrt(2)*Ye(3,1);
    Curvature_Lepton_F2H1[8][1][1] = (2.0*II)*std::sqrt(2)*Ye(3,1);
    Curvature_Lepton_F2H1[8][4][0] = 2*std::sqrt(2)*Ye(3,2);
    Curvature_Lepton_F2H1[8][4][1] = (2.0*II)*std::sqrt(2)*Ye(3,2);
    Curvature_Lepton_F2H1[8][7][0] = 2*std::sqrt(2)*Ye(3,3);
    Curvature_Lepton_F2H1[8][7][1] = (2.0*II)*std::sqrt(2)*Ye(3,3);

    Curvature_Gauge_G2H2[0][0][0][0] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[0][0][1][1] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[0][0][2][2] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[0][0][3][3] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[0][3][0][2] = 2*(C_g)*(C_gs);
    Curvature_Gauge_G2H2[0][3][1][3] = 2*(C_g)*(C_gs);
    Curvature_Gauge_G2H2[0][3][2][0] = 2*(C_g)*(C_gs);
    Curvature_Gauge_G2H2[0][3][3][1] = 2*(C_g)*(C_gs);
    Curvature_Gauge_G2H2[1][1][0][0] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[1][1][1][1] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[1][1][2][2] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[1][1][3][3] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[1][2][0][3] = -2*std::pow(C_g,2);
    Curvature_Gauge_G2H2[1][2][1][2] = 2*std::pow(C_g,2);
    Curvature_Gauge_G2H2[1][2][2][1] = 2*std::pow(C_g,2);
    Curvature_Gauge_G2H2[1][2][3][0] = -2*std::pow(C_g,2);
    Curvature_Gauge_G2H2[2][1][0][3] = -2*std::pow(C_g,2);
    Curvature_Gauge_G2H2[2][1][1][2] = 2*std::pow(C_g,2);
    Curvature_Gauge_G2H2[2][1][2][1] = 2*std::pow(C_g,2);
    Curvature_Gauge_G2H2[2][1][3][0] = -2*std::pow(C_g,2);
    Curvature_Gauge_G2H2[2][2][0][0] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[2][2][1][1] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[2][2][2][2] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[2][2][3][3] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[2][3][0][0] = (C_g)*(C_gs);
    Curvature_Gauge_G2H2[2][3][1][1] = (C_g)*(C_gs);
    Curvature_Gauge_G2H2[2][3][2][2] = -((C_g)*(C_gs));
    Curvature_Gauge_G2H2[2][3][3][3] = -((C_g)*(C_gs));
    Curvature_Gauge_G2H2[3][0][0][2] = 2*(C_g)*(C_gs);
    Curvature_Gauge_G2H2[3][0][1][3] = 2*(C_g)*(C_gs);
    Curvature_Gauge_G2H2[3][0][2][0] = 2*(C_g)*(C_gs);
    Curvature_Gauge_G2H2[3][0][3][1] = 2*(C_g)*(C_gs);
    Curvature_Gauge_G2H2[3][2][0][0] = (C_g)*(C_gs);
    Curvature_Gauge_G2H2[3][2][1][1] = (C_g)*(C_gs);
    Curvature_Gauge_G2H2[3][2][2][2] = -((C_g)*(C_gs));
    Curvature_Gauge_G2H2[3][2][3][3] = -((C_g)*(C_gs));
    Curvature_Gauge_G2H2[3][3][0][0] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[3][3][1][1] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[3][3][2][2] = std::pow(C_g,2)/2;
    Curvature_Gauge_G2H2[3][3][3][3] = std::pow(C_g,2)/2;

    Curvature_Higgs_L2[0][0] = muSq;
    Curvature_Higgs_L2[1][1] = muSq;
    Curvature_Higgs_L2[2][2] = muSq;
    Curvature_Higgs_L2[3][3] = muSq;
    Curvature_Higgs_L2[4][4] = Mssq;

    Curvature_Higgs_L3[0][0][4] = 3*kappa1;
    Curvature_Higgs_L3[0][4][0] = 3*kappa1;
    Curvature_Higgs_L3[1][1][4] = 3*kappa1;
    Curvature_Higgs_L3[1][4][1] = 3*kappa1;
    Curvature_Higgs_L3[2][2][4] = 3*kappa1;
    Curvature_Higgs_L3[2][4][2] = 3*kappa1;
    Curvature_Higgs_L3[3][3][4] = 3*kappa1;
    Curvature_Higgs_L3[3][4][3] = 3*kappa1;
    Curvature_Higgs_L3[4][0][0] = 3*kappa1;
    Curvature_Higgs_L3[4][1][1] = 3*kappa1;
    Curvature_Higgs_L3[4][2][2] = 3*kappa1;
    Curvature_Higgs_L3[4][3][3] = 3*kappa1;
    Curvature_Higgs_L3[4][4][4] = 2*kappa;

    Curvature_Higgs_L4[0][0][0][0] = 3*lambda;
    Curvature_Higgs_L4[0][0][1][1] = 6*lambda;
    Curvature_Higgs_L4[0][0][2][2] = 6*lambda;
    Curvature_Higgs_L4[0][0][3][3] = 6*lambda;
    Curvature_Higgs_L4[0][0][4][4] = 6*kappa2;
    Curvature_Higgs_L4[0][1][0][1] = 6*lambda;
    Curvature_Higgs_L4[0][1][1][0] = 6*lambda;
    Curvature_Higgs_L4[0][2][0][2] = 6*lambda;
    Curvature_Higgs_L4[0][2][2][0] = 6*lambda;
    Curvature_Higgs_L4[0][3][0][3] = 6*lambda;
    Curvature_Higgs_L4[0][3][3][0] = 6*lambda;
    Curvature_Higgs_L4[0][4][0][4] = 6*kappa2;
    Curvature_Higgs_L4[0][4][4][0] = 6*kappa2;
    Curvature_Higgs_L4[1][0][0][1] = 6*lambda;
    Curvature_Higgs_L4[1][0][1][0] = 6*lambda;
    Curvature_Higgs_L4[1][1][0][0] = 6*lambda;
    Curvature_Higgs_L4[1][1][1][1] = 3*lambda;
    Curvature_Higgs_L4[1][1][2][2] = 6*lambda;
    Curvature_Higgs_L4[1][1][3][3] = 6*lambda;
    Curvature_Higgs_L4[1][1][4][4] = 6*kappa2;
    Curvature_Higgs_L4[1][2][1][2] = 6*lambda;
    Curvature_Higgs_L4[1][2][2][1] = 6*lambda;
    Curvature_Higgs_L4[1][3][1][3] = 6*lambda;
    Curvature_Higgs_L4[1][3][3][1] = 6*lambda;
    Curvature_Higgs_L4[1][4][1][4] = 6*kappa2;
    Curvature_Higgs_L4[1][4][4][1] = 6*kappa2;
    Curvature_Higgs_L4[2][0][0][2] = 6*lambda;
    Curvature_Higgs_L4[2][0][2][0] = 6*lambda;
    Curvature_Higgs_L4[2][1][1][2] = 6*lambda;
    Curvature_Higgs_L4[2][1][2][1] = 6*lambda;
    Curvature_Higgs_L4[2][2][0][0] = 6*lambda;
    Curvature_Higgs_L4[2][2][1][1] = 6*lambda;
    Curvature_Higgs_L4[2][2][2][2] = 3*lambda;
    Curvature_Higgs_L4[2][2][3][3] = 6*lambda;
    Curvature_Higgs_L4[2][2][4][4] = 6*kappa2;
    Curvature_Higgs_L4[2][3][2][3] = 6*lambda;
    Curvature_Higgs_L4[2][3][3][2] = 6*lambda;
    Curvature_Higgs_L4[2][4][2][4] = 6*kappa2;
    Curvature_Higgs_L4[2][4][4][2] = 6*kappa2;
    Curvature_Higgs_L4[3][0][0][3] = 6*lambda;
    Curvature_Higgs_L4[3][0][3][0] = 6*lambda;
    Curvature_Higgs_L4[3][1][1][3] = 6*lambda;
    Curvature_Higgs_L4[3][1][3][1] = 6*lambda;
    Curvature_Higgs_L4[3][2][2][3] = 6*lambda;
    Curvature_Higgs_L4[3][2][3][2] = 6*lambda;
    Curvature_Higgs_L4[3][3][0][0] = 6*lambda;
    Curvature_Higgs_L4[3][3][1][1] = 6*lambda;
    Curvature_Higgs_L4[3][3][2][2] = 6*lambda;
    Curvature_Higgs_L4[3][3][3][3] = 3*lambda;
    Curvature_Higgs_L4[3][3][4][4] = 6*kappa2;
    Curvature_Higgs_L4[3][4][3][4] = 6*kappa2;
    Curvature_Higgs_L4[3][4][4][3] = 6*kappa2;
    Curvature_Higgs_L4[4][0][0][4] = 6*kappa2;
    Curvature_Higgs_L4[4][0][4][0] = 6*kappa2;
    Curvature_Higgs_L4[4][1][1][4] = 6*kappa2;
    Curvature_Higgs_L4[4][1][4][1] = 6*kappa2;
    Curvature_Higgs_L4[4][2][2][4] = 6*kappa2;
    Curvature_Higgs_L4[4][2][4][2] = 6*kappa2;
    Curvature_Higgs_L4[4][3][3][4] = 6*kappa2;
    Curvature_Higgs_L4[4][3][4][3] = 6*kappa2;
    Curvature_Higgs_L4[4][4][0][0] = 6*kappa2;
    Curvature_Higgs_L4[4][4][1][1] = 6*kappa2;
    Curvature_Higgs_L4[4][4][2][2] = 6*kappa2;
    Curvature_Higgs_L4[4][4][3][3] = 6*kappa2;
    Curvature_Higgs_L4[4][4][4][4] = 12*lambdas;  

    SetCurvatureDone=true;

}


void Class_Potential_SMRSS::MinimizeOrderVEV(const std::vector<double>& vevMinimizer, std::vector<double>& vevFunction){
  /*
   * Here you do the conversion from the result vector from the minimizer (which has nVEV entries) to a
   * vector with NHiggs entries which you can use to call the member functions
   */

	VevOrder.resize(nVEV);
	// Here you have to tell which scalar field gets which VEV.
	VevOrder[0] = 2;
	VevOrder[1] = 4;

	int count=0;
	if(vevFunction.size() != 0)
	{
		for(int i=0;i<NHiggs;i++)
			{
				if(i==VevOrder[count]) {
					vevFunction[i] =  vevMinimizer.at(count);
					count++;
				}
				else vevFunction[i] = 0;

			}
	 }
	 else{
		 for(int i=0;i<NHiggs;i++)
		 {
			 if(i==VevOrder[count]) {
				 vevFunction.push_back(vevMinimizer.at(count));
				 count++;
			 }
			 else vevFunction.push_back(0);
		 }
	}
}

bool Class_Potential_SMRSS::CalculateDebyeSimplified(){
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass matrix and implement
   * your formula here and return true. The vector is given by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_Potential_SMRSS::CalculateDebyeGaugeSimplified()
{
    return false;
}
double Class_Potential_SMRSS::VTreeSimplified(const std::vector<double>& v){
	UseVTreeSimplified = false;
	return 0;
}

double Class_Potential_SMRSS::VCounterSimplified(const std::vector<double>& v)
{
	UseVCounterSimplified = false;
	return 0;
}

void Class_Potential_SMRSS::Debugging(const std::vector<double>& input, std::vector<double>& output){

}
