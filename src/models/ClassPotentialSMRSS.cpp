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

  nPar = 4; // number of parameters in the tree-Level Lagrangian
  nParCT = 7; // number of parameters in the counterterm potential

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
  std::string out = "dT\tdlambda\tdmsquared";
  return out;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and the critical temperature. Use this to
 * complement the legend of the given input file
 */
std::string Class_Potential_SMRSS::addLegendTemp(){
  std::string out = "T_c\tv_c";
  //out += "Your VEV order";
  out += "omega";
  return out;
}

/**
 * returns a string which tells the user the chronological order of the Triple Higgs couplings. Use this to
 * complement the legend of the given input file
 *
 */
std::string Class_Potential_SMRSS::addLegendTripleCouplings(){
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
	out = "omega";
  return out;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Potential_SMRSS::ReadAndSet(const std::string& linestr, std::vector<double>& par )
{
	std::stringstream ss(linestr);
	double tmp;

	double lmu_hs, lmu_s, llambda_s, llambda_m;

	for(int k=1;k<=4;k++)
	{
	      ss>>tmp;
	      if(k==1) lmu_hs = tmp;
	      else if (k==2) lmu_s = tmp;
	      else if (k==3) llambda_s = tmp;
	      else if (k==4) llambda_m = tmp;
	}
	par[0] = lmu_hs;
	par[1] = lmu_s;
	par[2] = llambda_s;
	par[3] = llambda_m;


	set_gen(par); // This you have to call so that everything will be set
	return ;
}


/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Potential_SMRSS::set_gen(const std::vector<double>& par) {

    // direct input params
    mu_hs = par[0];
    mu_s = par[1];
    lambda_s = par[2];
    lambda_m = par[3];

    // EWSB condition eliminates lambda_h:
    lambda_h = (2.0*mu_hs)/(std::pow(C_vev0,2));

    // This sets the renormalization scale
	scale = C_vev0;

	vevTreeMin.resize(nVEV);
	vevTree.resize(NHiggs);

	// Here you have to set the vector vevTreeMin. The vector vevTree will then be set by the function MinimizeOrderVEV
	vevTreeMin[0] = C_vev0;
	vevTreeMin[1] = std::sqrt(-1.0*(mu_s/lambda_s)); // from EWSB

	MinimizeOrderVEV(vevTreeMin,vevTree);
	if(!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the entries of Curvature_Higgs_CT_L1 to
 * Curvature_Higgs_CT_L4.
 */
void Class_Potential_SMRSS::set_CT_Pot_Par(const std::vector<double>& par){

    Dmu_hs = par[0];
    Dlambda_h = par[1];
    Dmu_s = par[2];
    Dlambda_s = par[3];
    Dlambda_m = par[4];
    DT3 = par[5];
    DT5 = par[6];

    Curvature_Higgs_CT_L1[2] = DT3;
    Curvature_Higgs_CT_L1[4] = DT5;

    Curvature_Higgs_CT_L2[0][0] = (Dmu_hs)/2;
    Curvature_Higgs_CT_L2[1][1] = (Dmu_hs)/2;
    Curvature_Higgs_CT_L2[2][2] = (Dmu_hs)/2;
    Curvature_Higgs_CT_L2[3][3] = (Dmu_hs)/2;
    Curvature_Higgs_CT_L2[4][4] = -(Dmu_s);

    Curvature_Higgs_CT_L4[0][0][0][0] = (-3*(Dlambda_h))/2;
    Curvature_Higgs_CT_L4[0][0][1][1] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[0][0][2][2] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[0][0][3][3] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[0][0][4][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[0][1][0][1] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[0][1][1][0] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[0][2][0][2] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[0][2][2][0] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[0][3][0][3] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[0][3][3][0] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[0][4][0][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[0][4][4][0] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[1][0][0][1] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[1][0][1][0] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[1][1][0][0] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[1][1][1][1] = (-3*(Dlambda_h))/2;
    Curvature_Higgs_CT_L4[1][1][2][2] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[1][1][3][3] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[1][1][4][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[1][2][1][2] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[1][2][2][1] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[1][3][1][3] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[1][3][3][1] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[1][4][1][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[1][4][4][1] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[2][0][0][2] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[2][0][2][0] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[2][1][1][2] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[2][1][2][1] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[2][2][0][0] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[2][2][1][1] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[2][2][2][2] = (-3*(Dlambda_h))/2;
    Curvature_Higgs_CT_L4[2][2][3][3] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[2][2][4][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[2][3][2][3] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[2][3][3][2] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[2][4][2][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[2][4][4][2] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[3][0][0][3] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[3][0][3][0] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[3][1][1][3] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[3][1][3][1] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[3][2][2][3] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[3][2][3][2] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[3][3][0][0] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[3][3][1][1] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[3][3][2][2] = -3*(Dlambda_h);
    Curvature_Higgs_CT_L4[3][3][3][3] = (-3*(Dlambda_h))/2;
    Curvature_Higgs_CT_L4[3][3][4][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[3][4][3][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[3][4][4][3] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][0][0][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][0][4][0] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][1][1][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][1][4][1] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][2][2][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][2][4][2] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][3][3][4] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][3][4][3] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][4][0][0] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][4][1][1] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][4][2][2] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][4][3][3] = -3*(Dlambda_m);
    Curvature_Higgs_CT_L4[4][4][4][4] = -6*(Dlambda_s);


}


/**
 * console output of all Parameters
 */
void Class_Potential_SMRSS::write() {

	std::cout << "The parameters are : " << std::endl;
	std::cout << "lambda = " << lambda << std::endl
			<< "\tm^2 = " << ms << std::endl;

	std::cout << "The counterterm parameters are : " << std::endl;
	std::cout << "dT = "<< dT << std::endl
			<< "dlambda = " << dlambda << std::endl
			<< "dm^2 = "<< dms << std::endl;

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

	VectorXd NablaWeinberg(NHiggs);
	MatrixXd HesseWeinberg(NHiggs,NHiggs),HiggsRot(NHiggs,NHiggs);
	for(int i=0;i<NHiggs;i++)
	{
		NablaWeinberg[i] = WeinbergNabla[i];
		for(int j=0;j<NHiggs;j++) HesseWeinberg(i,j) = WeinbergHesse.at(j*NHiggs+i);
	}

	// Here you have to use your formulas for the counterterm scheme
	double t = 0;
	dT = t;
	dlambda = 3.0*t/std::pow(C_vev0,3) + 3.0/std::pow(C_vev0,3) * NablaWeinberg(0) -3.0/std::pow(C_vev0,2) *HesseWeinberg(0,0);
	dms = -3.0/(2*std::pow(C_vev0,2)) *NablaWeinberg(0) + 1.0/2.0 *HesseWeinberg(0,0) -3.0*t/(2*C_vev0);

	par[0] = dT;
	par[1] = dlambda;
	par[2] = dms;

	set_CT_Pot_Par(par);

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

    Curvature_Higgs_L2[0][0] = (mu_hs)/2;
    Curvature_Higgs_L2[1][1] = (mu_hs)/2;
    Curvature_Higgs_L2[2][2] = (mu_hs)/2;
    Curvature_Higgs_L2[3][3] = (mu_hs)/2;
    Curvature_Higgs_L2[4][4] = -(mu_s);

    Curvature_Higgs_L4[0][0][0][0] = (-3*(lambda_h))/2;
    Curvature_Higgs_L4[0][0][1][1] = -3*(lambda_h);
    Curvature_Higgs_L4[0][0][2][2] = -3*(lambda_h);
    Curvature_Higgs_L4[0][0][3][3] = -3*(lambda_h);
    Curvature_Higgs_L4[0][0][4][4] = -3*(lambda_m);
    Curvature_Higgs_L4[0][1][0][1] = -3*(lambda_h);
    Curvature_Higgs_L4[0][1][1][0] = -3*(lambda_h);
    Curvature_Higgs_L4[0][2][0][2] = -3*(lambda_h);
    Curvature_Higgs_L4[0][2][2][0] = -3*(lambda_h);
    Curvature_Higgs_L4[0][3][0][3] = -3*(lambda_h);
    Curvature_Higgs_L4[0][3][3][0] = -3*(lambda_h);
    Curvature_Higgs_L4[0][4][0][4] = -3*(lambda_m);
    Curvature_Higgs_L4[0][4][4][0] = -3*(lambda_m);
    Curvature_Higgs_L4[1][0][0][1] = -3*(lambda_h);
    Curvature_Higgs_L4[1][0][1][0] = -3*(lambda_h);
    Curvature_Higgs_L4[1][1][0][0] = -3*(lambda_h);
    Curvature_Higgs_L4[1][1][1][1] = (-3*(lambda_h))/2;
    Curvature_Higgs_L4[1][1][2][2] = -3*(lambda_h);
    Curvature_Higgs_L4[1][1][3][3] = -3*(lambda_h);
    Curvature_Higgs_L4[1][1][4][4] = -3*(lambda_m);
    Curvature_Higgs_L4[1][2][1][2] = -3*(lambda_h);
    Curvature_Higgs_L4[1][2][2][1] = -3*(lambda_h);
    Curvature_Higgs_L4[1][3][1][3] = -3*(lambda_h);
    Curvature_Higgs_L4[1][3][3][1] = -3*(lambda_h);
    Curvature_Higgs_L4[1][4][1][4] = -3*(lambda_m);
    Curvature_Higgs_L4[1][4][4][1] = -3*(lambda_m);
    Curvature_Higgs_L4[2][0][0][2] = -3*(lambda_h);
    Curvature_Higgs_L4[2][0][2][0] = -3*(lambda_h);
    Curvature_Higgs_L4[2][1][1][2] = -3*(lambda_h);
    Curvature_Higgs_L4[2][1][2][1] = -3*(lambda_h);
    Curvature_Higgs_L4[2][2][0][0] = -3*(lambda_h);
    Curvature_Higgs_L4[2][2][1][1] = -3*(lambda_h);
    Curvature_Higgs_L4[2][2][2][2] = (-3*(lambda_h))/2;
    Curvature_Higgs_L4[2][2][3][3] = -3*(lambda_h);
    Curvature_Higgs_L4[2][2][4][4] = -3*(lambda_m);
    Curvature_Higgs_L4[2][3][2][3] = -3*(lambda_h);
    Curvature_Higgs_L4[2][3][3][2] = -3*(lambda_h);
    Curvature_Higgs_L4[2][4][2][4] = -3*(lambda_m);
    Curvature_Higgs_L4[2][4][4][2] = -3*(lambda_m);
    Curvature_Higgs_L4[3][0][0][3] = -3*(lambda_h);
    Curvature_Higgs_L4[3][0][3][0] = -3*(lambda_h);
    Curvature_Higgs_L4[3][1][1][3] = -3*(lambda_h);
    Curvature_Higgs_L4[3][1][3][1] = -3*(lambda_h);
    Curvature_Higgs_L4[3][2][2][3] = -3*(lambda_h);
    Curvature_Higgs_L4[3][2][3][2] = -3*(lambda_h);
    Curvature_Higgs_L4[3][3][0][0] = -3*(lambda_h);
    Curvature_Higgs_L4[3][3][1][1] = -3*(lambda_h);
    Curvature_Higgs_L4[3][3][2][2] = -3*(lambda_h);
    Curvature_Higgs_L4[3][3][3][3] = (-3*(lambda_h))/2;
    Curvature_Higgs_L4[3][3][4][4] = -3*(lambda_m);
    Curvature_Higgs_L4[3][4][3][4] = -3*(lambda_m);
    Curvature_Higgs_L4[3][4][4][3] = -3*(lambda_m);
    Curvature_Higgs_L4[4][0][0][4] = -3*(lambda_m);
    Curvature_Higgs_L4[4][0][4][0] = -3*(lambda_m);
    Curvature_Higgs_L4[4][1][1][4] = -3*(lambda_m);
    Curvature_Higgs_L4[4][1][4][1] = -3*(lambda_m);
    Curvature_Higgs_L4[4][2][2][4] = -3*(lambda_m);
    Curvature_Higgs_L4[4][2][4][2] = -3*(lambda_m);
    Curvature_Higgs_L4[4][3][3][4] = -3*(lambda_m);
    Curvature_Higgs_L4[4][3][4][3] = -3*(lambda_m);
    Curvature_Higgs_L4[4][4][0][0] = -3*(lambda_m);
    Curvature_Higgs_L4[4][4][1][1] = -3*(lambda_m);
    Curvature_Higgs_L4[4][4][2][2] = -3*(lambda_m);
    Curvature_Higgs_L4[4][4][3][3] = -3*(lambda_m);
    Curvature_Higgs_L4[4][4][4][4] = -6*(lambda_s);

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
	bool Debug = false;
  if(Debug) std::cout << "Debug turned on in Class_Template :: " << __func__ << std::endl;

  /*
     * Use this function if you calculated the Debye corrections to the gauge mass matrix and implement
     * your formula here and return true. The vector is given by DebyeGauge[NGauge][NGauge]
     */


  return false;
}
double Class_Potential_SMRSS::VTreeSimplified(const std::vector<double>& v){
	UseVTreeSimplified = false;
	UseVTreeSimplified = true; // To use the simplified version
	double res = 0;

	double vIn = v[0];
	res = 0.5*ms*std::pow(vIn,2) + 1.0/24.0 * lambda * std::pow(vIn,4);

	return res;
}

double Class_Potential_SMRSS::VCounterSimplified(const std::vector<double>& v)
{
	UseVCounterSimplified = false;
	UseVCounterSimplified = true; // To use the simplified version
	if(not UseVCounterSimplified) return 0;
	double res = 0;
	double vIn = v[0];
		res = 0.5*dms*std::pow(vIn,2) + 1.0/24.0 * dlambda * std::pow(vIn,4) + dT * vIn;
	return res;
}

void Class_Potential_SMRSS::Debugging(const std::vector<double>& input, std::vector<double>& output){

}
