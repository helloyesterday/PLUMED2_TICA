/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
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
#include <numeric>
#include <iterator>
#include "Analysis.h"
#include "tools/Matrix.h"
#include "reference/Direction.h"
#include "reference/MetricRegister.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/ReferenceValuePack.h"
#include "core/ActionRegister.h"

//+PLUMEDOC GRIDCALC AVERAGE 
/* 
\endverbatim

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class TICA : public Analysis {
private:
	unsigned narg;
	unsigned delta_tau;
	unsigned npoints;
	unsigned tot_steps;
	unsigned rwsize;
	unsigned eff_size;
	unsigned comp;
	
	double wnorm;
	
	std::string eigval_file;
	std::string eigvec_file;
	std::string corr_file;
	
	OFile oeigval;
	OFile oeigvec;
	OFile ocorr;
	
	std::vector<double> myaverages;
	std::vector<std::vector<double> > mydatas;
	std::vector<std::vector<double> > Rnew;
	std::vector<double> Wnew;
	std::vector<double> mywts;
	
	Matrix<double> Czero;
	double get_correlation_matrix_element(const std::vector<double>& X,const std::vector<double>& Y,const std::vector<double>& wts,unsigned _tau);
public:
	static void registerKeywords( Keywords& keys );
	explicit TICA(const ActionOptions&ao);
	~TICA();
	void performAnalysis();
	void accumulate();
	Matrix<double> build_correlation(unsigned _tau=0);
};

PLUMED_REGISTER_ACTION(TICA,"TICA")

void TICA::registerKeywords( Keywords& keys ){
	Analysis::registerKeywords( keys );
  	keys.add("compulsory","DELTA_TAU","the delta lagged time to calculate TICA");
	keys.add("compulsory","POINTS_NUMBER","how much point to output");
	keys.add("compulsory","EIGEN_NUMBER","0","which eigenvector to print to eigenvector output file (0 corresponds to largest eigenvalue)");
	keys.add("compulsory","OUTPUT_FILE","TICA.dat","the file to output the result");
	keys.add("compulsory","CORRELATION_FILE","correlation.dat","the file to output the correlation matrix");
	
	keys.add("optional","EIGEN_VALUE_FILE","the file to output the eigen value");
}

TICA::TICA(const ActionOptions&ao):
PLUMED_ANALYSIS_INIT(ao),
narg(0),delta_tau(0),npoints(0),tot_steps(0),wnorm(0)
{
	addValue(); // Create a value so that we can output the average
	narg=getNumberOfArguments();
	std::string instring; 
	if( getPntrToArgument(0)->isPeriodic() ){
	   std::string min, max; getPntrToArgument(0)->getDomain(min,max); 
	   instring = "PERIODIC=" + min + "," + max; setPeriodic( min, max );
	} else {
	   setNotPeriodic();
	}
	myaverages.resize(narg);
	mydatas.resize(narg);
	
	parse("DELTA_TAU",delta_tau);
	if(delta_tau==0)
		error("DELTA_TAU could not be 0");
	parse("POINTS_NUMBER",npoints);
	if(npoints==0)
		error("POINTS_NUMBER could not be 0");
		
	parse("EIGEN_NUMBER",comp);
	parse("CORRELATION_FILE",corr_file);
	parse("OUTPUT_FILE",eigvec_file);
	parse("EIGEN_VALUE_FILE",eigval_file);
	
	ocorr.link(*this);
	ocorr.open(corr_file.c_str());
	
	oeigvec.link(*this);
	oeigvec.open(eigvec_file.c_str());
	oeigvec.addConstantField("COMPONENT");
	
	if(eigval_file.size()>0)
	{
		oeigval.link(*this);
		oeigval.open(eigvec_file.c_str());
	}

	checkRead();
	
	log.printf("  with %d arguments.\n",narg);
	log.printf("  with delta lagged time: %d.\n",delta_tau);
	log.printf("  with number of points: %d.\n",npoints);
}

TICA::~TICA()
{
	ocorr.close();
	oeigvec.close();
	if(eigval_file.size()>0)
		oeigval.close();
}

void TICA::accumulate()
{
	for(unsigned i=0;i!=narg;++i)
	{
		double arg_now=getArgument(i);
		myaverages[i]+=cweight*arg_now;
		mydatas[i].push_back(arg_now);
	}
	mywts.push_back(cweight);
	wnorm+=cweight;
	++tot_steps;
}

void TICA::performAnalysis()
{
	for(unsigned i=0;i!=narg;++i)
	{
		double arg_aver=myaverages[i]/wnorm;
		for(unsigned j=0;j!=tot_steps;++j)
			mydatas[i][j]-=arg_aver;
	}
	log.printf("  All data reading finished.\n");
	log.printf("  with total steps: %d.\n",tot_steps);

	std::vector<double> weights2(mywts);
	weights2.insert(weights2.begin(),1.0);
	std::vector<double> neff(weights2.size());
	std::partial_sum(weights2.begin(),weights2.end(),neff.begin());
	neff.pop_back();
	
	unsigned nsize=neff.size();
	unsigned rwsize=floor(neff.back()+0.5)+1;
	log.printf("  with size of the array: %d.\n",nsize);
	log.printf("  with reweight time steps: %d.\n",rwsize);
	
	Rnew.assign(narg,std::vector<double>(rwsize));
	Wnew.assign(rwsize,0);
	
	for(unsigned i=0;i!=nsize;++i)
	{
		unsigned s=int(neff[i]);
		for(unsigned j=0;j!=narg;++j)
			Rnew[j][s]=mydatas[j][i];
		Wnew[s]=mywts[i];
	}
		
	for(unsigned i=0;i!=narg;++i)
		Rnew[i].erase(Rnew[i].begin());
	Wnew.erase(Wnew.begin());

	Czero=build_correlation(0);
	
	for(unsigned i=0;i!=narg;++i)
	{
		ocorr.printField("ROW",int(i));
		for(unsigned j=0;j!=narg;++j)
		{
			std::string id;
			Tools::convert(i,id);
			std::string colid="COL"+id;
			ocorr.printField(colid,Czero[i][j]);
		}
		ocorr.printField();
	}
	ocorr.flush();
	
	std::vector<double> wt;
	Matrix<double> vt;
	
	diagMat(Czero,wt,vt);

	std::map<double,std::vector<double> > eigs;

	for(unsigned i=0;i!=wt.size();++i)
		eigs[wt[i]]=vt.get_row(i);

	std::vector<double> eigval;
	std::vector<std::vector<double> > eigvec;
	for(std::map<double,std::vector<double> >::const_iterator i=eigs.begin();i!=eigs.end();++i)
	{
		eigval.push_back(i->first);
		eigvec.push_back(i->second);
	}
	log.printf("  Eigenvalues of <C(0),C(0)>.\n");
	for(unsigned i=0;i!=narg;++i)
		log.printf("    %d\t%f\n",int(i),eigval[i]);
	log.printf("\n");
	
	std::vector<std::vector<double> >& U=eigvec;
	
	unsigned eff_size=0;
	for(unsigned i=0;i!=narg;++i)
		if(eigval[i]>0)	++eff_size;
	
	std::vector<double> s_sqrt;
	for(unsigned i=0;i!=eff_size;++i)
		s_sqrt.push_back(1.0/sqrt(eigval[i]));
	
	Matrix<double> X(narg,eff_size);
	for(unsigned i=0;i!=narg;++i)
	{
		for(unsigned j=0;j!=eff_size;++j)
			X[i][j]=U[j][i]*s_sqrt[j];
	}
	Matrix<double> X_adj=X.trans();
	
	log.printf("  Running TICA ...\n");
	for(unsigned ip=0;ip!=npoints;++ip)
	{
		unsigned tau=ip*delta_tau;
		plumed_massert(tau<rwsize,"the lagged time is too large!");
		double time=tau*delta_tau;
		
		Matrix<double> Clag=build_correlation(tau);
		Matrix<double> C_sym = 0.5*(Clag+Clag.trans());
		Matrix<double> C_can = X_adj*(C_sym*X);
		
		std::vector<double> nwt;
		Matrix<double> nvt;
		diagMat(C_can,nwt,nvt);
		
		std::map<double,std::vector<double> > eigs2;
		for(unsigned i=0;i!=nwt.size();++i)
			eigs2[nwt[i]]=nvt.get_row(i);
			
		std::vector<double> eigval2;
		std::vector<std::vector<double> > eigvec2;
		for(std::map<double,std::vector<double> >::iterator miter=eigs2.begin();miter!=eigs2.end();++miter)
		{
			eigval2.push_back(miter->first);
			eigvec2.push_back(miter->second);
		}
		log.printf("  setp %f: eigenvalues:",time);
		for(unsigned i=0;i!=narg;++i)
			log.printf(" %f",eigval[i]);
		log.printf("\n");
		
		Matrix<double> eigvecs(narg,narg);
		for(unsigned i=0;i!=narg;++i)
		{
			for(unsigned j=0;j!=narg;++j)
				eigvecs[i][j]=eigvec2[i][j];
		}
		std::vector<double> meigvec=eigvec2.front();
		
		Matrix<double> mvt=X*eigvecs.trans();
		if(mvt[0][comp]<0)
			mvt*=-1.0;
		
		std::vector<std::vector<double> > tica_res;
		for(unsigned k=0;k!=comp+1;++k)
		{
			std::vector<double> evec=mvt.get_column(k);
			double mnorm=sqrt(norm(evec));
			oeigvec.printField("COMPONENT",int(k));
			oeigvec.printField("TIME",time);
			for(unsigned i=0;i!=narg;++i)
			{
				evec[i]/=mnorm;
				std::string id;
				Tools::convert(i,id);
				std::string eigid="EIGVEC"+id;
				oeigvec.printField(eigid,evec[i]);
			}
			oeigvec.printField();
			tica_res.push_back(evec);
		}
		oeigvec.flush();

		for(unsigned i=0;i!=eff_size;++i)
			if(nwt[i]<0) nwt[i]=0;
	}
}

Matrix<double> TICA::build_correlation(unsigned _tau)
{
	std::vector<std::vector<double> > X(Rnew);
	std::vector<std::vector<double> > Y(Rnew);
	std::vector<double> wts=(Wnew);
	for(unsigned i=0;i!=narg;++i)
	{
		X[i].erase(X[i].end()-_tau,X[i].end());
		Y[i].erase(Y[i].begin(),Y[i].begin()+_tau);
	}
	wts.erase(wts.end()-_tau,wts.end());
	
	Matrix<double> corrmat(narg,narg);
	for(unsigned i=0;i!=narg;++i)
	{
		for(unsigned j=0;j!=narg;++j)
			corrmat[i][j]=get_correlation_matrix_element(X[i],Y[j],wts,_tau);
	}
	return corrmat;
}

double TICA::get_correlation_matrix_element(const std::vector<double>& X,const std::vector<double>& Y,const std::vector<double>& wts,unsigned _tau)
{
	double corrmean=0;
	double norm=0;
	for(unsigned i=0;i!=X.size();++i)
	{
		double value=X[i]*Y[i];
		double wt=wts[i];
		corrmean+=value*wt;
		norm+=wt;
	}
	corrmean/=norm;

	return corrmean;
}


}
}
