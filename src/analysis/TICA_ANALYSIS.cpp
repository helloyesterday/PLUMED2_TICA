/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
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
#include "Analysis.h"
#include "tools/Matrix.h"
#include "tools/IFile.h"
#include "reference/Direction.h"
#include "reference/MetricRegister.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/ReferenceValuePack.h"
#include "core/ActionRegister.h"

//+PLUMEDOC DIMRED TICA
/* 
Perform principal component analysis (TICA) using either the positions of the atoms a large number of collective variables as input.

Principal component analysis is a statistical technique that uses an orthogonal transformation to convert a set of observations of 
poorly correlated variables into a set of linearly uncorrelated variables.  You can read more about the specifics of this technique
here: https://en.wikipedia.org/wiki/Principal_component_analysis

When used with molecular dynamics simulations a set of frames taken from the trajectory, \f$\{X_i\}\f$, or the values of 
a number of collective variables which are calculated from the trajectory frames are used as input.  In this second instance your 
input to the TICA analysis algorithm is thus a set of high-dimensional vectors of collective variables.  However, if
collective variables are calculated from the positions of the atoms or if the positions are used directly the assumption is that 
this input trajectory is a set of poorly correlated (high-dimensional) vectors.  After principal component analysis has been 
performed the output is a set of orthogonal vectors that describe the directions in which the largest motions have been seen.  
In other words, principal component analysis provides a method for lowering the dimensionality of the data contained in a trajectory.
These output directions are some linear combination of the \f$x\f$, \f$y\f$ and \f$z\f$ positions if the positions were used as input 
or some linear combination of the input collective variables if a high-dimensional vector of collective variables was used as input.

As explained on the Wikipedia page you must calculate the average and covariance for each of the input coordinates.  In other words, you must 
calculate the average structure and the amount the system fluctuates around this average structure.  The problem in doing so when the 
\f$x\f$, \f$y\f$ and \f$z\f$ coordinates of a molecule are used as input is that the majority of the changes in the positions of the 
atoms comes from the translational and rotational degrees of freedom of the molecule.  The first six principal components will thus, most likely,
be uninteresting.  Consequently, to remedy this problem PLUMED provides the functionality to perform an RMSD alignment of the all the structures 
to be analysed to the first frame in the trajectory.  This can be used to effectively remove translational and/or rotational motions from 
consideration.  The resulting principal components thus describe vibrational motions of the molecule. 

If you wish to calculate the projection of a trajectory on a set of principal components calculated from this TICA action then the output can be 
used as input for the \ref TICAVARS action.

\par Examples

The following input instructs PLUMED to perform a principal component analysis in which the covariance matrix is calculated from changes in the positions
of the first 22 atoms.  The TYPE=OPTIMAL instruction ensures that translational and rotational degrees of freedom are removed from consideration.
The first two principal components will be output to a file called pca-comp.pdb.  Trajectory frames will be collected on every step and the TICA calculation
will be performed at the end of the simulation.

\verbatim
TICA METRIC=OPTIMAL ATOMS=1-22 STRIDE=1 USE_ALL_DATA NLOW_DIM=2 OFILE=pca-comp.pdb
\endverbatim

The following input instructs PLUMED to perform a principal component analysis in which the covariance matrix is calculated from chnages in the six distances
seen in the previous lines.  Notice that here the TYPE=EUCLIDEAN keyword is used to indicate that no alighment has to be done when calculating the various
elements of the covariance matrix from the input vectors.  In this calculation the first two principal components will be output to a file called pca-comp.pdb.
Trajectory frames will be collected every five steps and the TICA calculation is performed every 1000 steps.  Consequently, if you run a 2000 step simulation the 
TICA analysis will be performed twice.  The REWEIGHT_BIAS keyword in this input tells PLUMED that rather that ascribing a weight of one to each of the frames
when calculating averages and covariances a reweighting should be performed based and each frames' weight in these calculations should be determined based on 
the current value of the instantaneous bias (see \ref REWEIGHT_BIAS).  

\verbatim
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,3
d3: DISTANCE ATOMS=1,4
d4: DISTNACE ATOMS=2,3
d5: DISTANCE ATOMS=2,4
d6: DISTANCE ATOMS=3,4

TICA ARG=d1,d2,d3,d4,d5,d6 METRIC=EUCLIDEAN STRIDE=5 RUN=1000 NLOW_DIM=2 REWEIGHT_BIAS OFILE=pca-comp.pdb
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
	double dt;
	double rnorm;
	
	bool is_rescaled;
	bool is_read_corr;
	bool is_debug;
	
	std::string eigval_file;
	std::string eigvec_file;
	std::string corr_file;
	std::string corr_input;
	std::string rescale_file;
	std::string debug_file;
	
	OFile oeigval;
	OFile oeigvec;
	OFile ocorr;
	OFile orescale;
	OFile odebug;
	
	IFile icorr;
	
	std::vector<double> myaverages;
	std::vector<std::vector<double> > mydatas;
	std::vector<std::vector<double> > Rnew;
	std::vector<double> Wnew;
	std::vector<double> mywts;
	std::vector<bool> nandata;
	std::vector<std::vector<double>::size_type> effid;
	std::vector<std::string> strcolid;
	std::vector<double> atau;
	std::vector<Matrix<double> > Cread;
	
	Matrix<double> Czero;
	double get_correlation_matrix_element(const std::vector<double>& X,const std::vector<double>& Y,const std::vector<bool>& ynan,unsigned _tau);
public:
	static void registerKeywords( Keywords& keys );
	explicit TICA(const ActionOptions&ao);
	~TICA();
	void performAnalysis();
	void performTask( const unsigned& , const unsigned& , MultiValue& ) const { plumed_error(); }
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
	keys.add("optional","READ_CORR_FILE","read the correlations from file");
	keys.add("optional","STEP_SIZE","the simulation time step size");
	keys.add("optional","EIGEN_VALUE_FILE","the file to output the eigen value");
	keys.add("optional","RESCALE_FILE","the file that output rescaled trajectory");
	keys.add("optional","DEBUG_FILE","the file that debug information");
}

TICA::TICA(const ActionOptions&ao):
PLUMED_ANALYSIS_INIT(ao),
narg(0),delta_tau(0),npoints(0),tot_steps(0),wnorm(0),dt(1.0),rnorm(0.0),
is_rescaled(false),is_read_corr(false),is_debug(false)
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
	
	for(unsigned i=0;i!=narg;++i)
	{
		std::string id;
		Tools::convert(i,id);
		std::string colid="COL"+id;
		strcolid.push_back(colid);
	}
	
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
	parse("STEP_SIZE",dt);
	
	parse("RESCALE_FILE",rescale_file);
	if(rescale_file.size()>0)
	{
		is_rescaled=true;
		orescale.link(*this);
		orescale.open(rescale_file.c_str());
	}
	parse("DEBUG_FILE",debug_file);
	if(debug_file.size()>0)
	{
		is_debug=true;
		odebug.link(*this);
		odebug.open(debug_file.c_str());
		odebug.fmtField(" %f");
	}
	parse("READ_CORR_FILE",corr_input);
	if(corr_input.size()>0)
	{
		is_read_corr=true;
		icorr.link(*this);
		icorr.open(corr_input.c_str());
		icorr.allowIgnoredFields();
		
		int _narg=0;
		if(icorr.FieldExist("DIMENSION"))
			icorr.scanField("DIMENSION",_narg);
		plumed_massert(unsigned(_narg)==narg,"the number of CVs read from the correlation file mismatch!");
		
		int _npoints=0;
		if(icorr.FieldExist("NPOINTS"))
			icorr.scanField("NPOINTS",_npoints);
		plumed_massert(unsigned(_npoints)==npoints,"the number of points read from the correlation file mismatch!");
		
		if(icorr.FieldExist("TOTAL_WEIGHT"))
			icorr.scanField("TOTAL_WEIGHT",rnorm);
			
		for(unsigned ip=0;ip!=npoints;++ip)
		{
			int tt=0;
			icorr.scanField("LAGGED_TIME",tt);
			plumed_massert(unsigned(tt)==ip*delta_tau,"the lagged time read from the correlation file mismatch!");
			atau.push_back(unsigned(tt));
			Matrix<double> cc(narg,narg);
			for(unsigned i=0;i!=narg;++i)
			{
				int rr;
				icorr.scanField("ROW",rr);
				plumed_massert(unsigned(rr)==i,"the number of rows read from the correlation file mismatch!");
				for(unsigned j=0;j!=narg;++j)
				{
					double mm;
					icorr.scanField(strcolid[j],mm);
					cc[i][j]=mm;
				}
				icorr.scanField();
			}
			Cread.push_back(cc);
		}
		icorr.close();
	}
	
	ocorr.link(*this);
	ocorr.open(corr_file.c_str());
	
	oeigvec.link(*this);
	oeigvec.open(eigvec_file.c_str());
	oeigvec.addConstantField("COMPONENT");
	
	if(eigval_file.size()>0)
	{
		oeigval.link(*this);
		oeigval.open(eigval_file.c_str());
	}

	checkRead();
	
	log.printf("  with %d arguments.\n",narg);
	log.printf("  with delta lagged time: %d.\n",delta_tau);
	log.printf("  with time step size: %f.\n",dt);
	log.printf("  with number of points: %d.\n",npoints);
	if(eigval_file.size()>0)
		log.printf("  with eigen values output file: %s\n",eigval_file.c_str());
	if(is_rescaled)
		log.printf("  with rescaled output file: %s\n",rescale_file.c_str());
	if(is_debug)
		log.printf("  with debug file: %s\n",debug_file.c_str());
	if(is_read_corr)
	{
		log.printf("  read correlation from file: %s.\n",corr_input.c_str());
		log.printf("  with total weight: %f\n",rnorm);
		log.printf("  with lagged times:\n");
		for(unsigned i=0;i!=atau.size();++i)
			log.printf("    %d\t%f\n",int(i),atau[i]*dt);
	}
}

TICA::~TICA()
{
	ocorr.close();
	oeigvec.close();
	if(eigval_file.size()>0)
		oeigval.close();
	if(is_rescaled)
		orescale.close();
	if(is_debug)
		odebug.close();
}

void TICA::performAnalysis()
{
	tot_steps=getNumberOfDataPoints();

	for(unsigned j=0;j!=tot_steps;++j)
	{
		std::vector<double> points(narg);
		double ww;
		getDataPoint(j,points,ww);
		mywts.push_back(ww);
		wnorm+=ww;
		if(is_debug)
			odebug.printField("Index",int(j));
		for(unsigned i=0;i!=narg;++i)
		{
			mydatas[i].push_back(points[i]);
			myaverages[i]+=points[i]*ww;
			if(is_debug)
			{
				std::string id;
				Tools::convert(i,id);
				std::string strid="CV"+id;
				odebug.printField(strid,points[i]);
			}
		}
		if(is_debug)
		{
			odebug.printField("Weight",ww);
			odebug.printField();
		}
	}
	if(is_debug)
		odebug.flush();
	for(unsigned i=0;i!=narg;++i)
	{
		myaverages[i]/=wnorm;
		for(unsigned j=0;j!=tot_steps;++j)
			mydatas[i][j]-=myaverages[i];
	}
	log.printf("Fininished reading.\n");
	log.printf("  with total steps: %d.\n",tot_steps);
	
	double cr0,cr1;
	if(is_read_corr)
	{
		cr0=rnorm/(rnorm+wnorm);
		cr1=wnorm/(rnorm+wnorm);
	}

	std::vector<double> weights2(mywts);
	weights2.insert(weights2.begin(),0);
	std::vector<double> neff(weights2.size());
	std::partial_sum(weights2.begin(),weights2.end(),neff.begin());
	neff.pop_back();

	unsigned nsize=neff.size();
	unsigned rwsize=floor(neff.back()+0.5)+1;
	log.printf("  with size of the array: %d.\n",nsize);
	log.printf("  with reweight time steps: %d.\n",rwsize);

	Rnew.assign(narg,std::vector<double>(rwsize));
	Wnew.assign(rwsize,0);
	nandata.assign(rwsize,true);
	
	for(unsigned i=0;i!=nsize;++i)
	{
		unsigned s=int(neff[i]);
		effid.push_back(s);
		for(unsigned j=0;j!=narg;++j)
			Rnew[j][s]=mydatas[j][i];
		Wnew[s]=mywts[i];
		nandata[s]=false;
	}

	for(unsigned i=0;i!=narg;++i)
		Rnew[i].erase(Rnew[i].begin());
	Wnew.erase(Wnew.begin());
	
	if(is_rescaled)
	{
		orescale.fmtField(" %f");
		for(unsigned i=0;i!=rwsize;++i)
		{
			orescale.printField("Index",int(i));
			orescale.printField("NEFF",neff[i]);
			for(unsigned j=0;j!=narg;++j)
			{
				std::string id;
				Tools::convert(j,id);
				std::string strid="CV"+id;
				orescale.printField(strid,Rnew[j][i]);
			}
			orescale.printField("Weight",Wnew[i]);
			orescale.printField();
		}
		orescale.flush();
	}

	Czero=build_correlation(0);

	ocorr.fmtField(" %f");
	ocorr.addConstantField("DIMENSION");
	ocorr.printField("DIMENSION",int(narg));
	ocorr.addConstantField("NPOINTS");
	ocorr.printField("NPOINTS",int(npoints));
	ocorr.addConstantField("TOTAL_WEIGHT");
	ocorr.printField("TOTAL_WEIGHT",wnorm+rnorm);
	ocorr.addConstantField("LAGGED_TIME");
	
	//~ ocorr.printField("LAGGED_TIME",0);
	//~ for(unsigned i=0;i!=narg;++i)
	//~ {
		//~ ocorr.printField("ROW",int(i));
		//~ for(unsigned j=0;j!=narg;++j)
		//~ {
			//~ std::string id;
			//~ Tools::convert(j,id);
			//~ std::string colid="COL"+id;
			//~ ocorr.printField(colid,Czero[i][j]);
		//~ }
		//~ ocorr.printField();
	//~ }
	//~ ocorr.flush();
	
	std::vector<double> wt0;
	Matrix<double> vt0;
	
	diagMat(Czero,wt0,vt0);

	std::multimap<double,std::vector<double> > eigs0;

	for(unsigned i=0;i!=narg;++i)
	{
		std::vector<double> rowvec;
		for(unsigned j=0;j!=narg;++j)
			rowvec.push_back(vt0[i][j]);
		//~ eigs0[wt0[i]]=vt0.get_row(i);
		//~ eigs0[wt0[i]]=rowvec;
		eigs0.insert(std::pair<double,std::vector<double> >(wt0[i],rowvec));
	}

	std::vector<double> eigval0;
	std::vector<std::vector<double> > eigvec0;
	for(std::multimap<double,std::vector<double> >::reverse_iterator i=eigs0.rbegin();i!=eigs0.rend();++i)
	{
		eigval0.push_back(i->first);
		eigvec0.push_back(i->second);
	}
	log.printf("  with Eigenvalues of <C(0),C(0)>.\n");
	for(unsigned i=0;i!=narg;++i)
		log.printf("    %d\t%f\n",int(i),eigval0[i]);
	log.printf("\n");

	std::vector<std::vector<double> >& U=eigvec0;
	
	unsigned eff_size=0;
	for(unsigned i=0;i!=narg;++i)
		if(eigval0[i]>0)	++eff_size;
	
	std::vector<double> s_sqrt;
	for(unsigned i=0;i!=eff_size;++i)
		s_sqrt.push_back(1.0/sqrt(eigval0[i]));
	
	Matrix<double> X(narg,eff_size);
	for(unsigned i=0;i!=narg;++i)
	{
		for(unsigned j=0;j!=eff_size;++j)
			X[i][j]=U[j][i]*s_sqrt[j];
	}
	//~ Matrix<double> X_adj=X.trans();
	Matrix<double> X_adj;
	transpose(X,X_adj);
	
	log.printf("Running TICA ...\n");
	for(unsigned ip=0;ip!=npoints;++ip)
	{
		unsigned tau=ip*delta_tau;
		
		double time=tau*dt;
		plumed_massert(tau<rwsize,"the lagged time is too large!");
		
		Matrix<double> Clag=build_correlation(tau);
		if(is_read_corr)
		{
			for(unsigned i=0;i!=narg;++i)
			{
				for(unsigned j=0;j!=narg;++j)
					Clag[i][j]=cr0*Cread[ip][i][j]+cr1*Clag[i][j];
			}
		}
		
		ocorr.printField("LAGGED_TIME",int(tau));
		for(unsigned i=0;i!=narg;++i)
		{
			ocorr.printField("ROW",int(i));
			for(unsigned j=0;j!=narg;++j)
				ocorr.printField(strcolid[j],Clag[i][j]);
			ocorr.printField();
		}
		ocorr.flush();
		
		//~ Matrix<double> C_sym = 0.5*(Clag+Clag.trans());
		Matrix<double> Clag_T;
		transpose(Clag,Clag_T);
		Matrix<double> C_sym(narg,narg);
		for(unsigned i=0;i!=narg;++i)
		{
			for(unsigned j=0;j!=narg;++j)
				C_sym[i][j]=(Clag[i][j]+Clag_T[i][j])/2;
		}
		
		//~ Matrix<double> C_can = X_adj*(C_sym*X);
		Matrix<double> CsymX;
		mult(C_sym,X,CsymX);
		Matrix<double> C_can;
		mult(X_adj,CsymX,C_can);
		
		std::vector<double> nwt;
		Matrix<double> nvt;
		diagMat(C_can,nwt,nvt);
		
		std::multimap<double,std::vector<double> > eigs2;
				for(unsigned i=0;i!=nwt.size();++i)
		{
			//~ eigs2[nwt[i]]=nvt.get_row(i);
			std::vector<double> rowvec;
			for(unsigned j=0;j!=narg;++j)
				rowvec.push_back(nvt[i][j]);
			//~ eigs2[nwt[i]]=rowvec;
			eigs2.insert(std::pair<double,std::vector<double> >(nwt[i],rowvec));
		}

		std::vector<double> eigval2;
		std::vector<std::vector<double> > eigvec2;
		for(std::multimap<double,std::vector<double> >::reverse_iterator miter=eigs2.rbegin();miter!=eigs2.rend();++miter)
		{
			eigval2.push_back(miter->first);
			eigvec2.push_back(miter->second);
		}
		
		log.printf("  Setp %f:\n",time);
		log.printf("  Eigenvalues:");
		for(unsigned i=0;i!=narg;++i)
			log.printf(" %f",eigval2[i]);
		log.printf("\n\n");
		if(eigval_file.size()>0)
		{
			oeigval.fmtField(" %f");
			oeigval.printField("Time",time);
			for(unsigned i=0;i!=narg;++i)
			{
				std::string id;
				Tools::convert(i,id);
				std::string eigid="EIGVAL"+id;
				oeigval.printField(eigid,eigval2[i]);
			}
			oeigval.printField();
			oeigval.flush();
		}
		
		Matrix<double> eigvecs(narg,narg);
		for(unsigned i=0;i!=narg;++i)
		{
			for(unsigned j=0;j!=narg;++j)
				eigvecs[i][j]=eigvec2[i][j];
		}
		std::vector<double> meigvec=eigvec2.front();
		
		//~ Matrix<double> mvt=X*eigvecs.trans();
		Matrix<double> eigvecs_T;
		transpose(eigvecs,eigvecs_T);
		Matrix<double> mvt;
		mult(X,eigvecs_T,mvt);
		
		if(mvt[0][comp]<0)
		{
			//~ mvt*=-1.0;
			for(unsigned i=0;i!=narg;++i)
			{
				for(unsigned j=0;j!=narg;++j)
					mvt[i][j]*=-1;
			}
		}
		
		std::vector<std::vector<double> > tica_res;
		oeigvec.fmtField(" %f");
		for(unsigned k=0;k!=comp+1;++k)
		{
			//~ std::vector<double> evec=mvt.get_column(k);
			std::vector<double> evec;
			for(unsigned j=0;j!=narg;++j)
				evec.push_back(mvt[j][k]);
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
	std::vector<std::vector<double> > Y(Rnew);
	std::vector<bool> ynan;
	if(_tau!=0)
	{
		for(unsigned i=0;i!=narg;++i)
			Y[i].erase(Y[i].begin(),Y[i].begin()+_tau);
		ynan=nandata;
		ynan.erase(ynan.begin(),ynan.begin()+_tau);
	}
	
	Matrix<double> corrmat(narg,narg);
	for(unsigned i=0;i!=narg;++i)
	{
		if(_tau==0)
		{
			for(unsigned j=0;j!=narg;++j)
				corrmat[i][j]=get_correlation_matrix_element(Rnew[i],Y[j],nandata,_tau);
		}
		else
		{
			for(unsigned j=0;j!=narg;++j)
				corrmat[i][j]=get_correlation_matrix_element(Rnew[i],Y[j],ynan,_tau);
		}
	}
	return corrmat;
}

double TICA::get_correlation_matrix_element(const std::vector<double>& X,const std::vector<double>& Y,const std::vector<bool>& ynan,unsigned _tau)
{
	double corrmean=0;
	double norm=0;
	if(_tau==0)
	{
		for(unsigned i=0;i!=effid.size();++i)
		{
			std::vector<double>::size_type id=effid[i];
			double value=X[id]*Y[id];
			corrmean+=value*Wnew[id];
			norm+=Wnew[id];
		}
	}
	else
	{
		for(unsigned i=0;i!=effid.size();++i)
		{
			std::vector<double>::size_type id=effid[i];
			if(id>=Y.size())
				break;
			if(ynan[id])
				continue;
			double value=X[id]*Y[id];
			corrmean+=value*Wnew[id];
			norm+=Wnew[id];
		}
	}
	corrmean/=norm;

	return corrmean;
}


}
}
