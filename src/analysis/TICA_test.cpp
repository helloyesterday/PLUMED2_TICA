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
#include "vesselbase/ActionWithAveraging.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "AverageVessel.h"
#include "tools/Matrix.h"
#include "tools/IFile.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/ReferenceValuePack.h"

//+PLUMEDOC DIMRED TICA_TEST 
/* 
Time-lagged independent component analysis (TICA_TEST) using a large number of collective variables as input.


\par Examples

A simple example is as below, it will use the TICA_TEST method to analyze the CVs cv1 and cv2 from metadynamics simulaiton (COLVAR file colvar.0.data)
It will analyze the data using 20 different lag times (from 0 to 100*20=2000).

\verbatim
cv1: READ FILE=colvar.0.data VALUES=cv1 IGNORE_TIME
cv2: READ FILE=colvar.0.data VALUES=cv1 IGNORE_TIME
rbias: READ FILE=colvar.0.data VALUES=metad.rbias IGNORE_TIME
rw: REWEIGHT_METAD TEMP=330

TICA_TEST ...
 ARG=cv1,cv2
 LAG_TIME=400
 TAU_NUMBER=100
 STEP_SIZE=0.2
 LOGWEIGHTS=rw
... TICA_TEST
\endverbatim

Using plumed \ref driver to perform the analysis, then you will get the eigenvalues file (default name eigenvalue.data) and largest eigenvectors (default name eigenvector0.data) at different lag times.

After using TICA_TEST method, you will also get the correlation file (default name correlation.data).
If you have several parallel trajectories (such as using multiple walkers) to be analyze,
this correlation file can be read as the restart file to analyze the next trajectory:

\verbatim
cv1: READ FILE=colvar.1.data VALUES=cv1 IGNORE_TIME
cv2: READ FILE=colvar.1.data VALUES=cv2 IGNORE_TIME
rbias: READ FILE=colvar.1.data VALUES=metad.rbias IGNORE_TIME
rw: REWEIGHT_METAD TEMP=330

TICA_TEST ...
 ARG=cv1,cv2
 LAG_TIME=400
 TAU_NUMBER=100
 STEP_SIZE=0.2
 LOGWEIGHTS=rw
 READ_CORR_FILE=correlation.data
... TICA_TEST
\endverbatim 

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class TICA_TEST : public vesselbase::ActionWithAveraging {
private:
	// If we are reusing data are we ignoring the reweighting in that data
	bool ignore_reweight;
	bool is_corr_info;
	bool is_rescale;
	bool is_read_corr;
	bool is_debug;
	bool is_int_time;
	bool use_int_calc;
	bool use_mark;

	// The piece of data we are inserting
	unsigned idata;
	unsigned narg;
	unsigned npoints;
	unsigned tot_steps;
	//~ unsigned rwsize;
	unsigned eff_size;
	unsigned ncomp;
	//~ unsigned nsize;
	unsigned dsize;
	
	unsigned mark_size;
	unsigned mark_now;
	unsigned mark_old;
	
	double lag_time;
	double delta_tau;
	double wnorm;
	double mnorm;
	double dt;
	double rw_time;
	double rw_temp;
	double mark_min;
	double mark_max;
	double mark_avg;
	double mark_sd;

	std::string eigval_file;
	std::string eigvec_file;
	std::string corr_file;
	std::string corr_file2;
	std::string corr_input;
	std::string corr_info_file;
	std::string rescale_file;
	std::string debug_file;
	std::string mark_name;
	
	OFile oeigval;
	OFile ocorr;
	OFile ocorr2;
	OFile ocorrinfo;
	OFile orescale;
	OFile odebug;
	
	IFile icorr;
	
	// The weights of all the data points
	std::vector<double> logweights;
	// The weights of all the data points
	std::vector<double> weights;
	// Tempory vector to store values of arguments
	std::vector<double> current_args;
	std::vector<double> myaverages;
	std::vector<double> Wnew;
	std::vector<double> atau;
	std::vector<double> neff;
	std::vector<double> point_weights;
	
	// List of argument names 
	std::vector<std::string> argument_names;
	std::vector<std::string> row_args;
	std::vector<std::string> col_args;
	
	// The data we are going to analyze
	std::vector<std::vector<double> > data;
	std::vector<std::vector<double> > avgdata;
	std::vector<std::vector<double> > Rnew;
	
	std::vector<double> mark_norm;
	std::vector<std::vector<double> > mark_weights;
	std::vector<std::vector<double> > mark_neff;
	std::vector<std::vector<std::vector<double> > > mark_data;
	
	//~ std::vector<double> mywts;
	std::vector<bool> nandata;
	std::vector<std::vector<double>::size_type> effid;
	std::vector<Matrix<double> > Cread;
	
	Matrix<double> Czero;
	double get_neff(const std::vector<double>& _weights,std::vector<double>& _neff);
	double get_correlation_integer(const std::vector<double>& X,
		const std::vector<double>& Y,unsigned _tau,double unorm);
	double get_correlation_float(const std::vector<double>& X,
		const std::vector<double>& Y,const std::vector<unsigned>& xid,
		const std::vector<unsigned>& yid,const std::vector<double>& twt,
		double fnorm);
public:
	static void registerKeywords( Keywords& keys );
	explicit TICA_TEST( const ActionOptions& );
	~TICA_TEST();
	void calculate(){}
	void apply(){}
	void performOperations( const bool& from_update );
	void performTask( const unsigned& , const unsigned& , MultiValue& ) const { plumed_error(); }
	
	void accumulate();
	void performAnalysis();
	void runFinalJobs();
	void runAnalysis();
	bool isPeriodic(){ plumed_error(); return false; }
	
	Matrix<double> build_correlation_integer(const std::vector<std::vector<double> >& _data,double& unorm,unsigned _tau=0);
	Matrix<double> build_correlation_float(const std::vector<std::vector<double> >& _data,
	const std::vector<double>& _weights,const std::vector<double>& _neff,double& fnorm,double _tau=0);
};

PLUMED_REGISTER_ACTION(TICA_TEST,"TICA_TEST")

void TICA_TEST::registerKeywords( Keywords& keys ){
	vesselbase::ActionWithAveraging::registerKeywords( keys );
	keys.remove("SERIAL"); keys.remove("LOWMEM"); 
	keys.remove("ARG");
	keys.add("compulsory","ARG","the CVs for the input to calculate TICA_TEST"); 
  	keys.add("compulsory","LAG_TIME","the total lag time to calculate TICA_TEST");
	keys.add("compulsory","TAU_NUMBER","100","how much points of lag time to output");
	keys.add("compulsory","STEP_SIZE","1.0","the simulation time step size");
	keys.add("compulsory","EIGENVECTOR_FILE","eigenvector","the file to output the result");
	keys.add("compulsory","EIGENVALUE_FILE","eigenvalue.data","the file to output the eigen value");
	keys.add("compulsory","CORRELATION_FILE","correlation.data","the file to output the correlation matrix");
	keys.add("compulsory","CORR_CAN_FILE","correlation2.data","the file to output the correlation canonical matrix");
	keys.addFlag("UNIFORM_WEIGHTS",false,"make all the weights equal to one");
	keys.addFlag("USE_MARK",false,"use the mark in a discontinuational trajectory");
	keys.add("optional","EIGEN_NUMBERS","how many eigenvectors to be output (from large to small)");
	keys.add("optional","READ_CORR_FILE","read the correlations from file");
	keys.add("optional","CORR_INFO_FILE","the file that output the information of CVs correlation");
	keys.add("optional","RESCALE_FILE","the file that output rescaled trajectory");
	keys.add("optional","DEBUG_FILE","the file that debug information");
}

TICA_TEST::TICA_TEST(const ActionOptions&ao):
Action(ao),
ActionWithAveraging(ao),
ignore_reweight(false),is_corr_info(false),
is_rescale(false),is_read_corr(false),
is_debug(false),is_int_time(false),use_int_calc(false),
idata(0),narg(0),npoints(0),tot_steps(0),dsize(0),mark_old(0),
delta_tau(0.0),wnorm(0.0),mnorm(0.0),
mark_min(1e308),mark_max(0),mark_avg(0),mark_sd(0)
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
	
	parseFlag("USE_MARK",use_mark);
	
	if(use_mark)
	{
		--narg;
	}
	
	current_args.resize(narg);
	myaverages.resize(narg);
	argument_names.resize(narg);
	data.resize(narg);
	
	for(unsigned i=0;i!=narg;++i)
	{
		argument_names[i]=getPntrToArgument(i)->getName();
		std::string id;
		Tools::convert(i,id);
		std::string colid="COL"+id+"_"+argument_names[i];
		std::string rowid="ROW"+id+"_"+argument_names[i];
		col_args.push_back(colid);
		row_args.push_back(rowid);
	}
	if(use_mark)
		mark_name=getPntrToArgument(narg)->getName();
	
	//~ parse("DELTA_TAU",delta_tau);
	parse("LAG_TIME",lag_time);
	if(lag_time<=0)
		error("LAG_TIME could be larger than 0");
	parse("TAU_NUMBER",npoints);
	if(npoints==0)
		error("TAU_NUMBER could not be 0");
	delta_tau=lag_time/npoints;
	++npoints;
	point_weights.assign(npoints,0);
	
	parse("STEP_SIZE",dt);
	if(dt<=0)
		error("STEP_SIZE could be larger than 0");
		
	double tsize=delta_tau/dt;
	dsize=floor(tsize+0.5);
	if(fabs(tsize-dsize)<1e-6)
	{
		is_int_time=true;
		//~ nsize=floor(lag_time/dt+0.5);
	}

	ncomp=narg;
	parse("EIGEN_NUMBERS",ncomp);
	plumed_massert(ncomp>0,"the EIGEN_NUMBER must be larger than 0!");
	plumed_massert(ncomp<=narg,"the EIGEN_NUMBER cannot be larger than the number of CVs!");
	
	parse("CORRELATION_FILE",corr_file);
	parse("CORR_CAN_FILE",corr_file2);
	parse("EIGENVECTOR_FILE",eigvec_file);
	parse("EIGENVALUE_FILE",eigval_file);

	parseFlag("UNIFORM_WEIGHTS",ignore_reweight);
	

	if(ignore_reweight&&is_int_time)
		use_int_calc=true;

	parse("CORR_INFO_FILE",corr_info_file);
	if(corr_info_file.size()>0)
	{
		is_corr_info=true;
		ocorrinfo.link(*this);
		ocorrinfo.open(rescale_file.c_str());
		ocorrinfo.fmtField(" %e");
		ocorrinfo.addConstantField("Lag_Time");	
	}
	
	parse("RESCALE_FILE",rescale_file);
	if(rescale_file.size()>0)
	{
		is_rescale=true;
		orescale.link(*this);
		orescale.open(rescale_file.c_str());
		orescale.fmtField(" %e");
	}
	
	parse("DEBUG_FILE",debug_file);
	if(debug_file.size()>0)
	{
		is_debug=true;
		odebug.link(*this);
		odebug.open(debug_file.c_str());
		odebug.fmtField(" %e");
		odebug.addConstantField("Lag_Time");
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
			
		for(unsigned ip=0;ip!=npoints;++ip)
		{
			double lag_time=0;
			icorr.scanField("LAG_TIME",lag_time);
			double rnorm=0;
			if(icorr.FieldExist("NORMALIZATION"))
				icorr.scanField("NORMALIZATION",rnorm);
			point_weights[ip]=rnorm;
			plumed_massert(fabs(ip*delta_tau-lag_time)<1.0e-6,"the lag time read from the correlation file mismatch!");
			atau.push_back(unsigned(lag_time));
			Matrix<double> cc(narg,narg);
			for(unsigned i=0;i!=narg;++i)
			{
				std::string rowarg;
				icorr.scanField("MATRIX",rowarg);
				plumed_massert(rowarg==row_args[i],"the label of rows read from the correlation file mismatch! (\""+rowarg+"\" vs \""+row_args[i]+"\")");
				for(unsigned j=0;j!=narg;++j)
				{
					double mm;
					if(icorr.FieldExist(col_args[j]))
						icorr.scanField(col_args[j],mm);
					else
						error("could not found the label \""+col_args[j]+"\"");
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
	ocorr.fmtField(" %e");
	ocorr.addConstantField("DIMENSION");
	ocorr.printField("DIMENSION",int(narg));
	ocorr.addConstantField("NPOINTS");
	ocorr.printField("NPOINTS",int(npoints));
	ocorr.addConstantField("LAG_TIME");
	ocorr.addConstantField("NORMALIZATION");
	
	ocorr2.link(*this);
	ocorr2.open(corr_file2.c_str());
	ocorr2.fmtField(" %e");
	ocorr2.addConstantField("DIMENSION");
	ocorr2.printField("DIMENSION",int(narg));
	ocorr2.addConstantField("NPOINTS");
	ocorr2.printField("NPOINTS",int(npoints));
	ocorr2.addConstantField("LAG_TIME");
	ocorr2.addConstantField("NORMALIZATION");
	
	oeigval.link(*this);
	oeigval.open(eigval_file.c_str());

	checkRead();
	
	log.printf("  with %d arguments.\n",narg);
	log.printf("  with toal lag time: %f.\n",lag_time);
	log.printf("  with delta lag time: %f.\n",delta_tau);
	log.printf("  with number of time points: %d.\n",npoints);
	log.printf("  with time step size: %f.\n",dt);
	if(use_mark)
		log.printf("  with mark name: %s \n",mark_name.c_str());
	if(ignore_reweight)
		log.printf("  use uniform weights to calculate\n");
	if(use_int_calc)
		log.printf("  use integal algorithm to calculate all correlation\n");
	log.printf("  with eigen values output file: %s\n",eigval_file.c_str());
	log.printf("  with %d eigen vector output files started with name: %s\n",ncomp,eigvec_file.c_str());
	if(is_rescale)
		log.printf("  with rescaled output file: %s\n",rescale_file.c_str());
	if(is_debug)
		log.printf("  with debug file: %s\n",debug_file.c_str());
	if(is_read_corr)
	{
		log.printf("  read correlation from file: %s.\n",corr_input.c_str());
		log.printf("  with lag times:\n");
		for(unsigned i=0;i!=atau.size();++i)
			log.printf("    %d\t%f\t%f\n",int(i),atau[i],point_weights[i]);
	}
	log<<"Bibliography "<<plumed.cite("McCarty and Parrinello, J. Chem. Phys. 147, 204109 (2017)");
}

TICA_TEST::~TICA_TEST()
{
	ocorr.close();
	ocorr2.close();
	//~ for(unsigned i=0;i!=oeigvecs.size();++i)
		//~ oeigvecs[i].close();
	oeigval.close();
	if(is_corr_info)
		ocorrinfo.close();
	if(is_rescale)
		orescale.close();
	if(is_debug)
		odebug.close();
}

void TICA_TEST::performAnalysis()
{
	tot_steps=logweights.size();
	if(ignore_reweight)
		rw_time=idata*dt;
		
	for(unsigned i=0;i!=narg;++i)
		myaverages[i]/=wnorm;

	if(use_mark)
	{
		if(weights.size()>0)
		{
			mark_data.push_back(data);
			for(unsigned i=0;i<narg;++i)
			{
				data[i].clear();
				//~ data[i].swap(std::vector<double>());
			}
			mark_weights.push_back(weights);
			weights.clear();
			//~ weights.swap(std::vector<double>());
			mark_norm.push_back(mnorm);
			mark_avg+=mnorm;
			mark_sd+=mnorm*mnorm;
			if(mnorm<mark_min)
				mark_min=mnorm;
			if(mnorm>mark_max)
				mark_max=mnorm;
			mnorm=0;
		}
		plumed_massert(mark_data.size()==mark_weights.size(),"the number of the marked sub-trajectories is different from the number of weights");
		mark_size=mark_data.size();
		mark_avg/=mark_size;
		mark_sd=sqrt(mark_sd/mark_size-mark_avg*mark_avg);
		for(unsigned k=0;k!=mark_data.size();++k)
		{
			for(unsigned i=0;i!=narg;++i)
			{
				for(unsigned j=0;j!=mark_data[k][i].size();++j)
					mark_data[k][i][j]-=myaverages[i];
			}
		}
	}
	else
	{
		avgdata.resize(narg,std::vector<double>(tot_steps));
		for(unsigned i=0;i!=narg;++i)
		{
			for(unsigned j=0;j!=tot_steps;++j)
				avgdata[i][j]=data[i][j]-myaverages[i];
		}
	}
	
	log.printf("Fininished reading.\n");
	log.printf("  with reading steps: %d.\n",idata);
	log.printf("  with total steps: %d.\n",tot_steps);
	if(use_mark)
	{
		log.printf("  with number of sub-trajectories: %d.\n",int(mark_data.size()));
		log.printf("  with reweight times from %f to %f.\n",mark_min,mark_max);
		log.printf("  with average reweight time %f and standard deviation %f.\n",mark_avg,mark_sd);
	}


	if(!use_int_calc)
	{
		if(use_mark)
		{
			rw_time=0;
			mark_neff.resize(mark_weights.size());
			for(unsigned i=0;i!=mark_weights.size();++i)
			{
				double rwt=get_neff(mark_weights[i],mark_neff[i]);
				rw_time+=rwt;
			}
			if(is_rescale)
			{
				int count=0;
				for(unsigned i=0;i!=mark_size;++i)
				{
					for(unsigned j=0;i!=mark_weights[i].size();++j)
					{
						orescale.printField("INDEX",count);
						orescale.printField("ID-I",int(i));
						orescale.printField("ID-J",int(j));
						orescale.printField("RESCALED_TIME",mark_neff[i][j]);
						for(unsigned k=0;k!=narg;++k)
							orescale.printField(argument_names[k],mark_data[i][k][j]);
						orescale.printField("WEIGHT",mark_weights[i][j]);
						orescale.printField();
						++count;
					}
				}
				orescale.flush();
			}
		}
		else
		{
			rw_time=get_neff(weights,neff);
			if(is_rescale)
			{
				for(unsigned i=0;i!=tot_steps;++i)
				{
					orescale.printField("INDEX",int(i));
					orescale.printField("RESCALED_TIME",neff[i]);
					for(unsigned j=0;j!=narg;++j)
						orescale.printField(argument_names[j],data[j][i]);
					orescale.printField("WEIGHT",weights[i]);
					orescale.printField();
				}
				orescale.flush();
			}
		}
	}
	log.printf("  with reweight time steps: %f.\n",rw_time);

	double xnorm=0;
	if(use_mark)
	{
		Czero.resize(narg,narg);
		for(unsigned i=0;i!=narg;++i)
		{
			for(unsigned j=0;j!=narg;++j)
				Czero[i][j]=0;
		}
		for(unsigned i=0;i!=mark_size;++i)
		{
			Matrix<double> Cm0;
			mnorm=0;
			if(ignore_reweight)
				Cm0=build_correlation_integer(mark_data[i],mnorm,0);
			else
				Cm0=build_correlation_float(mark_data[i],mark_weights[i],mark_neff[i],mnorm,0);
			for(unsigned i=0;i!=narg;++i)
			{
				for(unsigned j=0;j!=narg;++j)
					Czero[i][j]+=Cm0[i][j]*mnorm;
			}
			xnorm+=mnorm;
		}
		for(unsigned i=0;i!=narg;++i)
		{
			for(unsigned j=0;j!=narg;++j)
				Czero[i][j]/=xnorm;
		}
	}
	else
	{
		if(ignore_reweight)
			Czero=build_correlation_integer(avgdata,xnorm,0);
		else
			Czero=build_correlation_float(avgdata,weights,neff,xnorm,0);
	}

	if(is_read_corr)
	{
		double cr0=point_weights[0]/(point_weights[0]+xnorm);
		double cr1=xnorm/(point_weights[0]+xnorm);
		for(unsigned i=0;i!=narg;++i)
		{
			for(unsigned j=0;j!=narg;++j)
				Czero[i][j]=cr0*Cread[0][i][j]+cr1*Czero[i][j];
		}
	}
	
	log.printf("  with the covariance matrix:.\n");
	for(unsigned i=0;i!=narg;++i)
	{
		log.printf("   ");
		for(unsigned j=0;j!=narg;++j)
			log.printf(" %f",Czero[i][j]);
		log.printf("\n");
	}
	
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

	std::vector<std::vector<double> >& U=eigvec0;
	
	unsigned eff_size=0;
	for(unsigned i=0;i!=narg;++i)
		if(eigval0[i]>0)	++eff_size;
	log.printf("  with %d effective eigenvalues and eigenvectors.\n",int(eff_size));
	
	plumed_massert(eff_size>0,"All the eigenvalues are less than 0, TICA_TEST could not be performed!");
	
	if(ncomp>eff_size)
	{
		log.printf("  as the EIGEN_NUMBER (%d) is larger than the effective one, it will be reduced.\n",int(ncomp));
		ncomp=eff_size;
	}
	
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
	
	std::vector<std::vector<std::vector<double> > > eigvec_points;
	log.printf("Running TICA_TEST ...\n");
	//~ std::vector<unsigned> compid(ncomp,0);
	std::vector<std::vector<double> > pre_evec(ncomp);
	std::vector<bool> resign(ncomp,false);
	for(unsigned ip=0;ip!=npoints;++ip)
	{
		double tau=ip*delta_tau;
		
		double ntau=tau/dt;
		unsigned int_tau=floor(ntau+0.5);
		bool is_use_int_calc=false;
		if(use_int_calc||(ignore_reweight&&fabs(int_tau-ntau)<1e-6))
		{
			is_use_int_calc=true;
			log.printf("  Setp: %d\tLag Time: %f (%f * %d)\n",ip,tau,dt,int_tau);
		}
		else
			log.printf("  Setp: %d\tLag Time: %f\n",ip,tau);
		
		if(tau>rw_time)
		{
			log.printf("  the lag time is too large. The calculation will be stop.\n");
			break;
		}

		xnorm=0;
		Matrix<double> Clag(narg,narg);
		if(use_mark)
		{
			for(unsigned i=0;i!=narg;++i)
			{
				for(unsigned j=0;j!=narg;++j)
					Clag[i][j]=0;
			}
			unsigned mcount=0;
			for(unsigned i=0;i!=mark_size;++i)
			{
				if(mark_norm[i]>tau||(use_int_calc&&mark_weights.size()>int_tau))
				{
					Matrix<double> Cml;
					mnorm=0;
					if(is_use_int_calc)
						Cml=build_correlation_integer(mark_data[i],mnorm,int_tau);
					else
						Cml=build_correlation_float(mark_data[i],mark_weights[i],mark_neff[i],mnorm,tau);

					for(unsigned i=0;i!=narg;++i)
					{
						for(unsigned j=0;j!=narg;++j)
							Clag[i][j]+=Cml[i][j]*mnorm;
					}
					xnorm+=mnorm;
					++mcount;
				}
			}
			if(xnorm>0)
			{
				for(unsigned i=0;i!=narg;++i)
				{
					for(unsigned j=0;j!=narg;++j)
						Clag[i][j]/=xnorm;
				}
				log.printf("  with %d effective sub-trajectories\n",mcount);
			}
			else
			{
				log.printf("  the lag time is too large. The calculation will be stop.\n");
				break;
			}
		}
		else
		{
			if(is_use_int_calc)
				Clag=build_correlation_integer(avgdata,xnorm,int_tau);
			else
				Clag=build_correlation_float(avgdata,weights,neff,xnorm,tau);
		}
		if(is_read_corr)
		{
			double ncr0=point_weights[ip]/(point_weights[ip]+xnorm);
			double ncr1=xnorm/(point_weights[ip]+xnorm);
			for(unsigned i=0;i!=narg;++i)
			{
				for(unsigned j=0;j!=narg;++j)
					Clag[i][j]=ncr0*Cread[ip][i][j]+ncr1*Clag[i][j];
			}
			xnorm+=point_weights[ip];
		}
		
		ocorr.printField("LAG_TIME",tau);
		ocorr2.printField("LAG_TIME",tau);

		for(unsigned i=0;i!=narg;++i)
		{
			ocorr.printField("NORMALIZATION",xnorm);
			ocorr.printField("MATRIX",row_args[i]);
			for(unsigned j=0;j!=narg;++j)
				ocorr.printField(col_args[j],Clag[i][j]);
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
		
		//~ C_can = 0.5*(C_can+C_can.trans());
		Matrix<double> C_can_T;
		transpose(C_can,C_can_T);
		for(unsigned i=0;i!=narg;++i)
		{
			for(unsigned j=0;j!=narg;++j)
				C_can[i][j]=(C_can[i][j]+C_can_T[i][j])/2;
		}
		
		for(unsigned i=0;i!=narg;++i)
		{
			ocorr2.printField("NORMALIZATION",xnorm);
			ocorr2.printField("MATRIX",row_args[i]);
			for(unsigned j=0;j!=narg;++j)
				ocorr2.printField(col_args[j],C_can[i][j]);
			ocorr2.printField();
		}
		ocorr2.flush();
		
		if(is_debug)
		{
			odebug.printField("LAG_TIME",tau);
			for(unsigned i=0;i!=narg;++i)
			{
				odebug.printField("ROW",int(i));
				for(unsigned j=0;j!=narg;++j)
					odebug.printField(col_args[j],C_can[i][j]);
				odebug.printField();
			}
			odebug.flush();
		}
		
		std::vector<double> nwt;
		Matrix<double> nvt;
		diagMat(C_can,nwt,nvt);
		
		std::vector<double> eigval2;
		std::vector<std::vector<double> > eigvec2;
		if(ip==0)
		{
			for(unsigned i=0;i!=nwt.size();++i)
			{
				std::vector<double> rowvec;
				for(unsigned j=0;j!=narg;++j)
					rowvec.push_back(nvt[i][j]);
				eigval2.push_back(nwt[i]);
				eigvec2.push_back(rowvec);
			}
		}
		else
		{
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
			for(std::multimap<double,std::vector<double> >::reverse_iterator miter=eigs2.rbegin();miter!=eigs2.rend();++miter)
			{
				eigval2.push_back(miter->first);
				eigvec2.push_back(miter->second);
			}
		}
		
		log.printf("  Eigenvalues:");
		for(unsigned i=0;i!=narg;++i)
			log.printf(" %f",eigval2[i]);
		log.printf("\n\n");

		oeigval.fmtField(" %f");
		oeigval.printField("LAG_TIME",tau);
		for(unsigned i=0;i!=narg;++i)
		{
			std::string id;
			Tools::convert(i,id);
			std::string eigid="EIGVAL"+id;
			oeigval.printField(eigid,eigval2[i]);
		}
		oeigval.printField();
		oeigval.flush();
		
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
		
		std::vector<std::vector<double> > tica_res;
		for(unsigned k=0;k!=ncomp;++k)
		{
			//~ std::vector<double> evec=mvt.get_column(k);
			std::vector<double> evec;
			
			double mnorm=0;
			for(unsigned j=0;j!=narg;++j)
			{
				double val=mvt[j][k];
				evec.push_back(val);
				mnorm+=val*val;
			}
			mnorm=sqrt(mnorm);
			
			double vmax=0;
			unsigned cpid=0;
			for(unsigned j=0;j!=narg;++j)
			{
				if(fabs(mvt[j][k])>vmax)
				{
					vmax=fabs(mvt[j][k]);
					cpid=j;
				}
			}
			
			if(ip>0)
			{
				if(evec[cpid]*pre_evec[k][cpid]<0)
					mnorm*=-1;
				if(evec[cpid]<0)
					resign[k]=true;
				else
					resign[k]=false;
			}

			for(unsigned i=0;i!=narg;++i)
				evec[i]/=mnorm;

			tica_res.push_back(evec);
			
			pre_evec[k]=evec;
		}
		eigvec_points.push_back(tica_res);

		for(unsigned i=0;i!=eff_size;++i)
			if(nwt[i]<0) nwt[i]=0;
	}
	
	for(unsigned i=0;i!=ncomp;++i)
	{
		OFile oeigvec;
		oeigvec.link(*this);
		std::string id;
		Tools::convert(i,id);
		std::string name=eigvec_file+id+".data";
		oeigvec.open(name.c_str());
		oeigvec.addConstantField("COMPONENT");
		oeigvec.printField("COMPONENT",int(i));
		oeigvec.fmtField(" %e");
		for(unsigned j=0;j!=npoints;++j)
		{
			double tau=j*delta_tau;
			oeigvec.printField("LAG_TIME",tau);
			for(unsigned k=0;k!=narg;++k)
			{
				double vvec=eigvec_points[j][i][k];
				if(resign[i])
					vvec*=-1;
				std::string id;
				Tools::convert(k,id);
				std::string eigid="EIGVEC"+id+"_"+argument_names[k];
				oeigvec.printField(eigid,vvec);
			}
			oeigvec.printField();
		}
		oeigvec.flush();
		oeigvec.close();
	}
}

double TICA_TEST::get_neff(const std::vector<double>& _weights,std::vector<double>& _neff)
{
	std::vector<double> weights2(_weights);
	weights2.insert(weights2.begin(),0);
	_neff.resize(weights2.size());
	std::partial_sum(weights2.begin(),weights2.end(),_neff.begin());
	double _time=_neff.back();
	_neff.pop_back();
	return _time;
}

void TICA_TEST::accumulate(){
	// Get the arguments ready to transfer to reference configuration
	for(unsigned i=0;i<narg;++i)
		current_args[i]=getArgument(i);

	if(use_mark)
	{
		mark_now=getArgument(narg);
		if(mark_now!=mark_old&&mark_old!=0)
		{
			mark_data.push_back(data);
			for(unsigned i=0;i<narg;++i)
			{
				data[i].clear();
				//~ data[i].swap(std::vector<double>());
			}
			mark_weights.push_back(weights);
			weights.clear();
			//~ weights.swap(std::vector<double>());
			mark_norm.push_back(mnorm);
			
			mark_avg+=mnorm;
			mark_sd+=mnorm*mnorm;
			if(mnorm<mark_min)
				mark_min=mnorm;
			if(mnorm>mark_max)
				mark_max=mnorm;
			
			mnorm=0;
		}
		mark_old=mark_now;
	}

	if(!use_mark||mark_now!=0)
	{
		double cdw=cweight*dt;
		double ldw=lweight+std::log(dt);
		for(unsigned i=0;i<narg;++i)
		{
			data[i].push_back(current_args[i]);
			if(ignore_reweight)
				myaverages[i]+=current_args[i]*dt;
			else
				myaverages[i]+=current_args[i]*cdw;
		}
		// Get the arguments and store them in a vector of vectors
		logweights.push_back(ldw);

		double ww=cdw;
		if(ignore_reweight)
			ww=dt;

		weights.push_back(ww);
		wnorm+=ww;

		if(use_mark)
			mnorm+=ww;
	}
	// Increment data counter
	++idata;
}

void TICA_TEST::performOperations( const bool& from_update ){
  accumulate();
}

void TICA_TEST::runAnalysis(){
	performAnalysis(); idata=0;
}

void TICA_TEST::runFinalJobs() {
  runAnalysis(); 
}

Matrix<double> TICA_TEST::build_correlation_integer(const std::vector<std::vector<double> >& _data,double& unorm,unsigned _tau)
{
	if(!ignore_reweight)
		error("The integer algorithm can be only used with the flag IGNORE_REWEIGHT");
	const std::vector<std::vector<double> >& X=_data;
	std::vector<std::vector<double> > Y(_data);
	if(_tau!=0)
	{
		for(unsigned i=0;i!=narg;++i)
			Y[i].erase(Y[i].begin(),Y[i].begin()+_tau);
	}
	unorm=Y[0].size()*dt;

	Matrix<double> corrmat(narg,narg);
	for(unsigned i=0;i!=narg;++i)
	{
		for(unsigned j=0;j!=narg;++j)
			corrmat[i][j]=get_correlation_integer(X[i],Y[j],_tau,unorm);
	}
	return corrmat;
}

Matrix<double> TICA_TEST::build_correlation_float(const std::vector<std::vector<double> >& _data,
	const std::vector<double>& _weights,const std::vector<double>& _neff,
	double& fnorm,double _tau)
{
	std::vector<unsigned> xid;
	std::vector<unsigned> yid;
	std::vector<double> twt;
	fnorm=0;
	if(_tau==0)
	{
		for(unsigned i=0;i!=_weights.size();++i)
		{
			xid.push_back(i);
			yid.push_back(i);
			fnorm+=_weights[i];
		}
		twt=_weights;
	}
	else
	{
		// the beginnig time;
		double tbeg=0;
		// the end time of lag time (begin time plus lag time);
		double tend=_tau;
		
		// the index of the beginning time;
		unsigned ibeg=0;
		// the index of the end time;
		//~ unsigned iend=_neff.size();
		//~ for(unsigned i=0;i!=_neff.size()-1;++i)
		//~ {
			//~ if(_tau>=_neff[i]&&_tau<_neff[i+1])
			//~ {
				//~ iend=i;
				//~ break;
			//~ }
		//~ }
		unsigned iend=_neff.size();
		for(unsigned i=0;i!=_neff.size();++i)
		{
			if(_tau>=_neff[i]&&_tau<(_neff[i]+_weights[i]))
			{
				iend=i;
				break;
			}
		}
		
		if(iend==_neff.size())
			error("Cannot find the index iend: the lag time is too large!");
		// the different time to the next index
		double dbeg=_neff[ibeg]+_weights[ibeg]-tbeg;
		double dend=_neff[iend]+_weights[iend]-tend;

		while(iend<_neff.size())
		{
			xid.push_back(ibeg);
			yid.push_back(iend);
			
			if(dbeg<dend)
			{
				twt.push_back(dbeg);
				fnorm+=dbeg;
				tend+=dbeg;
				dend-=dbeg;

				++ibeg;
				tbeg=_neff[ibeg];
				dbeg=_weights[ibeg];
			}
			else if(dbeg>dend)
			{
				twt.push_back(dend);
				fnorm+=dend;
				tbeg+=dend;
				dbeg-=dend;

				++iend;
				if(iend==_neff.size())
					break;
				tend=_neff[iend];
				dend=_weights[iend];
			}
			else
			{
				twt.push_back(dbeg);
				fnorm+=dbeg;
				++ibeg;
				tbeg=_neff[ibeg];
				dbeg=_weights[ibeg];
				
				++iend;
				if(iend==_neff.size())
					break;
				tend=_neff[iend];
				dend=_weights[iend];
			}
		}
	}
	if(is_corr_info)
	{	
		ocorrinfo.printField("Lag_Time",_tau);		
		for(unsigned i=0;i!=twt.size();++i)
		{
			ocorrinfo.printField("Index",int(i));
			ocorrinfo.printField("XID",int(xid[i]));
			ocorrinfo.printField("YID",int(yid[i]));
			ocorrinfo.printField("Weight",twt[i]);
			ocorrinfo.printField();
		}
		ocorrinfo.flush();
	}
	
	Matrix<double> corrmat(narg,narg);
	for(unsigned i=0;i!=narg;++i)
	{
		for(unsigned j=0;j!=narg;++j)
			corrmat[i][j]=get_correlation_float(_data[i],_data[j],xid,yid,twt,fnorm);
	}
	return corrmat;
}

double TICA_TEST::get_correlation_integer(const std::vector<double>& X,const std::vector<double>& Y,unsigned _tau,double unorm)
{
	double corrmean=0;
	
	for(unsigned id=0;id!=Y.size();++id)
	{
		double value=X[id]*Y[id];
		corrmean+=value*dt;
	}
	
	corrmean/=unorm;

	return corrmean;
}

double TICA_TEST::get_correlation_float(const std::vector<double>& X,const std::vector<double>& Y,const std::vector<unsigned>& xid,const std::vector<unsigned>& yid,const std::vector<double>& twt,double fnorm)
{
	double corrmean=0;
	
	plumed_massert(xid.size()==yid.size(),"the input indices mismatch");
	plumed_massert(xid.size()==twt.size(),"the indices and weigths mismatch");
	
	for(unsigned i=0;i!=xid.size();++i)
	{
		unsigned ix=xid[i];
		unsigned iy=yid[i];
		double ww=twt[i];
		double value=X[ix]*Y[iy];
		corrmean+=value*ww;
	}

	corrmean/=fnorm;
	return corrmean;
}

}
}
