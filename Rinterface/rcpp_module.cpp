#include <matvec/geneticdist.h>   
#include <matvec/mim.h>
#include <matvec/pedigree.h>
#include <matvec/safe_vectors.h>

#include <Rcpp.h>
using namespace Rcpp;

class RPed {
public: 
	RPed(string pedFile){
		std::string colNames = "individual sire  dam";
		std::string colTypes = "CLASS      CLASS CLASS";
		ped.putColNames(colNames);
		ped.inputPed(pedFile);
		mim.pedPtr = &ped;
		mim.setupPopulation0(ped , G);
	}
    CharacterVector getNeighborsFor(string pivot, unsigned radius){
		std::vector<string> neighbors = mim.getNeighborsFor(pivot, radius);
        return wrap(neighbors); 
	}
    double getAij(string iID, string jID){
        matvec::PNode *iPNodem = ped[iID];
        matvec::PNode *jPNodem = ped[jID];
        int i = iPNodem->ind;
        int j = jPNodem->ind;
        return 2*ped.get_rij(i,j);
    }
	
	List getPedList(){
		std::vector<int> indVector, sireVector, damVector;
		std::vector<std::string> idVector;
		std::vector<double> fVector;
		unsigned n = ped.size();
		indVector.resize(n);
		sireVector.resize(n);
		damVector.resize(n);
		idVector.resize(n);
		fVector.resize(n);
		SafeSTLVector<matvec::PNode*>::iterator vecit;
		unsigned i = 0;
		for (vecit=ped.pedVector.begin();vecit!=ped.pedVector.end();vecit++){
			indVector[i]  = (*vecit)->ind; 
			sireVector[i] = (*vecit)->sire;
			damVector[i]  = (*vecit)->dam;
			idVector[i]   = (*vecit)->ind_str;
			fVector[i]    = (*vecit)->f;
			i++;
		}
		return DataFrame::create(Named("id")   = idVector,
								 Named("ind")  = indVector,
								 Named("sire") = sireVector,
								 Named("dam")  = damVector,
								 Named("F")    = fVector);
	}
	

	matvec::RPedigree ped;
	matvec::MIM mim;
	matvec::GeneticDist G;
};


RCPP_MODULE(pedigree){
	using namespace Rcpp ;
    class_<RPed>( "RPed")
	.constructor<std::string>()
	.method ( "getNeighborsFor", &RPed::getNeighborsFor, 
    "returns neighbors within radius")
    .method ("getAij", &RPed::getAij,
    "returns additive relationship between iID and jID")
	.method ("dataFrame", &RPed::getPedList,
	"returns pedigree as a data frame")		 
	;
}                     

