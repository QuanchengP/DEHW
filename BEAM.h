
#include "VEM.h"

class BEAM: public VEM{
public:
	//
	Matrix<double,3,1> leng;//geometric length
	double lengFact;
	double angl;//geometric angle
	Matrix<long,3,1> diviNumb;//number of elements along different directions
	//
	double loadInte;//load intensity
	BEAM();
	bool MESH();//generate mesh
	bool SOLVE();//solve for displacement
};

BEAM::BEAM(){
	//
	mkdir("Beam", 0755);
	outpDire = "Beam/";
	//
	leng << 1.0, 0.12, 0.06;
	lengFact = 1.0 / 3.0;
	angl = 45.0 * PI / 180.0;
	diviNumb << 100, 12, 6;
	loadInte = -8000.0;
	MESH();
}

bool BEAM::MESH(){
	//points
	for(long ti = 0; ti <= diviNumb(0); ti ++){
		double angl_ti = angl * (double)ti / diviNumb(0);
		double widt_ti = leng(2) * (1.0- (double)ti / diviNumb(0) * lengFact);
		double heig_ti = leng(1) * (1.0- (double)ti / diviNumb(0) * lengFact);
		for(long tj = 0; tj <= diviNumb(1); tj ++){
			double yCoo = -heig_ti/2.0+heig_ti/diviNumb(1)*tj;
			for(long tk = 0; tk <= diviNumb(2); tk ++){
				double zCoo = -widt_ti/2.0+widt_ti/diviNumb(2)*tk;
				COORDINATE tempCoor;
				tempCoor <<
					leng(0)/diviNumb(0)*ti,
					cos(angl_ti)*yCoo - sin(angl_ti)*zCoo,
					sin(angl_ti)*yCoo + cos(angl_ti)*zCoo;
				TRY_ADD_POINT(tempCoor);
			}
		}
	}
	//volumes
	long layeNumb = (diviNumb(1)+1)*(diviNumb(2)+1);
	for(long ti = 0; ti < diviNumb(0); ti ++){
        if(ti % 10 == 0){
            cout<<"Element layer "<<ti<<"/"<<diviNumb(0)<<endl;
        }
		long tempI = ti*layeNumb;
		for(long tj = 0; tj < diviNumb(1); tj ++){
			long tempJ = tempI+tj*(diviNumb(2)+1);
			for(long tk = 0; tk < diviNumb(2); tk ++){
				Matrix<long,8,1> tempCorn;
				tempCorn <<
					tempJ+tk,
					tempJ+tk+(diviNumb(2)+1),
					tempJ+tk+1+(diviNumb(2)+1),
					tempJ+tk+1,
					tempJ+tk+layeNumb,
					tempJ+tk+(diviNumb(2)+1)+layeNumb,
					tempJ+tk+1+(diviNumb(2)+1)+layeNumb,
					tempJ+tk+1+layeNumb;
				ADD_BLOCK_FROM_POINT(tempCorn);
			}
		}
	}
	VOLUME_2_ELEMENT();
	OUTPUT_ELEMENT(0);
	return true;
}

bool BEAM::SOLVE(){
	freeTran.resize(nodeCoor.rows(),1);
	for(long ti = 0; ti < nodeCoor.rows(); ti ++){
		freeTran(ti) << 1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0;
	}
	//
	vector<long> consFree;
	for(long ti = 0; ti < nodeCoor.rows(); ti ++){
		if(nodeCoor(ti,0) <= 1.0E-10){
			consFree.push_back(3*ti+0);
			consFree.push_back(3*ti+1);
			consFree.push_back(3*ti+2);
		}
	}
	CONSTRAINT(consFree);
	//
	ELEMENT_STIFF();
	SparseMatrix<double> stif;
	stif.resize(freeNumb, freeNumb);
	stif.setFromTriplets(stifList.begin(), stifList.end());
    //
	OUTPUT_TIME("Load");
	MatrixXd forc = MatrixXd::Zero(freeNumb,1);
	for(long ti = 0; ti < nodeCoor.rows(); ti ++){
		if(abs(nodeCoor(ti,1)) < 1.0E-10 
			&& abs(nodeCoor(ti,2)) < 1.0E-10 && consDOF(3*ti + 2) != -1){
            forc(consDOF(3*ti + 2)) = loadInte*leng(0) / diviNumb(0);
		}
	}
	//solve
	OUTPUT_TIME("Solver");
	ConjugateGradient<SparseMatrix<double>, Lower|Upper> solv;
	solv.compute(stif);
	MatrixXd disp = solv.solve(forc);
	OUTPUT_DISPLACEMENT(disp, 0);
	OUTPUT_TIME("Done");
	return true;
}
