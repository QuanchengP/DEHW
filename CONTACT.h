
#ifndef _CONTACT_H
#define _CONTACT_H

#include "VEM.h"

/*************************************************************************************************/
/*******************************************DECLARATION*******************************************/
/*************************************************************************************************/

//Two polygons in 2d are intersecting, the minimum value of the intersecting area
//If the area of a triangle is less than miniArea, the triangle is degenerated to a line-segment
const double miniArea = 1.0E-14;

//quadrature over non-mortar side, one quadrature point = one CONTACT_ELEMENT
//project mortar triangle onto non-mortar triangle, a polygon is obtained
//divide the polygon into sub triangles by its centroid
//2*2 Guass quadrature over one sub triangle
//one sub triangle = four quadrature points = four CONTACT_ELEMENTs
class CONTACT_ELEMENT{
public:
	Matrix<long,3,2> node;//three nodes: one element face
	Matrix<double,3,2> shap;//three shape/basis functions
	Matrix<double,3,1> dual;//three dual basis functions on non-mortar side
	Matrix<double,3,2> contPoin;
	Vector3d normVect;//normal vector of non-mortar side
	double initGap;//= normVect * (contPoin 1 - contPoin 0)
	double quadWeig;//quadrature weight
};

class CONTACT{
public:
	CONTACT();
	//
	Matrix<VEM,2,1> vem;//0 - non-mortar side, 1 - mortar side
	//n*3, n segments/triangles/element faces, 3 nodes per segment
	Matrix<long,Dynamic,Dynamic> nonmSegm;
	Matrix<long,Dynamic,Dynamic> mortSegm;
	bool OUTPUT_CONTACT();//output nonmSegm and mortSegm
	//
	//global contact search, bucket search
	//2d local coordinate, 0 - infimum, 1 - supremum, 2 - increment
	Matrix<double,2,3> buckLoca;
	Matrix<long,2,1> buckNumb;//number of buckets along 2 directions
	Matrix<list<long>,Dynamic,Dynamic> buck;//bucket
	//sort non-mortar segments into different buckets
	bool BUCKET_SORT(const Matrix<double,Dynamic,Dynamic> &locaCoor, 
		const MatrixXd featSize, long diviNumb_0, long diviNumb_1);
	//
	//Guass quadrature over a triangle, shape function = triangle area coordinate
	//c1~2: two shape functions N1~N2 of quadrature point, N3 = 1 - N1 - N2
	//c3: quadrature weight
	Matrix<double,4,3> quadTabl;
	Matrix<CONTACT_ELEMENT,Dynamic,Dynamic> contElem;
	//local contact search
	//sort mortar segment into a bucket, only test the non-mortar segments already in this bucket
	bool CONTACT_SEARCH(const Matrix<double,Dynamic,Dynamic> &locaCoor);
	//intersection of two triangle, project onto 2d non-mortar plane, recover into 3d
	bool TRIANGLE_INTERSECT(long tempNonm, long tempMort, vector<CONTACT_ELEMENT> &tempCont);
	//sub procedure of TRIANGLE_INTERSECT
	bool TI_SUB(const Matrix<Vector3d,3,1> &nonmTria, const Matrix<Vector3d,3,1> &mortTria, 
		Vector3d &normVect, double &origArea, 
		Matrix<Vector2d,3,1> &nonmProj, Matrix<Vector2d,3,1> &mortProj, 
		vector<Vector2d> &tempQuad, vector<double> &tempWeig);
	//Guass quadrature over a triangle
	bool TRIANGLE_QUADRATURE(Vector2d tempPoin_0, Vector2d tempPoin_1, Vector2d tempPoin_2, 
		vector<Vector2d> &tempQuad, vector<double> &tempWeig, vector<double> &tempArea);
	//one line-segment is cross another line
	bool IS_CROSS_2D(Vector2d tempPoin_0, Vector2d tempPoin_1, 
		Vector2d tempPoin_2, Vector2d tempPoin_3);
	//sort tempPoin_0 and tempPoin_1 by their tempIden-th coordinate component
	bool SORT_BY_2D(Vector2d &tempPoin_0, Vector2d &tempPoin_1, long tempIden);
	//the intersect of two line-segments: none, one point, one line-segment (collinear)
	bool LINE_INTERSECT_2D(Vector2d tempPoin_0, Vector2d tempPoin_1, 
		Vector2d tempPoin_2, Vector2d tempPoin_3, vector<Vector2d> &tempInte);
	//is the point in the triangle
	bool IN_TRIANGLE_2D(Vector2d tempPoin, Matrix<Vector2d,3,1> tempTria);
	//the area of the triangle
	double TRIANGLE_AREA_2D(Vector2d tempPoin_0, Vector2d tempPoin_1, Vector2d tempPoin_2);
	//output contElem
	bool OUTPUT_CONTACT_ELEMENT();
	//
	long nonmIden;//always be 0
	//penalty parameter in semi-smooth Newton algorithm
	double contFact;//the order of Young's modulus of contacting bodies
	bool CONTACT_ANALYSIS(const Matrix<vector<pair<long,double>>,2,1> &loadList, 
		const Matrix<Vector3d,Dynamic,Dynamic> &fricDire);
};

//compare two CONTACT_ELEMENTs by CONTACT_ELEMENT.initGap
int COMPARE_CE(const void *a, const void *b);

/*************************************************************************************************/
/*****************************************IMPLEMENTATION******************************************/
/*************************************************************************************************/

CONTACT::CONTACT(){
	//quadrature point = - sqrt(3.0) / 3.0, quadrature weight = 1.0
	double legeGaus[2][2] = {{- sqrt(3.0) / 3.0, 1.0},{sqrt(3.0) / 3.0, 1.0}};
	//intergral boundary transformation: from [-1, 1] * [-1, 1] to [0, 1] * [0, 1-x]
	//H.T. Rathod, K.V. Nagaraja, B. Venkatesudu, N.L. Ramesh.
	//Gauss Legendre quadrature over a triangle. IIS, 2004, 84: 183-188.
	for(long ti = 0; ti < 2; ti ++){
		for(long tj = 0; tj < 2; tj ++){
			quadTabl(ti*2+tj,0) = (1.0 + legeGaus[ti][0])/2.0;
			quadTabl(ti*2+tj,1) = (1.0 - legeGaus[ti][0]) * (1.0 + legeGaus[tj][0]) / 4.0;
			quadTabl(ti*2+tj,2) = (1.0 - legeGaus[ti][0]) 
				/ 8.0 * legeGaus[ti][1] * legeGaus[tj][1];
		}
	}
	nonmIden = 0;
	contFact = 1.0E11;
}

bool CONTACT::OUTPUT_CONTACT(){
	ofstream tempOfst(DIRECTORY("resuNonmSegm.txt"), ios::out);
	for(long ti = 0; ti < nonmSegm.rows(); ti ++){
		tempOfst << setw(10) << nonmSegm(ti,0)
			<< setw(10) << nonmSegm(ti,1)
			<< setw(10) << nonmSegm(ti,2) << endl;
	}
	tempOfst.close();
	tempOfst.open(DIRECTORY("resuMortSegm.txt"), ios::out);
	for(long ti = 0; ti < mortSegm.rows(); ti ++){
		tempOfst <<setw(10) << mortSegm(ti,0)
			<< setw(10) << mortSegm(ti,1)
			<< setw(10) << mortSegm(ti,2) << endl;
	}
	tempOfst.close();
	return true;
}

bool CONTACT::BUCKET_SORT(const Matrix<double,Dynamic,Dynamic> &locaCoor, 
	const MatrixXd featSize, long diviNumb_0, long diviNumb_1){
	buckLoca(0,0) = locaCoor.block(0,0,locaCoor.rows(),1).minCoeff() - featSize(0);
	buckLoca(0,1) = locaCoor.block(0,0,locaCoor.rows(),1).maxCoeff() + featSize(0);
	buckLoca(1,0) = locaCoor.block(0,1,locaCoor.rows(),1).minCoeff() - featSize(1);
	buckLoca(1,1) = locaCoor.block(0,1,locaCoor.rows(),1).maxCoeff() + featSize(1);
	buckNumb << diviNumb_0, diviNumb_1;
	buck.resize(buckNumb(0), buckNumb(1));
	buckLoca(0,2) = (buckLoca(0,1) - buckLoca(0,0)) / buckNumb(0);
	buckLoca(1,2) = (buckLoca(1,1) - buckLoca(1,0)) / buckNumb(1);
	for(long ti = 0; ti < nonmSegm.rows(); ti ++){
		long buck_r = (locaCoor(ti,0) - buckLoca(0,0)) / buckLoca(0,2);
		long buck_c = (locaCoor(ti,1) - buckLoca(1,0)) / buckLoca(1,2);
		buck(buck_r, buck_c).push_back(ti);
	}
	return true;
}

bool CONTACT::CONTACT_SEARCH(const Matrix<double,Dynamic,Dynamic> &locaCoor){
	vector<CONTACT_ELEMENT> tempCont;
	for(long ti = 0; ti < mortSegm.rows(); ti ++){
		if(ti%1000 == 0){
			cout<<"The "<<ti<<"/"<<mortSegm.rows()<<"-th segment."<<endl;
		}
		long buck_r = (locaCoor(ti,0) - buckLoca(0,0)) / buckLoca(0,2);
		long buck_c = (locaCoor(ti,1) - buckLoca(1,0)) / buckLoca(1,2);
		if(buck_r < 0 || buck_r > buckNumb(0) - 1 || buck_c < 0 || buck_c > buckNumb(1) - 1){
			continue;
		}
		for(list<long>::iterator iter_0 = buck(buck_r, buck_c).begin(); 
			iter_0 != buck(buck_r, buck_c).end(); iter_0 ++){
			TRIANGLE_INTERSECT(*iter_0, ti, tempCont);
		}
		for(long tj = max(buck_r - 1, (long)0); tj <= min(buck_r + 1, buckNumb(0) - 1); tj ++){
			for(long tk = max(buck_c - 1, (long)0); tk <= min(buck_c + 1,buckNumb(1) - 1); tk ++){
				if(tj == buck_r && tk == buck_c){
					continue;
				}
				for(list<long>::iterator iter_0 = buck(tj, tk).begin(); 
					iter_0 != buck(tj, tk).end(); iter_0 ++){
					TRIANGLE_INTERSECT(*iter_0, ti, tempCont);
				}
			}
		}
	}
	//sort by initial gap
	OUTPUT_TIME("Contact segment sort");
	CONTACT_ELEMENT * tempCont_1 = new CONTACT_ELEMENT[tempCont.size()];
	for(long ti = 0; ti < tempCont.size(); ti ++){
		tempCont_1[ti] = tempCont[ti];
	}
	qsort(tempCont_1, tempCont.size(), sizeof(CONTACT_ELEMENT), COMPARE_CE);
	contElem.resize(tempCont.size(),1);
	for(long ti = 0; ti < tempCont.size(); ti ++){
		contElem(ti) = tempCont_1[ti];
	}
	return true;
}

bool CONTACT::TRIANGLE_INTERSECT(long tempNonm, long tempMort, 
	vector<CONTACT_ELEMENT> &tempCont){
	Matrix<Vector3d,3,1> nonmTria, mortTria;
	for(long ti = 0; ti < 3; ti ++){
		nonmTria(ti) = vem(0).nodeCoor.block(nonmSegm(tempNonm,ti),0,1,3).transpose();
		mortTria(ti) = vem(1).nodeCoor.block(mortSegm(tempMort,ti),0,1,3).transpose();
	}
	Vector3d normVect;
	double origArea;
	Matrix<Vector2d,3,1> nonmProj, mortProj;
	vector<Vector2d> tempQuad;
	vector<double> tempWeig;
	TI_SUB(nonmTria, mortTria, normVect, origArea, nonmProj, mortProj, tempQuad, tempWeig);
	//dual basis
	Matrix<double,3,3> matrA, matrD, matrM;
	matrD << 1.0/3.0 * origArea, 0.0, 0.0,
		0.0, 1.0/3.0 * origArea, 0.0,
		0.0, 0.0, 1.0/3.0 * origArea;
	matrM = MatrixXd::Zero(3,3);
	for(long ti = 0; ti < quadTabl.rows(); ti ++){
		matrM(0,0) = matrM(0,0) + quadTabl(ti,2) 
			* quadTabl(ti,0) * quadTabl(ti,0);
		matrM(0,1) = matrM(0,1) + quadTabl(ti,2) 
			* quadTabl(ti,0) * quadTabl(ti,1);
		matrM(0,2) = matrM(0,2) + quadTabl(ti,2) 
			* quadTabl(ti,0) * (1.0 - quadTabl(ti,0) - quadTabl(ti,1));
		matrM(1,0) = matrM(1,0) + quadTabl(ti,2) 
			* quadTabl(ti,1) * quadTabl(ti,0);
		matrM(1,1) = matrM(1,1) + quadTabl(ti,2) 
			* quadTabl(ti,1) * quadTabl(ti,1);
		matrM(1,2) = matrM(1,2) + quadTabl(ti,2) 
			* quadTabl(ti,1) * (1.0 - quadTabl(ti,0) - quadTabl(ti,1));
		matrM(2,0) = matrM(2,0) + quadTabl(ti,2) 
			* (1.0 - quadTabl(ti,0) - quadTabl(ti,1)) * quadTabl(ti,0);
		matrM(2,1) = matrM(2,1) + quadTabl(ti,2) 
			* (1.0 - quadTabl(ti,0) - quadTabl(ti,1)) * quadTabl(ti,1);
		matrM(2,2) = matrM(2,2) + quadTabl(ti,2) 
			* (1.0 - quadTabl(ti,0) - quadTabl(ti,1)) * (1.0 - quadTabl(ti,0) - quadTabl(ti,1));
	}
	matrM = 2.0 * origArea * matrM;
	matrA = matrD * matrM.inverse();
	//shape function, initial gap
	Matrix<double,2,1> triaArea;
	triaArea << TRIANGLE_AREA_2D(nonmProj(0), nonmProj(1), nonmProj(2)), 
		TRIANGLE_AREA_2D(mortProj(0), mortProj(1), mortProj(2));
	for(long ti = 0; ti < tempQuad.size(); ti ++){
		Matrix<double,3,2> tempShap, tempPoin;
		for(long tj = 0; tj < 3; tj ++){
			tempShap((tj+2)%3,0) = 
				TRIANGLE_AREA_2D(nonmProj(tj), nonmProj((tj + 1) % 3), 
					tempQuad[ti]) / triaArea(0);
			tempShap((tj+2)%3,1) = 
				TRIANGLE_AREA_2D(mortProj(tj), mortProj((tj + 1) % 3), 
					tempQuad[ti]) / triaArea(1);
		}
		tempPoin = MatrixXd::Zero(3,2);
		for(long tj = 0; tj < 3; tj ++){
			tempPoin.block(0,0,3,1) = tempPoin.block(0,0,3,1).eval() 
				+ tempShap(tj,0) * nonmTria(tj);
			tempPoin.block(0,1,3,1) = tempPoin.block(0,1,3,1).eval() 
				+ tempShap(tj,1) * mortTria(tj);
		}
		Vector3d nonmMort = tempPoin.block(0,1,3,1).eval() 
			- tempPoin.block(0,0,3,1).eval();
		CONTACT_ELEMENT tempElem;
		for(long tj = 0; tj < 3; tj ++){
			tempElem.node(tj,0) = nonmSegm(tempNonm,tj);
			tempElem.node(tj,1) = mortSegm(tempMort,tj);
		}
		tempElem.shap = tempShap;
		tempElem.dual = matrA * tempShap.block(0,nonmIden,3,1).eval();
		tempElem.contPoin = tempPoin;
		tempElem.normVect = normVect;
		tempElem.initGap = nonmMort.dot(normVect);
		tempElem.quadWeig = tempWeig[ti];
		tempCont.push_back(tempElem);
	}
	return true;
}

bool CONTACT::TI_SUB(const Matrix<Vector3d,3,1> &nonmTria, const Matrix<Vector3d,3,1> &mortTria, 
	Vector3d &normVect, double &origArea, 
	Matrix<Vector2d,3,1> &nonmProj, Matrix<Vector2d,3,1> &mortProj, 
	vector<Vector2d> &tempQuad, vector<double> &tempWeig){
	//projection
	Matrix<Vector3d,2,1> tangVect;
	tangVect(0) = nonmTria(1) - nonmTria(0);
	tangVect(1) = nonmTria(2) - nonmTria(0);
	normVect = tangVect(0).cross(tangVect(1));
	origArea = normVect.norm() / 2.0;
	normVect.normalize();
	Matrix<Vector3d,2,1> unitVect;
	unitVect(0) = tangVect(0).normalized();
	unitVect(1) = normVect.cross(unitVect(0));
	nonmProj(0) << 0.0, 0.0;
	nonmProj(1) << tangVect(0).dot(unitVect(0)), tangVect(0).dot(unitVect(1));
	nonmProj(2) << tangVect(1).dot(unitVect(0)), tangVect(1).dot(unitVect(1));
	for(long ti = 0; ti < 3; ti ++){
		Vector3d tempVect = mortTria(ti) - nonmTria(0);
		mortProj(ti) << tempVect.dot(unitVect(0)), tempVect.dot(unitVect(1));
	}
	//intersection
	vector<Vector2d> tempInte_0, tempInte_1, tempInte_2;
	for(long ti = 0; ti < 3; ti ++){
		if(IN_TRIANGLE_2D(mortProj(ti), nonmProj)){
			tempInte_0.push_back(mortProj(ti));
		}
		if(IN_TRIANGLE_2D(nonmProj(ti), mortProj)){
			tempInte_0.push_back(nonmProj(ti));
		}
	}
	for(long ti = 0; ti < 3; ti ++){
		for(long tj = 0; tj < 3; tj ++){
			LINE_INTERSECT_2D(nonmProj(ti), nonmProj((ti+1)%3), 
				mortProj(tj), mortProj((tj+1)%3), tempInte_0);
		}
	}
	if(tempInte_0.size() < 3){
		return false;
	}
	//no repeat
	Matrix<long,Dynamic,Dynamic> inde;
	inde.resize(tempInte_0.size(),1);
	for(long ti = 0; ti < tempInte_0.size(); ti ++){
		inde(ti) = ti;
	}
	auto rule = [tempInte_0](long ti, long tj)->bool{
		if(tempInte_0[ti](0) < tempInte_0[tj](0) - 1.0E-10){
			return true;
		}
		else if(tempInte_0[ti](0) <= tempInte_0[tj](0) + 1.0E-10){
			if(tempInte_0[ti](1) < tempInte_0[tj](1) - 1.0E-10){
				return true;
			}
			else if(tempInte_0[ti](1) <= tempInte_0[tj](1) + 1.0E-10){
				return false;
			}
			else{
				return false;
			}
		}
		else{
			return false;
		}
	};
	sort(inde.data(), inde.data() + inde.size(), rule);
	tempInte_1.push_back(tempInte_0[inde(0)]);
	for(long ti = 1; ti < tempInte_0.size(); ti ++){
		if(abs(tempInte_0[inde(ti)](0) - tempInte_0[inde(ti-1)](0)) > 1.0E-10
			|| abs(tempInte_0[inde(ti)](1) - tempInte_0[inde(ti-1)](1)) > 1.0E-10){
			tempInte_1.push_back(tempInte_0[inde(ti)]);
		}
	}
	//sort
	Vector2d tempCent = MatrixXd::Zero(2,1);
	for(vector<Vector2d>::iterator iter_0 = tempInte_1.begin(); 
		iter_0 != tempInte_1.end(); iter_0 ++){
		tempCent = tempCent + *iter_0;
	}
	tempCent = tempCent / tempInte_1.size();
	MatrixXd tempAngl;
	tempAngl.resize(tempInte_1.size(),1);
	long tempInde = 0;
	for(vector<Vector2d>::iterator iter_0 = tempInte_1.begin(); 
		iter_0 != tempInte_1.end(); iter_0 ++){
		Vector2d tempVect = *iter_0 - tempCent;
		tempAngl(tempInde) = atan2(tempVect(1), tempVect(0));
		tempInde ++;
	}
	inde.resize(tempInte_1.size(),1);
	for(long ti = 0; ti < tempInte_1.size(); ti ++){
		inde(ti) = ti;
	}
	auto rule_1 = [tempAngl](long ti, long tj)->bool{
		return tempAngl(ti) < tempAngl(tj);
	};
	sort(inde.data(), inde.data()+inde.size(), rule_1);
	for(long ti = 0; ti < tempInte_1.size(); ti ++){
		tempInte_2.push_back(tempInte_1[inde(ti)]);
	}
	vector<double> tempArea;
	if(tempInte_2.size() == 3){
		double area = TRIANGLE_AREA_2D(tempInte_2[0], tempInte_2[1], tempInte_2[2]);
		if(abs(area) <= miniArea){
			return false;
		}
		TRIANGLE_QUADRATURE(tempInte_2[0], tempInte_2[1], tempInte_2[2], 
			tempQuad, tempWeig, tempArea);
	}
	else{
		//centroid and area
		//R. Nurnberg. Calculating the area and centroid of a polygon in 2d.
		//https://paulbourke.net/geometry/polygonmesh/centroid.pdf
		double area = 0.0;
		tempCent = MatrixXd::Zero(2,1);
		long tempSize = tempInte_2.size();
		for(long ti = 0; ti < tempInte_2.size(); ti ++){
			area = area + tempInte_2[ti](0) * tempInte_2[(ti+1)%tempSize](1) 
				- tempInte_2[(ti+1)%tempSize](0) * tempInte_2[ti](1);
			tempCent(0) = tempCent(0) + (tempInte_2[ti](0) + tempInte_2[(ti+1)%tempSize](0)) 
				* (tempInte_2[ti](0) * tempInte_2[(ti+1)%tempSize](1) 
				- tempInte_2[(ti+1)%tempSize](0) * tempInte_2[ti](1));
			tempCent(1) = tempCent(1) + (tempInte_2[ti](1) + tempInte_2[(ti+1)%tempSize](1)) 
				* (tempInte_2[ti](0) * tempInte_2[(ti+1)%tempSize](1) 
				- tempInte_2[(ti+1)%tempSize](0) * tempInte_2[ti](1));
		}
		area /= 2.0;
		if(abs(area) <= miniArea){
			return false;
		}
		tempCent = tempCent / 6.0 / area;
		for(long ti = 0; ti < tempInte_2.size(); ti ++){
			TRIANGLE_QUADRATURE(tempCent, tempInte_2[ti], 
				tempInte_2[(ti + 1) % tempInte_2.size()], tempQuad, tempWeig, tempArea);
		}
	}
	return true;
}

bool CONTACT::TRIANGLE_QUADRATURE(Vector2d tempPoin_0, Vector2d tempPoin_1, Vector2d tempPoin_2, 
		vector<Vector2d> &tempQuad, vector<double> &tempWeig, vector<double> &tempArea){
	double area = TRIANGLE_AREA_2D(tempPoin_0, tempPoin_1, tempPoin_2);
	for(long ti = 0; ti < quadTabl.rows(); ti ++){
		Vector2d tempPoin = quadTabl(ti,0) * tempPoin_0 
			+ quadTabl(ti,1) * tempPoin_1 
			+ (1.0 - quadTabl(ti,0) - quadTabl(ti,1)) * tempPoin_2;
		tempQuad.push_back(tempPoin);
		//integral variable transformation
		//from Cartesian coordinate to shape function
		tempWeig.push_back(2.0 * area * quadTabl(ti,2));
		tempArea.push_back(area);
	}
	return true;
}

bool CONTACT::IS_CROSS_2D(Vector2d tempPoin_0, Vector2d tempPoin_1, 
		Vector2d tempPoin_2, Vector2d tempPoin_3){
	if(max(tempPoin_0(0), tempPoin_1(0)) < min(tempPoin_2(0), tempPoin_3(0))
		|| max(tempPoin_0(1), tempPoin_1(1)) < min(tempPoin_2(1), tempPoin_3(1))
		|| min(tempPoin_0(0), tempPoin_1(0)) > max(tempPoin_2(0), tempPoin_3(0))
		|| min(tempPoin_0(1), tempPoin_1(1)) > max(tempPoin_2(1), tempPoin_3(1))){
		return false;
	}
	if(((tempPoin_2(0) - tempPoin_0(0)) * (tempPoin_2(1) - tempPoin_3(1)) 
		- (tempPoin_2(1) - tempPoin_0(1)) * (tempPoin_2(0) - tempPoin_3(0))) 
		* ((tempPoin_2(0) - tempPoin_1(0)) * (tempPoin_2(1) - tempPoin_3(1)) 
		- (tempPoin_2(1) - tempPoin_1(1)) * (tempPoin_2(0) - tempPoin_3(0))) <= 0
		&& ((tempPoin_0(0) - tempPoin_2(0)) * (tempPoin_0(1) - tempPoin_1(1)) 
		- (tempPoin_0(1) - tempPoin_2(1)) * (tempPoin_0(0) - tempPoin_1(0))) 
		* ((tempPoin_0(0) - tempPoin_3(0)) * (tempPoin_0(1) - tempPoin_1(1)) 
		- (tempPoin_0(1) - tempPoin_3(1)) * (tempPoin_0(0) - tempPoin_1(0))) <= 0){
		return true;
	}
	else{
		return false;
	}
}

bool CONTACT::SORT_BY_2D(Vector2d &tempPoin_0, Vector2d &tempPoin_1, long tempIden){
	if(tempPoin_0(tempIden) > tempPoin_1(tempIden)){
		Vector2d tempPoin;
		tempPoin = tempPoin_0;
		tempPoin_0 = tempPoin_1;
		tempPoin_1 = tempPoin;
	}
	return true;
}

bool CONTACT::LINE_INTERSECT_2D(Vector2d tempPoin_0, Vector2d tempPoin_1, 
		Vector2d tempPoin_2, Vector2d tempPoin_3, vector<Vector2d> &tempInte){
	if(IS_CROSS_2D(tempPoin_0, tempPoin_1, tempPoin_2, tempPoin_3)){
		double area_2 = TRIANGLE_AREA_2D(tempPoin_2, tempPoin_0, tempPoin_1);
		double area_3 = TRIANGLE_AREA_2D(tempPoin_3, tempPoin_0, tempPoin_1);
		if(abs(area_2) < miniArea && abs(area_3) < miniArea){//collinear
			if(abs(tempPoin_0(0) - tempPoin_1(0)) > 1.0E-10){
				SORT_BY_2D(tempPoin_0, tempPoin_1, 0);
				SORT_BY_2D(tempPoin_2, tempPoin_3, 0);
				double from_x = tempPoin_0(0);
				double from_y = tempPoin_0(1);
				if(tempPoin_0(0) < tempPoin_2(0)){
					from_x = tempPoin_2(0);
					from_y = tempPoin_2(1);
				}
				double to_x = tempPoin_1(0);
				double to_y = tempPoin_1(1);
				if(tempPoin_1(0) > tempPoin_3(0)){
					to_x = tempPoin_3(0);
					to_y = tempPoin_3(1);
				}
				Vector2d tempPoin;
				if(abs(from_x - to_x)<1.0E-10){
					tempPoin << from_x, from_y;
					tempInte.push_back(tempPoin);
				}
				else{
					tempPoin << from_x, from_y;
					tempInte.push_back(tempPoin);
					tempPoin << to_x, to_y;
					tempInte.push_back(tempPoin);
				}
			}
			else{
				SORT_BY_2D(tempPoin_0, tempPoin_1, 1);
				SORT_BY_2D(tempPoin_2, tempPoin_3, 1);
				double from_x = tempPoin_0(0);
				double from_y = tempPoin_0(1);
				if(tempPoin_0(1) < tempPoin_2(1)){
					from_x = tempPoin_2(0);
					from_y = tempPoin_2(1);
				}
				double to_x = tempPoin_1(0);
				double to_y = tempPoin_1(1);
				if(tempPoin_1(1) > tempPoin_3(1)){
					to_x = tempPoin_3(0);
					to_y = tempPoin_3(1);
				}
				Vector2d tempPoin;
				if(abs(from_x - to_x)<1.0E-10){
					tempPoin << from_x, from_y;
					tempInte.push_back(tempPoin);
				}
				else{
					tempPoin << from_x, from_y;
					tempInte.push_back(tempPoin);
					tempPoin << to_x, to_y;
					tempInte.push_back(tempPoin);
				}
			}
		}
		else if(abs(area_2) < miniArea){//one endpoint lies on the another line-segment
			tempInte.push_back(tempPoin_2);
		}
		else if(abs(area_3) < miniArea){//one endpoint lies on the another line-segment
			tempInte.push_back(tempPoin_3);
		}
		else{//true intersect
			double tempFact = area_2 / area_3;
			Vector2d tempPoin;
			tempPoin << (tempPoin_2(0) + tempFact*tempPoin_3(0)) / (1.0+tempFact), 
				(tempPoin_2(1) + tempFact*tempPoin_3(1)) / (1.0+tempFact);
			tempInte.push_back(tempPoin);
		}
	}
	return true;
}

bool CONTACT::IN_TRIANGLE_2D(Vector2d tempPoin, Matrix<Vector2d,3,1> tempTria){
	double area_0 = TRIANGLE_AREA_2D(tempPoin, tempTria(0), tempTria(1));
	double area_1 = TRIANGLE_AREA_2D(tempPoin, tempTria(1), tempTria(2));
	double area_2 = TRIANGLE_AREA_2D(tempPoin, tempTria(2), tempTria(0));
	double areaTota = TRIANGLE_AREA_2D(tempTria(0), tempTria(1), tempTria(2));
	if(area_0+area_1+area_2 <= areaTota+miniArea){
		return true;
	}
	else{
		return false;
	}
}

double CONTACT::TRIANGLE_AREA_2D(Vector2d tempPoin_0, Vector2d tempPoin_1, Vector2d tempPoin_2){
	Vector2d tempVect_0 = tempPoin_1 - tempPoin_0;
	Vector2d tempVect_1 = tempPoin_2 - tempPoin_0;
	return abs(tempVect_0(0) * tempVect_1(1) - tempVect_0(1) * tempVect_1(0)) / 2.0;
}

bool CONTACT::OUTPUT_CONTACT_ELEMENT(){
	ofstream tempOfst(DIRECTORY("resuContElem.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	for(long ti = 0; ti < contElem.rows(); ti ++){
		tempOfst<<setw(30)<<contElem(ti).contPoin(0,0)
			<<setw(30)<<contElem(ti).contPoin(1,0)
			<<setw(30)<<contElem(ti).contPoin(2,0)
			<<setw(30)<<contElem(ti).contPoin(0,1)
			<<setw(30)<<contElem(ti).contPoin(1,1)
			<<setw(30)<<contElem(ti).contPoin(2,1)
			<<setw(30)<<contElem(ti).initGap
			<<setw(30)<<contElem(ti).quadWeig<<endl;
	}
	tempOfst.close();
	return true;
}

bool CONTACT::CONTACT_ANALYSIS(const Matrix<vector<pair<long,double>>,2,1> &loadList, 
	const Matrix<Vector3d,Dynamic,Dynamic> &fricDire){
	//stiffness matrix
	for(long tv = 0; tv < 2; tv ++){
		vem(tv).ELEMENT_STIFF();
	}
	//all non-mortar nodes: active
	map<long,long> isAI;
	map<long,long> contAAIN;//contact, all active and inactive
	Matrix<Vector3d,Dynamic,Dynamic> tracVect;
	tracVect.resize(contElem.rows(),1);
	for(long ti = 0; ti < contElem.rows(); ti ++){
		for(long tj = 0; tj <3; tj ++){
			isAI.insert(map<long,long>::value_type(contElem(ti).node(tj,nonmIden), 1));
			map<long,long>::iterator iter_O = contAAIN.find(contElem(ti).node(tj,nonmIden));
			if(iter_O == contAAIN.end()){
				long tempSize = contAAIN.size();
				contAAIN.insert(map<long,long>::value_type(
					contElem(ti).node(tj,nonmIden), tempSize));
			}
		}
	}
	//non-mortar segment: contact integration point
	OUTPUT_TIME("Boundary consistent treatment");
	map<MATRIX,vector<long>> segmPoin;
	for(long ti = 0; ti < contElem.rows(); ti ++){
		MATRIX tempMatr;
		tempMatr.resize(3,1);
		//MATRIX: no need to sort the components in tempMatr
		tempMatr(0) = contElem(ti).node(0,nonmIden);
		tempMatr(1) = contElem(ti).node(1,nonmIden);
		tempMatr(2) = contElem(ti).node(2,nonmIden);
		map<MATRIX,vector<long>>::iterator iter_S = segmPoin.find(tempMatr);
		if(iter_S == segmPoin.end()){
			vector<long> tempVect;
			tempVect.push_back(ti);
			segmPoin.insert(map<MATRIX,vector<long>>::value_type(tempMatr, tempVect));
		}
		else{
			(iter_S->second).push_back(ti);
		}
	}
	for(map<MATRIX,vector<long>>::iterator iter_S = segmPoin.begin(); 
		iter_S != segmPoin.end(); iter_S ++){
		//boundary consistent treatment, partially integrated non-mortar segment
		Matrix3d tempD = MatrixXd::Zero(3,3);
		Matrix3d tempM = MatrixXd::Zero(3,3);
		for(long tj = 0; tj < (iter_S->second).size(); tj ++){
			CONTACT_ELEMENT tempElem = contElem((iter_S->second)[tj]);
			Matrix3d tempD_j, tempM_j;
			tempD_j = MatrixXd::Zero(3,3);
			for(long tk = 0; tk < 3; tk ++){
				tempD_j(tk,tk) = tempElem.shap(tk,nonmIden);
				for(long tm = 0; tm < 3; tm ++){
					tempM_j(tk,tm) = tempElem.shap(tk,nonmIden) * tempElem.shap(tm,nonmIden);
				}
			}
			tempD_j = tempElem.quadWeig * tempD_j;
			tempM_j = tempElem.quadWeig * tempM_j;
			tempD = tempD + tempD_j;
			tempM = tempM + tempM_j;
		}
		Matrix3d tempA = tempD * tempM.inverse();
		for(long tj = 0; tj < (iter_S->second).size(); tj ++){
			contElem((iter_S->second)[tj]).dual = tempA 
				* contElem((iter_S->second)[tj]).shap.block(0,nonmIden,3,1);
		}
		//average traction vector
		if(fricDire.rows() > 0){
			Vector3d tempTrac = MatrixXd::Zero(3,1);
			double quadSumm = 0.0;
			for(long tj = 0; tj < (iter_S->second).size(); tj ++){
				CONTACT_ELEMENT tempElem = contElem((iter_S->second)[tj]);
				tempTrac = tempTrac + tempElem.quadWeig * fricDire((iter_S->second)[tj]);
				quadSumm += tempElem.quadWeig;
			}
			tempTrac = tempTrac / quadSumm;
			for(long tj = 0; tj < (iter_S->second).size(); tj ++){
				tracVect((iter_S->second)[tj]) = tempTrac;
			}
		}
	}
	segmPoin.clear();
	//
	set<long> oppoReco;
	while(true){
		map<long,long> nnan;//non-mortar node, active contact number
		for(long ti = 0; ti < contElem.rows(); ti ++){
			for(long tj = 0; tj < 3; tj ++){
				map<long,long>::iterator iter_A = isAI.find(contElem(ti).node(tj,nonmIden));
				map<long,long>::iterator iter_C = nnan.find(contElem(ti).node(tj,nonmIden));
				if(iter_A->second == 1 && iter_C == nnan.end()){
					long tempSize = nnan.size();
					nnan.insert(map<long,long>::value_type(
						contElem(ti).node(tj,nonmIden), tempSize));
				}
			}
		}
		//total non-mortar and mortar DOF
		long tonmFree = vem(0).freeNumb + vem(1).freeNumb + nnan.size();
		//load
		MatrixXd origForc = MatrixXd::Zero(tonmFree,1);
		for(vector<pair<long,double>>::const_iterator iter_L = loadList(0).begin();
			iter_L != loadList(0).end(); iter_L ++){
			origForc(vem(0).consDOF(iter_L->first)) = 
				origForc(vem(0).consDOF(iter_L->first)) + iter_L->second;
		}
		for(vector<pair<long,double>>::const_iterator iter_L = loadList(1).begin();
			iter_L != loadList(1).end(); iter_L ++){
			origForc(vem(0).freeNumb + vem(1).consDOF(iter_L->first)) = 
				origForc(vem(0).freeNumb + vem(1).consDOF(iter_L->first)) + iter_L->second;
		}
		//system matrix, n, m, c
		SparseMatrix<double, RowMajor> origStif;
		//displacement TO contact penetration
		//Lagrange multiplier to sliding friction force
		vector<Triplet<double>> tranD2cp, tranM2sf, origList;
		MatrixXd initClea = MatrixXd::Zero(contAAIN.size(),1);
		for(long ti = 0; ti < 2; ti ++){
			long baseInde_0 = 0, baseInde_1 = 0;
			if(ti == 1){
				baseInde_0 = vem(0).freeNumb;
				baseInde_1 = vem(0).freeNumb;
			}
			for(vector<Triplet<double>>::const_iterator iter_0 = vem(ti).stifList.begin(); 
				iter_0 != vem(ti).stifList.end(); iter_0 ++){
				origList.emplace_back(baseInde_0 + iter_0->row(), baseInde_1 + iter_0->col(),
					iter_0->value());
			}
		}
		for(long ti = 0; ti < contElem.rows(); ti ++){
			Matrix<double,3,9> T_e_0, T_e_1;
			Matrix<double,3,1> G_e;
			Matrix<double,3,1> M_e = contElem(ti).dual;
			Matrix<double,3,2> tempShap;
			tempShap << contElem(ti).shap(0,nonmIden), contElem(ti).shap(0,1-nonmIden), 
				contElem(ti).shap(1,nonmIden), contElem(ti).shap(1,1-nonmIden), 
				contElem(ti).shap(2,nonmIden), contElem(ti).shap(2,1-nonmIden);
			Matrix<double,3,9> N_e_0, N_e_1;
			N_e_0 << 
				tempShap(0,0), 0.0, 0.0, tempShap(1,0), 0.0, 0.0, tempShap(2,0), 0.0, 0.0,
				0.0, tempShap(0,0), 0.0, 0.0, tempShap(1,0), 0.0, 0.0, tempShap(2,0), 0.0, 
				0.0, 0.0, tempShap(0,0), 0.0, 0.0, tempShap(1,0), 0.0, 0.0, tempShap(2,0);
			N_e_1 << 
				tempShap(0,1), 0.0, 0.0, tempShap(1,1), 0.0, 0.0, tempShap(2,1), 0.0, 0.0,
				0.0, tempShap(0,1), 0.0, 0.0, tempShap(1,1), 0.0, 0.0, tempShap(2,1), 0.0, 
				0.0, 0.0, tempShap(0,1), 0.0, 0.0, tempShap(1,1), 0.0, 0.0, tempShap(2,1);
			//
			T_e_0 = M_e * contElem(ti).normVect.transpose() * N_e_0 * contElem(ti).quadWeig;
			T_e_1 = M_e * contElem(ti).normVect.transpose() * N_e_1 * contElem(ti).quadWeig;
			G_e = M_e * contElem(ti).initGap * contElem(ti).quadWeig;
			Matrix<double,9,3> F_e_0 = MatrixXd::Zero(9,3);
			Matrix<double,9,3> F_e_1 = MatrixXd::Zero(9,3);
			if(fricDire.rows() > 0){
				F_e_0 = N_e_0.transpose() * tracVect(ti) 
					* M_e.transpose() * contElem(ti).quadWeig;
				F_e_1 = N_e_1.transpose() * tracVect(ti) 
					* M_e.transpose() * contElem(ti).quadWeig;
			}
			for(long tj = 0; tj < 3; tj ++){
				//all active and inactive
				map<long,long>::iterator iter_O = contAAIN.find(contElem(ti).node(tj,nonmIden));
				long row_tj = iter_O->second;
				for(long tk = 0; tk < 3; tk ++){
					long nonm_tk = contElem(ti).node(tk,nonmIden);
					long mort_tk = contElem(ti).node(tk,1-nonmIden);
					for(long tm = 0; tm < 3; tm ++){
						if(tj == tk){
							if(vem(0).consDOF(3*nonm_tk+tm) != -1){
								tranD2cp.emplace_back(row_tj, 
									0 + vem(0).consDOF(3 * nonm_tk + tm), 
									- T_e_0(tj, 3 * tk + tm));
								tranM2sf.emplace_back(vem(0).consDOF(3 * nonm_tk + tm), 
									row_tj, F_e_0(3 * tk + tm, tj));
							}
						}
						if(vem(1).consDOF(3*mort_tk+tm) != -1){
							tranD2cp.emplace_back(row_tj, 
								vem(0).freeNumb + vem(1).consDOF(3 * mort_tk + tm), 
								T_e_1(tj, 3 * tk + tm));
							tranM2sf.emplace_back(
								vem(0).freeNumb + vem(1).consDOF(3 * mort_tk + tm), 
								row_tj, - F_e_1(3 * tk + tm, tj));
						}
					}
				}
				initClea(row_tj) = initClea(row_tj) + G_e(tj);
				//only active
				map<long,long>::iterator iter_C = nnan.find(contElem(ti).node(tj,nonmIden));
				if(iter_C == nnan.end()){
					continue;
				}
				row_tj = iter_C->second;
				for(long tk = 0; tk < 3; tk ++){
					long nonm_tk = contElem(ti).node(tk,nonmIden);
					long mort_tk = contElem(ti).node(tk,1-nonmIden);
					for(long tm = 0; tm < 3; tm ++){
						if(tj == tk){
							if(vem(0).consDOF(3*nonm_tk+tm) != -1){
								//displacement TO contact penetration
								origList.emplace_back(vem(0).freeNumb + vem(1).freeNumb + row_tj, 
									0 + vem(0).consDOF(3 * nonm_tk + tm), T_e_0(tj, 3 * tk + tm));
								//multiplier TO normal contact force and sliding friction force
								origList.emplace_back(0 + vem(0).consDOF(3 * nonm_tk + tm), 
									vem(0).freeNumb + vem(1).freeNumb + row_tj, 
									T_e_0(tj, 3 * tk + tm) + F_e_0(3 * tk + tm, tj));
							}
						}
						if(vem(1).consDOF(3*mort_tk+tm) != -1){
							//displacement TO contact penetration
							origList.emplace_back(vem(0).freeNumb + vem(1).freeNumb + row_tj, 
								vem(0).freeNumb + vem(1).consDOF(3 * mort_tk + tm), 
								- T_e_1(tj, 3 * tk + tm));
							//multiplier TO normal contact force and sliding friction force
							origList.emplace_back(
								vem(0).freeNumb + vem(1).consDOF(3 * mort_tk + tm), 
								vem(0).freeNumb + vem(1).freeNumb + row_tj, 
								- T_e_1(tj, 3 * tk + tm) - F_e_1(3 * tk + tm, tj));
						}
					}
				}
				origForc(vem(0).freeNumb + vem(1).freeNumb + row_tj) = 
					origForc(vem(0).freeNumb + vem(1).freeNumb + row_tj) 
					+ G_e(tj);
			}
		}
		origStif.resize(tonmFree, tonmFree);
		origStif.setFromTriplets(origList.begin(), origList.end());
		origList.clear();
		//system matrix, 0, 1, c
		vector<Triplet<double>> listK_00, listK_01, listK_10, listK_11, 
			listT_0, listT_f0, inveListT_1, inveListT_f1;
		map<long,long> rman;//row maximum, active contact number
		for(map<long,long>::iterator iter_C = nnan.begin(); 
			iter_C != nnan.end(); iter_C ++){
			long tempNode = iter_C->first;
			long tempIden = vem(0).freeNumb + vem(1).freeNumb + iter_C->second;
			double maxiCoef = 0.0;
			long maxiIden;
			for(long ti = 0; ti < 3; ti ++){
				long tempNdof = vem(nonmIden).consDOF(3 * tempNode + ti);
				if(tempNdof != -1 && abs(origStif.coeff(tempIden, tempNdof)) > abs(maxiCoef)){
					maxiCoef = origStif.coeff(tempIden, tempNdof);
					maxiIden = tempNdof;
				}
			}
			rman.insert(map<long,long>::value_type(maxiIden, iter_C->second));
			inveListT_1.emplace_back(iter_C->second, iter_C->second, 1.0 / maxiCoef);
			inveListT_f1.emplace_back(iter_C->second, iter_C->second, 
				1.0 / origStif.coeff(maxiIden, tempIden));
		}
		Matrix<long,Dynamic,Dynamic> freeFlag, nm00Cast;//non-mortar,mortar TO 0,0
		freeFlag.resize(vem(0).freeNumb + vem(1).freeNumb, 1);
		nm00Cast.resize(vem(0).freeNumb + vem(1).freeNumb, 1);
		for(long ti = 0; ti < vem(0).freeNumb + vem(1).freeNumb; ti ++){
			freeFlag(ti) = 1;
			nm00Cast(ti) = 0;
		}
		for(map<long,long>::iterator iter_D = rman.begin(); iter_D != rman.end(); iter_D ++){
			freeFlag(iter_D->first) = 0;
		}
		for(long ti = 1; ti < vem(0).freeNumb + vem(1).freeNumb; ti ++){
			nm00Cast(ti) = nm00Cast(ti-1) + freeFlag(ti-1);
		}
		long tota00 = nm00Cast(vem(0).freeNumb + vem(1).freeNumb - 1) 
			+ freeFlag(vem(0).freeNumb + vem(1).freeNumb - 1);
		MatrixXd resuForc = MatrixXd::Zero(tota00,1);
		MatrixXd rmanForc = MatrixXd::Zero(nnan.size(),1);
		MatrixXd lagrForc = origForc.block(vem(0).freeNumb + vem(1).freeNumb,0,nnan.size(),1);
		for(long ti = 0; ti < vem(0).freeNumb + vem(1).freeNumb; ti ++){
			map<long,long>::iterator iter_D = rman.find(ti);
			for(SparseMatrix<double,RowMajor>::InnerIterator iter_S(origStif, ti); 
				iter_S; ++ iter_S){
				if(iter_S.col() >= vem(0).freeNumb + vem(1).freeNumb){
					break;
				}
				map<long,long>::iterator iter_D1 = rman.find(iter_S.col());
				if(iter_D != rman.end()){
					if(iter_D1 != rman.end()){
						listK_11.emplace_back(iter_D->second, 
							iter_D1->second, iter_S.value());
					}
					else{
						listK_10.emplace_back(iter_D->second, 
							nm00Cast(iter_S.col()), iter_S.value());
					}
				}
				else{
					if(iter_D1 != rman.end()){
						listK_01.emplace_back(nm00Cast(ti), 
							iter_D1->second, iter_S.value());
					}
					else{
						listK_00.emplace_back(nm00Cast(ti), 
							nm00Cast(iter_S.col()), iter_S.value());
					}
				}
			}
			if(iter_D != rman.end()){
				rmanForc(iter_D->second) = origForc(ti);
			}
			else{
				resuForc(nm00Cast(ti)) = origForc(ti);
			}
		}
		for(long ti = 0; ti < nnan.size(); ti ++){
			long row_ti = vem(0).freeNumb + vem(1).freeNumb + ti;
			for(SparseMatrix<double,RowMajor>::InnerIterator iter_S(origStif, row_ti); 
				iter_S; ++ iter_S){
				if(iter_S.col() >= vem(0).freeNumb + vem(1).freeNumb){
					break;
				}
				map<long,long>::iterator iter_D = rman.find(iter_S.col());
				if(iter_D == rman.end()){
					listT_0.emplace_back(ti, 
						nm00Cast(iter_S.col()), iter_S.value());
					listT_f0.emplace_back(nm00Cast(iter_S.col()), 
						ti, origStif.coeff(iter_S.col(), row_ti));
				}
			}
		}
		origStif.resize(0,0);
		//only RowMajor supports parallelization of BiCGSTAB (the colMajor not)
		SparseMatrix<double, RowMajor> K_00, K_01, K_10, K_11, T_0, T_f0, 
			inveT_1, inveT_f1, resuStif;
		K_00.resize(tota00, tota00);
		K_01.resize(tota00, nnan.size());
		K_10.resize(nnan.size(), tota00);
		K_11.resize(nnan.size(), nnan.size());
		T_0.resize(nnan.size(), tota00);
		T_f0.resize(tota00, nnan.size());
		inveT_1.resize(nnan.size(), nnan.size());
		inveT_f1.resize(nnan.size(), nnan.size());
		K_00.setFromTriplets(listK_00.begin(), listK_00.end());
		listK_00.clear();
		K_01.setFromTriplets(listK_01.begin(), listK_01.end());
		listK_01.clear();
		K_10.setFromTriplets(listK_10.begin(), listK_10.end());
		listK_10.clear();
		K_11.setFromTriplets(listK_11.begin(), listK_11.end());
		listK_11.clear();
		T_0.setFromTriplets(listT_0.begin(), listT_0.end());
		listT_0.clear();
		T_f0.setFromTriplets(listT_f0.begin(), listT_f0.end());
		listT_f0.clear();
		inveT_1.setFromTriplets(inveListT_1.begin(), inveListT_1.end());
		inveListT_1.clear();
		inveT_f1.setFromTriplets(inveListT_f1.begin(), inveListT_f1.end());
		inveListT_f1.clear();
		resuStif = K_00 - K_01 * inveT_1 * T_0 - T_f0 * inveT_f1 * K_10
			+ T_f0 * inveT_f1 * K_11 * inveT_1 * T_0;
		resuStif.makeCompressed();
		K_00.resize(0,0);
		resuForc = resuForc - K_01 * inveT_1 * lagrForc - T_f0 * inveT_1 * rmanForc
			+ T_f0 * inveT_f1 * K_11 * inveT_1 * lagrForc;
		//solving
		Matrix<double,Dynamic,Dynamic> resuDisp;
		if(fricDire.rows() > 0){
			BiCGSTAB<SparseMatrix<double, RowMajor>> solv;
			solv.setMaxIterations(10000);
			// solv.setTolerance(1.0E-10);
			solv.compute(resuStif);
			OUTPUT_TIME("Solver");
			resuDisp = solv.solve(resuForc);
			cout<<"#Iterations:      "<<solv.iterations()<<endl;
			cout<<"#Estimated error: "<<solv.error()<<endl;
			if(solv.info() != Success){
#ifdef _POSIX_VERSION
				OUTPUT_TIME("Solution failed, try IncompleteLUT preconditioner in Eigen");
				BiCGSTAB<SparseMatrix<double, RowMajor>, IncompleteLUT<double>> solv_1;
				OUTPUT_TIME("Factorization (TIME CONSUMING!)");
				solv_1.compute(resuStif);
				OUTPUT_TIME("Solver");
				resuDisp = solv_1.solve(resuForc);
				cout<<"#Iterations:      "<<solv_1.iterations()<<endl;
				cout<<"#Estimated error: "<<solv_1.error()<<endl;
#elif __linux__
				OUTPUT_TIME("Solution failed, try KSPBICG in PETSc");
				Mat stif;
				Vec forc, disp;
				KSP ksp;
				PC pc;
				const PetscScalar *arra;
				PetscViewerAndFormat *vf;
				OUTPUT_TIME("Stiffness matrix assemble.");
				PetscCall(MatCreate(PETSC_COMM_WORLD, &stif));
				PetscCall(MatSetSizes(stif, PETSC_DECIDE, PETSC_DECIDE, 
					resuStif.rows(), resuStif.cols()));
				PetscCall(MatSetFromOptions(stif));
				PetscCall(MatSetUp(stif));
				ASSEMBLE(resuStif, stif);
				PetscCall(MatAssemblyBegin(stif, MAT_FINAL_ASSEMBLY));
				PetscCall(MatAssemblyEnd(stif, MAT_FINAL_ASSEMBLY));
				//load
				PetscCall(VecCreate(PETSC_COMM_WORLD, &forc));
				PetscCall(VecSetSizes(forc, PETSC_DECIDE, resuForc.rows()));
				PetscCall(VecSetFromOptions(forc));
				for(int ti = 0; ti < resuForc.rows(); ti ++){
					PetscInt row = ti;
					PetscScalar val = resuForc(ti);
					PetscCall(VecSetValues(forc, 1, &row, &val, INSERT_VALUES));
				}
				PetscCall(VecAssemblyBegin(forc));
				PetscCall(VecAssemblyEnd(forc));
				PetscCall(VecDuplicate(forc, &disp));
				//solve
				OUTPUT_TIME("Preconditioner.");
				PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
				PetscCall(KSPSetOperators(ksp, stif, stif));
				PetscCall(KSPGetPC(ksp, &pc));
				PetscCall(PCSetType(pc, PCILU));//PCJACOBI
				PetscCall(KSPSetType(ksp, KSPBICG));//KSPGMRES
				PetscCall(PCFactorSetReuseOrdering(pc, PETSC_TRUE));
				PetscCall(KSPSetTolerances(ksp, 2.0E-6, 
					PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
				PetscCall(PetscViewerAndFormatCreate(
					PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf));
				PetscCall(KSPMonitorSet(ksp,
					(PetscErrorCode (*)(KSP, PetscInt, PetscReal, void *))KSPMonitorResidual,
					vf, (PetscErrorCode(*)(void **))PetscViewerAndFormatDestroy));
				PetscCall(KSPSetFromOptions(ksp));
				PetscCall(KSPSetUp(ksp));
				OUTPUT_TIME("Solver");
				PetscCall(KSPSolve(ksp, forc, disp));
				PetscCall(VecGetArrayRead(disp, &arra));
				resuDisp.resize(resuForc.rows(), 1);
				for(int ti = 0; ti < resuForc.rows(); ti ++){
					resuDisp(ti) = arra[ti];
				}
				PetscCall(VecRestoreArrayRead(disp, &arra));
				//
				PetscCall(KSPDestroy(&ksp));
				PetscCall(MatDestroy(&stif));
				PetscCall(VecDestroy(&forc));
				PetscCall(VecDestroy(&disp));
#else
				cout<<"Error: unknown system!!!"<<endl;
#endif
			}
		}
		else{
			ConjugateGradient<SparseMatrix<double>, Lower|Upper> solv;
			solv.compute(resuStif);
			OUTPUT_TIME("Solver");
			resuDisp = solv.solve(resuForc);
			cout<<"#Iterations:      "<<solv.iterations()<<endl;
			cout<<"#Estimated error: "<<solv.error()<<endl;
		}
		OUTPUT_TIME("Done");
		resuStif.resize(0,0);
		//original displacement
		MatrixXd rmanDisp = inveT_1 * lagrForc - inveT_1 * T_0 * resuDisp;
		MatrixXd origDisp = MatrixXd::Zero(vem(0).freeNumb + vem(1).freeNumb, 1);
		for(long ti = 0; ti < vem(0).freeNumb + vem(1).freeNumb; ti ++){
			map<long,long>::iterator iter_D = rman.find(ti);
			if(iter_D != rman.end()){
				origDisp(ti) = rmanDisp(iter_D->second);
			}
			else{
				origDisp(ti) = resuDisp(nm00Cast(ti));
			}
		}
		for(long tv = 0; tv < 2; tv ++){
			MatrixXd tempDisp;
			if(tv == 0){
				tempDisp = origDisp.block(0, 0, vem(0).freeNumb, 1);
			}
			else{
				tempDisp = origDisp.block(vem(0).freeNumb, 0, vem(1).freeNumb, 1);
			}
			vem(tv).OUTPUT_DISPLACEMENT(tempDisp, tv);
		}
		//all active and inactive non-mortar nodes
		SparseMatrix<double> d2cpMatr;
		d2cpMatr.resize(contAAIN.size(), vem(0).freeNumb + vem(1).freeNumb);
		d2cpMatr.setFromTriplets(tranD2cp.begin(), tranD2cp.end());
		tranD2cp.clear();
		MatrixXd newtDefo = d2cpMatr * origDisp;
		MatrixXd newtClea = newtDefo + initClea;
		MatrixXd newtLagr = MatrixXd::Zero(contAAIN.size(), 1);
		MatrixXd actiLagr = inveT_f1 * rmanForc - inveT_f1 * K_11 * inveT_1 * lagrForc 
			- inveT_f1 * K_10 * resuDisp + inveT_f1 *K_11 * inveT_1 * T_0 * resuDisp;
		for(map<long,long>::iterator iter_C = nnan.begin(); 
			iter_C != nnan.end(); iter_C ++){
			map<long,long>::iterator iter_O = contAAIN.find(iter_C->first);
			newtLagr(iter_O->second) = actiLagr(iter_C->second);
		}
		//contact force
		MatrixXd newtForc = d2cpMatr.transpose() * newtLagr;
		d2cpMatr.resize(0,0);
		//sliding friction force
		SparseMatrix<double> m2sfMatr;
		m2sfMatr.resize(vem(0).freeNumb + vem(1).freeNumb, contAAIN.size());
		m2sfMatr.setFromTriplets(tranM2sf.begin(), tranM2sf.end());
		tranM2sf.clear();
		MatrixXd slidFric = - m2sfMatr * newtLagr;
		m2sfMatr.resize(0,0);
		ofstream tempOfst(DIRECTORY("resuCont.txt"), ios::out);
		tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
		for(map<long,long>::iterator iter_O = contAAIN.begin(); 
			iter_O != contAAIN.end(); iter_O ++){
			tempOfst<<setw(10)<<iter_O->first
				<<setw(30)<<newtLagr(iter_O->second)
				<<setw(30)<<newtDefo(iter_O->second)
				<<setw(30)<<initClea(iter_O->second)
				<<setw(30)<<newtClea(iter_O->second);
			for(long ti = 0; ti < 3; ti ++){
				if(vem(0).consDOF(3 * (iter_O->first) +ti) != -1){
					tempOfst<<setw(30)<<newtForc(vem(0).consDOF(3 * (iter_O->first) +ti));
				}
				else{
					tempOfst<<setw(30)<<0.0;
				}
			}
			for(long ti = 0; ti < 3; ti ++){
				if(vem(0).consDOF(3 * (iter_O->first) +ti) != -1){
					tempOfst<<setw(30)<<slidFric(vem(0).consDOF(3 * (iter_O->first) +ti));
				}
				else{
					tempOfst<<setw(30)<<0.0;
				}
			}
			tempOfst<<endl;
		}
		tempOfst.close();
		//semi-smooth Newton algorithm
		long oppoNumb = 0;
		for(map<long,long>::iterator iter_O = contAAIN.begin(); 
			iter_O != contAAIN.end(); iter_O ++){
			double tempCrit = newtLagr(iter_O->second) - contFact * newtClea(iter_O->second);
			map<long,long>::iterator iter_A = isAI.find(iter_O->first);
			if(tempCrit > 0.0){
				if(iter_A->second == 0){
					oppoNumb ++;
				}
			}
			else{
				if(iter_A->second == 1){
					oppoNumb ++;
				}
			}
		}
		cout<<"There are "<<oppoNumb<<" opposite multipliers."<<endl;
		if(oppoNumb == 0){
			break;
		}
		for(map<long,long>::iterator iter_O = contAAIN.begin(); 
			iter_O != contAAIN.end(); iter_O ++){
			double tempCrit = newtLagr(iter_O->second) - contFact * newtClea(iter_O->second);
			map<long,long>::iterator iter_A = isAI.find(iter_O->first);
			if(tempCrit > 0.0){
					iter_A->second = 1;
			}
			else{
				iter_A->second = 0;
			}
		}
		oppoReco.insert(oppoNumb);
	}
	return true;
}

int COMPARE_CE(const void *a, const void *b){
	if((*(CONTACT_ELEMENT *)b).initGap - (*(CONTACT_ELEMENT *)a).initGap <0.0){
		return 1;
	}
	else if((*(CONTACT_ELEMENT *)b).initGap - (*(CONTACT_ELEMENT *)a).initGap >0.0){
		return -1;
	}
	else{
		return 0;
	}
}

#endif
