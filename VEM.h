
#ifndef _VEM_H
#define _VEM_H

#include "PREP.h"

/*************************************************************************************************/
/*******************************************DECLARATION*******************************************/
/*************************************************************************************************/

class COORDINATE;

class MATRIX;

class VEM{
public:
	double mateElas;//material elasticity
	double matePois;//material Poisson ratio
	//element
	//n*3: n nodes, 3 coordinates per node
	Matrix<double,Dynamic,Dynamic> nodeCoor;
	//f*3;e*1: f faces of the i-th lement, 3 nodes per face, e elements
	Matrix<Matrix<long,Dynamic,Dynamic>,Dynamic,Dynamic> elemFace;
	//3n*1: the 1..3-th degree of the 1..n-th node
	Matrix<long,Dynamic,Dynamic> consFlag;
	//3n*1: map from node DOF to matrix row/column numbering
	Matrix<long,Dynamic,Dynamic> consDOF;
	//n*1: 3*3 rotation matrix from global coordinate to local coordinate
	Matrix<Matrix3d,Dynamic,Dynamic> freeTran;
	//the total number of unconstrained DOF
	long freeNumb;
	//row numbering + column numbering + stiffness value
	vector<Triplet<double>> stifList;
	//volume
	//from point numbering to coordinate
	map<long,COORDINATE> poinCoor;
	//from coordinate to point numbering
	map<COORDINATE,long> coorPoin;
	long lineCoun;//total number of lines
	//two points forms one line: 
	//point numberings are ascending - positive line numbering
	//point numberings are descending - negative line numbering
	//from line numbering to two end point numberings
	map<long,MATRIX> linePoin;
	//from two end point numberings to line numbering
	map<MATRIX,long> poinLine;
	long faceCoun;//total number of faces
	//from face numbering to line numberings
	map<long,MATRIX> faceLine;
	//from line numberings to face numbering
	map<MATRIX,long> lineFace;
	//line has been used by which faces
	map<long,list<long>> lineUsed;
	long voluCoun;//total number of volumes
	//from volume numbering to face numberings
	map<long,MATRIX> voluFace;
	//face has been used by which volumes
	map<long,list<long>> faceUsed;
	//
	VEM();
	bool CLEAR();
	//try to add one point, the added point may already exist
	long TRY_ADD_POINT(COORDINATE tempCoor);
	//try to add one line, the added line may already exist
	long TRY_ADD_LINE(MATRIX tempMATR);
	//try to add one face, the added face may already exist
	long TRY_ADD_FACE(MATRIX tempMATR);
	long ADD_VOLUME(MATRIX tempMATR);//add one volume
	bool ADD_FACEUSED(long tempFace, long tempVolu);//adjust faceUsed
	bool ADD_LINEUSED(long tempLine, long tempFace);//adjust lineUsed
	//add one volume from 8 corner points
	long ADD_BLOCK_FROM_POINT(const Matrix<long,8,1> & tempCorn);
	bool DELETE_VOLUME(long tempVolu);//delete the volume
	//refine the volume in spliVolu_0
	//spliFlag: contains specified sub volumes after refinement
	//planSurf: Cartesian - curvilinear - Cartesian
	bool REFINE(list<long> &spliVolu_0, const Matrix<long,Dynamic,Dynamic> &spliFlag, 
		const map<MATRIX,COORDINATE> &planSurf);
	//output volume information
	bool OUTPUT_VOLUME(long fileIden);
	//transfer volume to element
	bool VOLUME_2_ELEMENT();
	//output element information
	bool OUTPUT_ELEMENT(long fileIden);
	//calculate element stiffness matrix into stifList
	bool ELEMENT_STIFF();
	//output displacement result
	bool OUTPUT_DISPLACEMENT(Matrix<double,Dynamic,Dynamic> tempDisp, long fileIden);
	//zero displacement constraint
	bool CONSTRAINT(const vector<long> & consFree);
	//DOF coupling: specified DOFs will have a same value
	bool COUPLE(set<long> coupFree);
};

/*************************************************************************************************/
/*****************************************IMPLEMENTATION******************************************/
/*************************************************************************************************/

//3d coordinate, two coordinates can be compared
class COORDINATE:public Vector3d{
public:
	bool operator<(const COORDINATE & coor_1) const{
		for(long ti = 0; ti < 3; ti ++){
			if((*this)(ti) < coor_1(ti) - 1.0E-10) return true;
			else if((*this)(ti) > coor_1(ti) + 1.0E-10) return false;
		}
		return false;
	}
	COORDINATE operator+(const COORDINATE & coor_1) const{
		COORDINATE coor_2;
		for(long ti = 0; ti < 3; ti ++){
			coor_2(ti) = (*this)(ti) + coor_1(ti);
		}
		return coor_2;
	}
	COORDINATE operator/(const double & divi) const{
		COORDINATE coor_2;
		for(long ti = 0; ti < 3; ti ++){
			coor_2(ti) = (*this)(ti) / divi;
		}
		return coor_2;
	}
	bool COPY(const Vector3d & tempVect){
		for(long ti = 0; ti < 3; ti ++){
			(*this)(ti) = tempVect(ti);
		}
		return true;
	}
	bool COPY_TO(Vector3d & tempVect){
		for(long ti = 0; ti < 3; ti ++){
			tempVect(ti) = (*this)(ti);
		}
		return true;
	}
};

//1d array, two arrays can be compared
class MATRIX:public Matrix<long,Dynamic,Dynamic>{
public:
	bool SORT(){
		MATRIX tempMatr_0 = *this;
		VectorXi inde = VectorXi::LinSpaced(tempMatr_0.rows(), 0, tempMatr_0.rows()-1);
		auto rule = [tempMatr_0](long ti, long tj)->bool{
			return tempMatr_0(ti) < tempMatr_0(tj);
		};
		sort(inde.data(), inde.data()+inde.size(), rule);
		for(long ti = 0; ti < tempMatr_0.rows(); ti ++){
			(*this)(ti) = tempMatr_0(inde(ti));
		}
		return true;
	}
	bool operator<(const MATRIX & matr_2) const{
		MATRIX matr_0 = *this;
		MATRIX matr_1 = matr_2;
		if(matr_0.rows() < matr_1.rows()) return true;
		else if(matr_0.rows() > matr_1.rows()) return false;
		for(long ti=0; ti< matr_0.rows(); ti++){
			matr_0(ti) = abs(matr_0(ti));
			matr_1(ti) = abs(matr_1(ti));
		}
		matr_0.SORT();
		matr_1.SORT();
		for(long ti = 0; ti < matr_0.rows(); ti ++){
			if(matr_0(ti) < matr_1(ti)) return true;
			else if(matr_0(ti) > matr_1(ti)) return false;
		}
		return false;
	}
	bool FLAG_SORT(long flag){
		if(flag < 0){
			MATRIX tempMATR = *this;
			if(this->rows()==2){
				for(long ti = 0; ti < this->rows(); ti ++){
					(*this)(ti) = tempMATR(this->rows() - 1 - ti);
				}
			}
			else{
				for(long ti = 0; ti < this->rows(); ti ++){
					(*this)(ti) = - tempMATR(this->rows() - 1 - ti);
				}
			}
		}
		return true;
	}
};

VEM::VEM(){
	//
	mateElas = 210.0E9;
	matePois = 0.3;
	lineCoun = 0;//for 0: +0=-0
	faceCoun = 0;
	voluCoun = 0;
}

bool VEM::CLEAR(){
	stifList.clear();
	return true;
}

long VEM::TRY_ADD_POINT(COORDINATE tempCoor){
	map<COORDINATE,long>::iterator iter_0 = coorPoin.find(tempCoor);
	if(iter_0 == coorPoin.end()){
		long tempSize = poinCoor.size();
		poinCoor.insert(map<long,COORDINATE>::value_type(tempSize, tempCoor));
		coorPoin.insert(map<COORDINATE,long>::value_type(tempCoor, tempSize));
		return tempSize;
	}
	else{
		return iter_0->second;
	}
}

long VEM::TRY_ADD_LINE(MATRIX tempMATR){
	long flag = 1;
	if(tempMATR(0) > tempMATR(1)){
		flag = -1;
		long temp = tempMATR(0);
		tempMATR(0) = tempMATR(1);
		tempMATR(1) = temp;
	}
	map<MATRIX,long>::iterator iter_0 = poinLine.find(tempMATR);
	if(iter_0 == poinLine.end()){
		lineCoun ++;
		linePoin.insert(map<long,MATRIX>::value_type(lineCoun, tempMATR));
		poinLine.insert(map<MATRIX,long>::value_type(tempMATR, lineCoun));
		return flag * lineCoun;
	}
	else{
		return flag * iter_0->second;
	}
}

long VEM::TRY_ADD_FACE(MATRIX tempMATR){
	map<MATRIX,long>::iterator iter_0 = lineFace.find(tempMATR);
	if(iter_0 == lineFace.end()){
		faceCoun ++;
		faceLine.insert(map<long,MATRIX>::value_type(faceCoun, tempMATR));
		lineFace.insert(map<MATRIX,long>::value_type(tempMATR, faceCoun));
		for(long ti = 0; ti < tempMATR.rows(); ti ++){
			map<long,list<long>>::iterator iter_1 = lineUsed.find(abs(tempMATR(ti)));
			if(iter_1 == lineUsed.end()){
				list<long> tempList;
				tempList.push_back(faceCoun);
				lineUsed.insert(map<long,list<long>>::value_type(abs(tempMATR(ti)), tempList));
			}
			else{
				list<long>::iterator iter_2 = 
					find((iter_1->second).begin(), (iter_1->second).end(), faceCoun);
				if(iter_2 == (iter_1->second).end()){
					(iter_1->second).push_back(faceCoun);
				}
			}
		}
		return faceCoun;//!!!@@@###$$$%%%flag*
	}
	else{
		for(long tj = 0; tj < tempMATR.rows(); tj ++){
			if(abs(tempMATR(tj)) == abs((iter_0->first)(0))){
				if(tempMATR(tj) == (iter_0->first)(0)){
					return iter_0->second;
				}
				else{
					return - iter_0->second;
				}
			}
		}
	}
	cout<<"Error in VEM::TRY_ADD_FACE!"<<endl;
	return -1;
}

long VEM::ADD_VOLUME(MATRIX tempMATR){
	voluCoun ++;
	voluFace.insert(map<long,MATRIX>::value_type(voluCoun, tempMATR));
	for(long ti = 0; ti < tempMATR.rows(); ti ++){
		map<long,list<long>>::iterator iter_0 = faceUsed.find(abs(tempMATR(ti)));
		if(iter_0 == faceUsed.end()){
			list<long> tempList;
			tempList.push_back(voluCoun);
			faceUsed.insert(map<long,list<long>>::value_type(abs(tempMATR(ti)), tempList));
		}
		else{
			list<long>::iterator iter_1 = 
				find((iter_0->second).begin(), (iter_0->second).end(), voluCoun);
			if(iter_1 == (iter_0->second).end()){
				(iter_0->second).push_back(voluCoun);
			}
		}
	}
	return voluCoun;
}

bool VEM::ADD_FACEUSED(long tempFace, long tempVolu){
	map<long,list<long>>::iterator iter_0 = faceUsed.find(abs(tempFace));
	list<long>::iterator iter_1 = 
		find((iter_0->second).begin(), (iter_0->second).end(), tempVolu);
	if(iter_1 == (iter_0->second).end()){
		(iter_0->second).push_back(tempVolu);
	}
	return true;
}

bool VEM::ADD_LINEUSED(long tempLine, long tempFace){
	map<long,list<long>>::iterator iter_0 = lineUsed.find(abs(tempLine));
	list<long>::iterator iter_1 = 
		find((iter_0->second).begin(), (iter_0->second).end(), abs(tempFace));
	if(iter_1 == (iter_0->second).end()){
		(iter_0->second).push_back(abs(tempFace));
	}
	return true;
}

long VEM::ADD_BLOCK_FROM_POINT(const Matrix<long,8,1> & tempCorn){
	Matrix<long,12,1> tempLine;
	MATRIX tempMATR;
	tempMATR.resize(2,1);
	tempMATR << tempCorn(0), tempCorn(1);
	tempLine(0) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(1), tempCorn(2);
	tempLine(1) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(2), tempCorn(3);
	tempLine(2) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(3), tempCorn(0);
	tempLine(3) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(0), tempCorn(4);
	tempLine(4) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(1), tempCorn(5);
	tempLine(5) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(2), tempCorn(6);
	tempLine(6) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(3), tempCorn(7);
	tempLine(7) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(4), tempCorn(5);
	tempLine(8) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(5), tempCorn(6);
	tempLine(9) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(6), tempCorn(7);
	tempLine(10) = TRY_ADD_LINE(tempMATR);
	tempMATR << tempCorn(7), tempCorn(4);
	tempLine(11) = TRY_ADD_LINE(tempMATR);
	MATRIX tempFace;
	tempFace.resize(6,1);
	tempMATR.resize(4,1);
	tempMATR << -tempLine(3), -tempLine(2), -tempLine(1), -tempLine(0);
	tempFace(0) = TRY_ADD_FACE(tempMATR);
	tempMATR << tempLine(8), tempLine(9), tempLine(10), tempLine(11);
	tempFace(1) = TRY_ADD_FACE(tempMATR);
	tempMATR << tempLine(3), tempLine(4), -tempLine(11), -tempLine(7);
	tempFace(2) = TRY_ADD_FACE(tempMATR);
	tempMATR << tempLine(0), tempLine(5), -tempLine(8), -tempLine(4);
	tempFace(3) = TRY_ADD_FACE(tempMATR);
	tempMATR << tempLine(1), tempLine(6), -tempLine(9), -tempLine(5);
	tempFace(4) = TRY_ADD_FACE(tempMATR);
	tempMATR << tempLine(2), tempLine(7), -tempLine(10), -tempLine(6);
	tempFace(5) = TRY_ADD_FACE(tempMATR);
	return ADD_VOLUME(tempFace);
}

bool VEM::DELETE_VOLUME(long tempVolu){
	map<long,MATRIX>::iterator iterVolu = voluFace.find(tempVolu);
	for(long ti = 0; ti < (iterVolu->second).rows(); ti ++){
		long face_i = abs((iterVolu->second)(ti));
		map<long,list<long>>::iterator iterFace = faceUsed.find(face_i);
		(iterFace->second).remove(tempVolu);//faceUsed-2
		if((iterFace->second).size() == 0){
			map<long,MATRIX>::iterator iter_0 = faceLine.find(face_i);
			for(long tj = 0; tj < (iter_0->second).rows(); tj ++){
				long line_j = abs((iter_0->second)(tj));
				map<long,list<long>>::iterator iterLine = lineUsed.find(line_j);
				(iterLine->second).remove(face_i);//lineUsed-2
				if((iterLine->second).size() == 0){
					map<long,MATRIX>::iterator iter_1 = linePoin.find(line_j);
					poinLine.erase(iter_1->second);
					linePoin.erase(line_j);
					lineUsed.erase(line_j);//lineUsed-1
				}
			}
			lineFace.erase(iter_0->second);
			faceLine.erase(face_i);
			faceUsed.erase(face_i);//faceUsed-1
		}
	}
	voluFace.erase(tempVolu);
	return true;
}

bool VEM::REFINE(list<long> & spliVolu_0, const Matrix<long,Dynamic,Dynamic> & spliFlag, 
	const map<MATRIX,COORDINATE> & planSurf){
	list<long> spliVolu_1;
	long tempNicn = 0;
	for(list<long>::iterator iter_0 = spliVolu_0.begin(); 
		iter_0 != spliVolu_0.end(); iter_0 ++){
		if(tempNicn%100 == 0){
			cout<<"The "<<tempNicn<<"/"<<spliVolu_0.size()<<"-th volume."<<endl;
		}
		//*iter_0 is always a hexahedron!!!
		long spliPoin[3][3][3];
		long spliLine[3*12+2*9];
		long spliFace[3*3*4];
		long spliVolu[8];
		//points
		map<long,MATRIX>::iterator iterVolu = voluFace.find(*iter_0);
		long bottFace = (iterVolu->second)(0);
		long uppeFace = (iterVolu->second)(1);
		map<long,MATRIX>::iterator iterFace = faceLine.find(abs(bottFace));
		MATRIX bottLine = iterFace->second;
		bottLine.FLAG_SORT(bottFace);
		iterFace = faceLine.find(abs(uppeFace));
		MATRIX uppeLine = iterFace->second;
		uppeLine.FLAG_SORT(uppeFace);
		map<long,MATRIX>::iterator iterLine = linePoin.find(abs(bottLine(0)));
		MATRIX line_0 = iterLine->second;
		line_0.FLAG_SORT(bottLine(0));
		spliPoin[0][0][0] = line_0(0);
		spliPoin[0][2][0] = line_0(1);
		iterLine = linePoin.find(abs(uppeLine(0)));
		line_0 = iterLine->second;
		line_0.FLAG_SORT(uppeLine(0));
		spliPoin[2][0][0] = line_0(0);
		spliPoin[2][0][2] = line_0(1);
		iterLine = linePoin.find(abs(bottLine(2)));
		line_0 = iterLine->second;
		line_0.FLAG_SORT(bottLine(2));
		spliPoin[0][2][2] = line_0(0);
		spliPoin[0][0][2] = line_0(1);
		iterLine = linePoin.find(abs(uppeLine(2)));
		line_0 = iterLine->second;
		line_0.FLAG_SORT(uppeLine(2));
		spliPoin[2][2][2] = line_0(0);
		spliPoin[2][2][0] = line_0(1);
		Matrix<COORDINATE,8,1> cornCoor;
		cornCoor << poinCoor[spliPoin[0][0][0]],
			poinCoor[spliPoin[0][0][2]],
			poinCoor[spliPoin[0][2][2]],
			poinCoor[spliPoin[0][2][0]],
			poinCoor[spliPoin[2][0][0]],
			poinCoor[spliPoin[2][0][2]],
			poinCoor[spliPoin[2][2][2]],
			poinCoor[spliPoin[2][2][0]];
		MATRIX surfMatr;
		surfMatr.resize(2,1);
		surfMatr << spliPoin[0][0][0], spliPoin[0][0][2];
		map<MATRIX,COORDINATE>::const_iterator iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[0][0][1] = TRY_ADD_POINT((cornCoor(0)+cornCoor(1))/2.0);
		}
		else{
			spliPoin[0][0][1] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[0][0][0], spliPoin[0][2][0];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[0][1][0] = TRY_ADD_POINT((cornCoor(0)+cornCoor(3))/2.0);
		}
		else{
			spliPoin[0][1][0] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(4,1);
		surfMatr << spliPoin[0][0][0], spliPoin[0][0][2], spliPoin[0][2][2], spliPoin[0][2][0];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[0][1][1] = 
				TRY_ADD_POINT((cornCoor(0)+cornCoor(1)+cornCoor(2)+cornCoor(3))/4.0);
		}
		else{
			spliPoin[0][1][1] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[0][0][2], spliPoin[0][2][2];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[0][1][2] = TRY_ADD_POINT((cornCoor(1)+cornCoor(2))/2.0);
		}
		else{
			spliPoin[0][1][2] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[0][2][2], spliPoin[0][2][0];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[0][2][1] = TRY_ADD_POINT((cornCoor(2)+cornCoor(3))/2.0);
		}
		else{
			spliPoin[0][2][1] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[0][0][0], spliPoin[2][0][0];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[1][0][0] = TRY_ADD_POINT((cornCoor(0)+cornCoor(4))/2.0);
		}
		else{
			spliPoin[1][0][0] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(4,1);
		surfMatr << spliPoin[0][0][0], spliPoin[0][0][2], spliPoin[2][0][2], spliPoin[2][0][0];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[1][0][1] = 
				TRY_ADD_POINT((cornCoor(0)+cornCoor(1)+cornCoor(5)+cornCoor(4))/4.0);
		}
		else{
			spliPoin[1][0][1] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[0][0][2], spliPoin[2][0][2];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[1][0][2] = TRY_ADD_POINT((cornCoor(1)+cornCoor(5))/2.0);
		}
		else{
			spliPoin[1][0][2] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(4,1);
		surfMatr << spliPoin[0][0][0], spliPoin[0][2][0], spliPoin[2][2][0], spliPoin[2][0][0];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[1][1][0] = 
				TRY_ADD_POINT((cornCoor(0)+cornCoor(3)+cornCoor(7)+cornCoor(4))/4.0);
		}
		else{
			spliPoin[1][1][0] = TRY_ADD_POINT(iterSurf->second);
		}
		spliPoin[1][1][1] = TRY_ADD_POINT((cornCoor(0)+cornCoor(1)+cornCoor(2)+cornCoor(3)
			+cornCoor(4)+cornCoor(5)+cornCoor(6)+cornCoor(7))/8.0);
		surfMatr.resize(4,1);
		surfMatr << spliPoin[0][0][2], spliPoin[0][2][2], spliPoin[2][2][2], spliPoin[2][0][2];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[1][1][2] = 
				TRY_ADD_POINT((cornCoor(1)+cornCoor(2)+cornCoor(6)+cornCoor(5))/4.0);
		}
		else{
			spliPoin[1][1][2] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[0][2][0], spliPoin[2][2][0];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[1][2][0] = TRY_ADD_POINT((cornCoor(3)+cornCoor(7))/2.0);
		}
		else{
			spliPoin[1][2][0] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(4,1);
		surfMatr << spliPoin[0][2][2], spliPoin[0][2][0], spliPoin[2][2][0], spliPoin[2][2][2];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[1][2][1] = 
				TRY_ADD_POINT((cornCoor(2)+cornCoor(3)+cornCoor(7)+cornCoor(6))/4.0);
		}
		else{
			spliPoin[1][2][1] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[0][2][2], spliPoin[2][2][2];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[1][2][2] = TRY_ADD_POINT((cornCoor(2)+cornCoor(6))/2.0);
		}
		else{
			spliPoin[1][2][2] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[2][0][0], spliPoin[2][0][2];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[2][0][1] = TRY_ADD_POINT((cornCoor(4)+cornCoor(5))/2.0);
		}
		else{
			spliPoin[2][0][1] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[2][0][0], spliPoin[2][2][0];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[2][1][0] = TRY_ADD_POINT((cornCoor(4)+cornCoor(7))/2.0);
		}
		else{
			spliPoin[2][1][0] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(4,1);
		surfMatr << spliPoin[2][0][0], spliPoin[2][0][2], spliPoin[2][2][2], spliPoin[2][2][0];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[2][1][1] = 
				TRY_ADD_POINT((cornCoor(4)+cornCoor(5)+cornCoor(6)+cornCoor(7))/4.0);
		}
		else{
			spliPoin[2][1][1] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[2][0][2], spliPoin[2][2][2];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[2][1][2] = TRY_ADD_POINT((cornCoor(5)+cornCoor(6))/2.0);
		}
		else{
			spliPoin[2][1][2] = TRY_ADD_POINT(iterSurf->second);
		}
		surfMatr.resize(2,1);
		surfMatr << spliPoin[2][2][2], spliPoin[2][2][0];
		iterSurf = planSurf.find(surfMatr);
		if(iterSurf == planSurf.end()){
			spliPoin[2][2][1] = TRY_ADD_POINT((cornCoor(6)+cornCoor(7))/2.0);
		}
		else{
			spliPoin[2][2][1] = TRY_ADD_POINT(iterSurf->second);
		}
		//lines
		MATRIX tempMATR;
		tempMATR.resize(2,1);
		tempMATR << spliPoin[0][0][0], spliPoin[0][0][1];
		spliLine[0] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][0][1], spliPoin[0][0][2];
		spliLine[1] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][0][0], spliPoin[0][1][0];
		spliLine[2] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][0][1], spliPoin[0][1][1];
		spliLine[3] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][0][2], spliPoin[0][1][2];
		spliLine[4] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][1][0], spliPoin[0][1][1];
		spliLine[5] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][1][1], spliPoin[0][1][2];
		spliLine[6] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][1][0], spliPoin[0][2][0];
		spliLine[7] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][1][1], spliPoin[0][2][1];
		spliLine[8] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][1][2], spliPoin[0][2][2];
		spliLine[9] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][2][0], spliPoin[0][2][1];
		spliLine[10] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][2][1], spliPoin[0][2][2];
		spliLine[11] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][0][0], spliPoin[1][0][0];
		spliLine[12] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][0][1], spliPoin[1][0][1];
		spliLine[13] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][0][2], spliPoin[1][0][2];
		spliLine[14] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][1][0], spliPoin[1][1][0];
		spliLine[15] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][1][1], spliPoin[1][1][1];
		spliLine[16] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][1][2], spliPoin[1][1][2];
		spliLine[17] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][2][0], spliPoin[1][2][0];
		spliLine[18] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][2][1], spliPoin[1][2][1];
		spliLine[19] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[0][2][2], spliPoin[1][2][2];
		spliLine[20] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][0][0], spliPoin[1][0][1];
		spliLine[21] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][0][1], spliPoin[1][0][2];
		spliLine[22] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][0][0], spliPoin[1][1][0];
		spliLine[23] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][0][1], spliPoin[1][1][1];
		spliLine[24] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][0][2], spliPoin[1][1][2];
		spliLine[25] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][1][0], spliPoin[1][1][1];
		spliLine[26] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][1][1], spliPoin[1][1][2];
		spliLine[27] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][1][0], spliPoin[1][2][0];
		spliLine[28] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][1][1], spliPoin[1][2][1];
		spliLine[29] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][1][2], spliPoin[1][2][2];
		spliLine[30] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][2][0], spliPoin[1][2][1];
		spliLine[31] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][2][1], spliPoin[1][2][2];
		spliLine[32] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][0][0], spliPoin[2][0][0];
		spliLine[33] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][0][1], spliPoin[2][0][1];
		spliLine[34] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][0][2], spliPoin[2][0][2];
		spliLine[35] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][1][0], spliPoin[2][1][0];
		spliLine[36] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][1][1], spliPoin[2][1][1];
		spliLine[37] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][1][2], spliPoin[2][1][2];
		spliLine[38] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][2][0], spliPoin[2][2][0];
		spliLine[39] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][2][1], spliPoin[2][2][1];
		spliLine[40] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[1][2][2], spliPoin[2][2][2];
		spliLine[41] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][0][0], spliPoin[2][0][1];
		spliLine[42] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][0][1], spliPoin[2][0][2];
		spliLine[43] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][0][0], spliPoin[2][1][0];
		spliLine[44] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][0][1], spliPoin[2][1][1];
		spliLine[45] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][0][2], spliPoin[2][1][2];
		spliLine[46] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][1][0], spliPoin[2][1][1];
		spliLine[47] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][1][1], spliPoin[2][1][2];
		spliLine[48] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][1][0], spliPoin[2][2][0];
		spliLine[49] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][1][1], spliPoin[2][2][1];
		spliLine[50] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][1][2], spliPoin[2][2][2];
		spliLine[51] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][2][0], spliPoin[2][2][1];
		spliLine[52] = TRY_ADD_LINE(tempMATR);
		tempMATR << spliPoin[2][2][1], spliPoin[2][2][2];
		spliLine[53] = TRY_ADD_LINE(tempMATR);
		//faces
		tempMATR.resize(4,1);
		tempMATR << spliLine[0], spliLine[3], -spliLine[5], -spliLine[2];
		spliFace[0] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[1], spliLine[4], -spliLine[6], -spliLine[3];
		spliFace[1] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[5], spliLine[8], -spliLine[10], -spliLine[7];
		spliFace[2] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[6], spliLine[9], -spliLine[11], -spliLine[8];
		spliFace[3] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[21], spliLine[24], -spliLine[26], -spliLine[23];
		spliFace[4] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[22], spliLine[25], -spliLine[27], -spliLine[24];
		spliFace[5] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[26], spliLine[29], -spliLine[31], -spliLine[28];
		spliFace[6] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[27], spliLine[30], -spliLine[32], -spliLine[29];
		spliFace[7] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[42], spliLine[45], -spliLine[47], -spliLine[44];
		spliFace[8] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[43], spliLine[46], -spliLine[48], -spliLine[45];
		spliFace[9] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[47], spliLine[50], -spliLine[52], -spliLine[49];
		spliFace[10] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[48], spliLine[51], -spliLine[53], -spliLine[50];
		spliFace[11] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[2], spliLine[15], -spliLine[23], -spliLine[12];
		spliFace[12] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[7], spliLine[18], -spliLine[28], -spliLine[15];
		spliFace[13] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[23], spliLine[36], -spliLine[44], -spliLine[33];
		spliFace[14] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[28], spliLine[39], -spliLine[49], -spliLine[36];
		spliFace[15] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[3], spliLine[16], -spliLine[24], -spliLine[13];
		spliFace[16] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[8], spliLine[19], -spliLine[29], -spliLine[16];
		spliFace[17] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[24], spliLine[37], -spliLine[45], -spliLine[34];
		spliFace[18] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[29], spliLine[40], -spliLine[50], -spliLine[37];
		spliFace[19] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[4], spliLine[17], -spliLine[25], -spliLine[14];
		spliFace[20] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[9], spliLine[20], -spliLine[30], -spliLine[17];
		spliFace[21] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[25], spliLine[38], -spliLine[46], -spliLine[35];
		spliFace[22] = TRY_ADD_FACE(tempMATR);
		tempMATR << spliLine[30], spliLine[41], -spliLine[51], -spliLine[38];
		spliFace[23] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[0], spliLine[12], spliLine[21], -spliLine[13];
		spliFace[24] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[1], spliLine[13], spliLine[22], -spliLine[14];
		spliFace[25] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[21], spliLine[33], spliLine[42], -spliLine[34];
		spliFace[26] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[22], spliLine[34], spliLine[43], -spliLine[35];
		spliFace[27] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[5], spliLine[15], spliLine[26], -spliLine[16];
		spliFace[28] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[6], spliLine[16], spliLine[27], -spliLine[17];
		spliFace[29] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[26], spliLine[36], spliLine[47], -spliLine[37];
		spliFace[30] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[27], spliLine[37], spliLine[48], -spliLine[38];
		spliFace[31] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[10], spliLine[18], spliLine[31], -spliLine[19];
		spliFace[32] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[11], spliLine[19], spliLine[32], -spliLine[20];
		spliFace[33] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[31], spliLine[39], spliLine[52], -spliLine[40];
		spliFace[34] = TRY_ADD_FACE(tempMATR);
		tempMATR << -spliLine[32], spliLine[40], spliLine[53], -spliLine[41];
		spliFace[35] = TRY_ADD_FACE(tempMATR);
		//volumes
		tempMATR.resize(6,1);
		tempMATR << 
			-spliFace[0], spliFace[4], -spliFace[12], -spliFace[24], spliFace[16], spliFace[28];
		spliVolu[0] = ADD_VOLUME(tempMATR);
		tempMATR << 
			-spliFace[1], spliFace[5], -spliFace[16], -spliFace[25], spliFace[20], spliFace[29];
		spliVolu[1] = ADD_VOLUME(tempMATR);
		tempMATR << 
			-spliFace[2], spliFace[6], -spliFace[13], -spliFace[28], spliFace[17], spliFace[32];
		spliVolu[2] = ADD_VOLUME(tempMATR);
		tempMATR << 
			-spliFace[3], spliFace[7], -spliFace[17], -spliFace[29], spliFace[21], spliFace[33];
		spliVolu[3] = ADD_VOLUME(tempMATR);
		tempMATR << 
			-spliFace[4], spliFace[8], -spliFace[14], -spliFace[26], spliFace[18], spliFace[30];
		spliVolu[4] = ADD_VOLUME(tempMATR);
		tempMATR << 
			-spliFace[5], spliFace[9], -spliFace[18], -spliFace[27], spliFace[22], spliFace[31];
		spliVolu[5] = ADD_VOLUME(tempMATR);
		tempMATR << 
			-spliFace[6], spliFace[10], -spliFace[15], -spliFace[30], spliFace[19], spliFace[34];
		spliVolu[6] = ADD_VOLUME(tempMATR);
		tempMATR << 
			-spliFace[7], spliFace[11], -spliFace[19], -spliFace[31], spliFace[23], spliFace[35];
		spliVolu[7] = ADD_VOLUME(tempMATR);
		if(spliFlag(tempNicn) ==0){
			spliVolu_1.push_back(spliVolu[0]);
			spliVolu_1.push_back(spliVolu[1]);
			spliVolu_1.push_back(spliVolu[2]);
			spliVolu_1.push_back(spliVolu[3]);
		}
		else if(spliFlag(tempNicn) == 1){
			spliVolu_1.push_back(spliVolu[4]);
			spliVolu_1.push_back(spliVolu[5]);
			spliVolu_1.push_back(spliVolu[6]);
			spliVolu_1.push_back(spliVolu[7]);
		}
		else if(spliFlag(tempNicn) == 2){
			spliVolu_1.push_back(spliVolu[2]);
			spliVolu_1.push_back(spliVolu[0]);
			spliVolu_1.push_back(spliVolu[6]);
			spliVolu_1.push_back(spliVolu[4]);
		}
		else if(spliFlag(tempNicn) == 3){
			spliVolu_1.push_back(spliVolu[0]);
			spliVolu_1.push_back(spliVolu[1]);
			spliVolu_1.push_back(spliVolu[4]);
			spliVolu_1.push_back(spliVolu[5]);
		}
		else if(spliFlag(tempNicn) == 4){
			spliVolu_1.push_back(spliVolu[1]);
			spliVolu_1.push_back(spliVolu[3]);
			spliVolu_1.push_back(spliVolu[5]);
			spliVolu_1.push_back(spliVolu[7]);
		}
		else if(spliFlag(tempNicn) == 5){
			spliVolu_1.push_back(spliVolu[3]);
			spliVolu_1.push_back(spliVolu[2]);
			spliVolu_1.push_back(spliVolu[7]);
			spliVolu_1.push_back(spliVolu[6]);
		}
		else if(spliFlag(tempNicn) == 6){
			spliVolu_1.push_back(spliVolu[0]);
			spliVolu_1.push_back(spliVolu[1]);
			spliVolu_1.push_back(spliVolu[2]);
			spliVolu_1.push_back(spliVolu[3]);
			spliVolu_1.push_back(spliVolu[4]);
			spliVolu_1.push_back(spliVolu[5]);
			spliVolu_1.push_back(spliVolu[6]);
			spliVolu_1.push_back(spliVolu[7]);
		}
		else if(spliFlag(tempNicn) == 7){
			spliVolu_1.push_back(spliVolu[1]);
			spliVolu_1.push_back(spliVolu[5]);
		}
		else if(spliFlag(tempNicn) == 8){
			spliVolu_1.push_back(spliVolu[0]);
			spliVolu_1.push_back(spliVolu[4]);
		}
		//
		for(long ti = 0; ti < (iterVolu->second).rows(); ti ++){
			long face_i = abs((iterVolu->second)(ti));
			map<long,list<long>>::iterator iter_i = faceUsed.find(face_i);
			list<long> tempList = iter_i->second;
			for(list<long>::iterator iter_1 = tempList.begin(); 
				iter_1 != tempList.end(); iter_1 ++){
				if(find(spliVolu_0.begin(), spliVolu_0.end(), *iter_1) == spliVolu_0.end()){
					map<long,MATRIX>::iterator iterVolu_1 = voluFace.find(*iter_1);
					tempMATR.resize((iterVolu_1->second).rows() + 3,1);
					long tempSize = -1;
					for(long tj = 0; tj < (iterVolu_1->second).rows(); tj ++){
						if(face_i == abs((iterVolu_1->second)(tj))){
							if(ti == 0){
								tempSize ++;
								tempMATR(tempSize) = spliFace[0];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = spliFace[1];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = spliFace[2];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = spliFace[3];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
							}
							else if(ti == 1){
								tempSize ++;
								tempMATR(tempSize) = -spliFace[8];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = -spliFace[9];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = -spliFace[10];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = -spliFace[11];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
							}
							else if(ti == 2){
								tempSize ++;
								tempMATR(tempSize) = spliFace[12];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = spliFace[13];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = spliFace[14];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = spliFace[15];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
							}
							else if(ti == 3){
								tempSize ++;
								tempMATR(tempSize) = spliFace[24];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = spliFace[25];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = spliFace[26];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = spliFace[27];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
							}
							else if(ti == 4){
								tempSize ++;
								tempMATR(tempSize) = -spliFace[20];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = -spliFace[21];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = -spliFace[22];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = -spliFace[23];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
							}
							else if(ti == 5){
								tempSize ++;
								tempMATR(tempSize) = -spliFace[32];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = -spliFace[33];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = -spliFace[34];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
								tempSize ++;
								tempMATR(tempSize) = -spliFace[35];
								ADD_FACEUSED(abs(tempMATR(tempSize)), *iter_1);
							}
						}
						else{
							tempSize ++;
							tempMATR(tempSize) = (iterVolu_1->second)(tj);
						}
					}
					iterVolu_1->second = tempMATR;//voluFace
					//faceUsed: at this time, (iter_i->second) will not become empty
					(iter_i->second).remove(*iter_1);
				}
			}
		}
		//
		Matrix<MATRIX,12,1> origLine_0;
		for(long ti = 0; ti < 12; ti ++){
			origLine_0(ti).resize(2,1);
		}
		origLine_0(0) << spliPoin[0][0][0], spliPoin[0][0][2];
		origLine_0(1) << spliPoin[0][0][2], spliPoin[0][2][2];
		origLine_0(2) << spliPoin[0][2][0], spliPoin[0][2][2];
		origLine_0(3) << spliPoin[0][0][0], spliPoin[0][2][0];
		origLine_0(4) << spliPoin[0][0][0], spliPoin[2][0][0];
		origLine_0(5) << spliPoin[0][0][2], spliPoin[2][0][2];
		origLine_0(6) << spliPoin[0][2][2], spliPoin[2][2][2];
		origLine_0(7) << spliPoin[0][2][0], spliPoin[2][2][0];
		origLine_0(8) << spliPoin[2][0][0], spliPoin[2][0][2];
		origLine_0(9) << spliPoin[2][0][2], spliPoin[2][2][2];
		origLine_0(10) << spliPoin[2][2][0], spliPoin[2][2][2];
		origLine_0(11) << spliPoin[2][0][0], spliPoin[2][2][0];
		Matrix<long,12,1> origLine;
		origLine(0) = poinLine[origLine_0(0)];
		origLine(1) = poinLine[origLine_0(1)];
		origLine(2) = poinLine[origLine_0(2)];
		origLine(3) = poinLine[origLine_0(3)];
		origLine(4) = poinLine[origLine_0(4)];
		origLine(5) = poinLine[origLine_0(5)];
		origLine(6) = poinLine[origLine_0(6)];
		origLine(7) = poinLine[origLine_0(7)];
		origLine(8) = poinLine[origLine_0(8)];
		origLine(9) = poinLine[origLine_0(9)];
		origLine(10) = poinLine[origLine_0(10)];
		origLine(11) = poinLine[origLine_0(11)];
		for(long ti = 0; ti < 12; ti ++){
			map<long,list<long>>::iterator iter_2 = lineUsed.find(origLine(ti));
			list<long> tempList = iter_2->second;
			for(list<long>::iterator iter_3 = tempList.begin(); 
				iter_3 != tempList.end(); iter_3 ++){
				bool isin = false;
				map<long,list<long>>::iterator iter_4 = faceUsed.find(abs(*iter_3));
				for(list<long>::iterator iter_5 = (iter_4->second).begin(); 
					iter_5 != (iter_4->second).end(); iter_5 ++){
					if(find(spliVolu_0.begin(),spliVolu_0.end(),*iter_5) != spliVolu_0.end()){
						//in E
						isin = true;
					}
				}
				if(!isin){//update face *iter_3
					map<long,MATRIX>::iterator iterFace = faceLine.find(abs(*iter_3));
					tempMATR.resize((iterFace->second).rows() + 1,1);
					long tempSize = -1;
					for(long tj = 0; tj < (iterFace->second).rows(); tj ++){
						if(abs((iterFace->second)(tj)) == origLine(ti)){
							long tempFlag_0 = 1;
							if(origLine_0(ti)(0) > origLine_0(ti)(1)){
								tempFlag_0 = -1;
							}
							long spliLine_i[2];
							if(ti==0){
								spliLine_i[0] = spliLine[0];
								spliLine_i[1] = spliLine[1];
							}
							else if(ti==1){
								spliLine_i[0] = spliLine[4];
								spliLine_i[1] = spliLine[9];
							}
							else if(ti==2){
								spliLine_i[0] = spliLine[10];
								spliLine_i[1] = spliLine[11];
							}
							else if(ti==3){
								spliLine_i[0] = spliLine[2];
								spliLine_i[1] = spliLine[7];
							}
							else if(ti==4){
								spliLine_i[0] = spliLine[12];
								spliLine_i[1] = spliLine[33];
							}
							else if(ti==5){
								spliLine_i[0] = spliLine[14];
								spliLine_i[1] = spliLine[35];
							}
							else if(ti==6){
								spliLine_i[0] = spliLine[20];
								spliLine_i[1] = spliLine[41];
							}
							else if(ti==7){
								spliLine_i[0] = spliLine[18];
								spliLine_i[1] = spliLine[39];
							}
							else if(ti==8){
								spliLine_i[0] = spliLine[42];
								spliLine_i[1] = spliLine[43];
							}
							else if(ti==9){
								spliLine_i[0] = spliLine[46];
								spliLine_i[1] = spliLine[51];
							}
							else if(ti==10){
								spliLine_i[0] = spliLine[52];
								spliLine_i[1] = spliLine[53];
							}
							else if(ti==11){
								spliLine_i[0] = spliLine[44];
								spliLine_i[1] = spliLine[49];
							}
							if((iterFace->second)(tj) == tempFlag_0*origLine(ti)){
								tempSize ++;
								tempMATR(tempSize) = spliLine_i[0];
								ADD_LINEUSED(abs(tempMATR(tempSize)), abs(*iter_3));
								tempSize ++;
								tempMATR(tempSize) = spliLine_i[1];
								ADD_LINEUSED(abs(tempMATR(tempSize)), abs(*iter_3));
							}
							else{
								tempSize ++;
								tempMATR(tempSize) = -spliLine_i[1];
								ADD_LINEUSED(abs(tempMATR(tempSize)), abs(*iter_3));
								tempSize ++;
								tempMATR(tempSize) = -spliLine_i[0];
								ADD_LINEUSED(abs(tempMATR(tempSize)), abs(*iter_3));
							}
						}
						else{
							tempSize ++;
							tempMATR(tempSize) = (iterFace->second)(tj);
						}
					}
					lineFace.erase(iterFace->second);
					iterFace->second = tempMATR;//faceLine
					lineFace.insert(map<MATRIX,long>::value_type(tempMATR, abs(*iter_3)));
					//lineUsed: at this time, (iter_2->second) will not become empty
					(iter_2->second).remove(abs(*iter_3));
				}
			}
		}
		DELETE_VOLUME(*iter_0);
		tempNicn ++;
	}
	spliVolu_0 = spliVolu_1;
	return true;
}

bool VEM::OUTPUT_VOLUME(long fileIden){
	stringstream tempStre;
	tempStre << "resuPoin_" << fileIden << ".txt";
	ofstream tempOfst(DIRECTORY(tempStre.str()), ios::out);
	tempStre.str("");
	tempStre.clear();
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	for(map<long,COORDINATE>::iterator iter_0 = poinCoor.begin(); 
		iter_0 != poinCoor.end(); iter_0 ++){
		tempOfst<<setw(30)<<(iter_0->second)(0)
			<<setw(30)<<(iter_0->second)(1)
			<<setw(30)<<(iter_0->second)(2)<<endl;
	}
	tempOfst.close();
	tempStre << "resuFace_" << fileIden << ".txt";
	tempOfst.open(DIRECTORY(tempStre.str()), ios::out);
	tempStre.str("");
	tempStre.clear();
	long maxiPoin = -1;
	for(map<long,MATRIX>::iterator iter_0 = voluFace.begin(); 
		iter_0 != voluFace.end(); iter_0 ++){
		for(long ti = 0; ti < (iter_0->second).rows(); ti ++){
			long face_i = abs((iter_0->second)(ti));
			map<long,MATRIX>::iterator iter_1 = faceLine.find(face_i);
			long tempNumb = (iter_1->second).rows();
			maxiPoin = max(maxiPoin, tempNumb);
		}
	}
	for(map<long,MATRIX>::iterator iter_0 = voluFace.begin(); 
		iter_0 != voluFace.end(); iter_0 ++){
		for(long ti = 0; ti < (iter_0->second).rows(); ti ++){
			long face_i = abs((iter_0->second)(ti));
			map<long,MATRIX>::iterator iter_1 = faceLine.find(face_i);
			long lastIden = 0;
			for(long tj = 0; tj < (iter_1->second).rows(); tj ++){
				long line_j = (iter_1->second)(tj);
				map<long,MATRIX>::iterator iter_2 = linePoin.find(abs(line_j));
				if(line_j > 0){
					tempOfst<<setw(10)<<(iter_2->second)(0);
					lastIden = (iter_2->second)(0);
				}
				else{
					tempOfst<<setw(10)<<(iter_2->second)(1);
					lastIden = (iter_2->second)(1);
				}
			}
			if((iter_1->second).rows() < maxiPoin){
				for(long tj = (iter_1->second).rows(); tj < maxiPoin; tj ++){
					tempOfst<<setw(10)<<lastIden;
				}
			}
			tempOfst<<endl;
		}
	}
	tempOfst.close();
	tempStre << "resuVolu_" << fileIden << ".txt";
	tempOfst.open(DIRECTORY(tempStre.str()), ios::out);
	tempStre.str("");
	tempStre.clear();
	for(map<long,MATRIX>::iterator iter_0 = voluFace.begin(); 
		iter_0 != voluFace.end(); iter_0 ++){
		for(long ti = 0; ti < (iter_0->second).rows(); ti ++){
			tempOfst << iter_0->first << endl;
		}
	}
	tempOfst.close();
	return true;
}

bool VEM::VOLUME_2_ELEMENT(){
	//
	long nodeNumb = poinCoor.size();
	long elemNumb = 0;
	elemFace.resize(voluFace.size(), 1);
	map<COORDINATE,long> newNode;
	for(map<long,MATRIX>::iterator iter_0 = voluFace.begin(); 
		iter_0 != voluFace.end(); iter_0 ++){
		long faceNumb_ = 0;//number of faces
		for(long ti = 0; ti < (iter_0->second).rows(); ti ++){
			long face_i = (iter_0->second)(ti);
			map<long,MATRIX>::iterator iter_1 = faceLine.find(abs(face_i));
			if((iter_1->second).rows() == 4){
				faceNumb_ += 2;
			}
			else if((iter_1->second).rows() > 4){
				faceNumb_ += (iter_1->second).rows();
			}
			else{
				cout<<"Error in VEM::VOLUME_2_ELEMENT!"<<endl;
			}
		}
		elemFace(elemNumb).resize(faceNumb_,3);
		faceNumb_ = 0;
		for(long ti = 0; ti < (iter_0->second).rows(); ti ++){
			long face_i = (iter_0->second)(ti);
			map<long,MATRIX>::iterator iter_1 = faceLine.find(abs(face_i));
			MATRIX tempMATR;//points of the face_i
			tempMATR.resize((iter_1->second).rows(),1);
			if(face_i > 0){
				for(long tj = 0; tj < (iter_1->second).rows(); tj ++){
					long line_j = (iter_1->second)(tj);
					map<long,MATRIX>::iterator iter_2 = linePoin.find(abs(line_j));
					if(line_j > 0){
						tempMATR(tj) = (iter_2->second)(0);
					}
					else{
						tempMATR(tj) = (iter_2->second)(1);
					}
				}
			}
			else{
				for(long tj = (iter_1->second).rows() - 1; tj >= 0; tj --){
					long line_j = (iter_1->second)(tj);
					map<long,MATRIX>::iterator iter_2 = linePoin.find(abs(line_j));
					if(line_j > 0){
						tempMATR((iter_1->second).rows()-1-tj) = (iter_2->second)(1);
					}
					else{
						tempMATR((iter_1->second).rows()-1-tj) = (iter_2->second)(0);
					}
				}
			}
			if(tempMATR.rows() == 4){
				elemFace(elemNumb).block(faceNumb_,0,1,3) 
					<< tempMATR(0), tempMATR(1), tempMATR(2);
				faceNumb_ ++;
				elemFace(elemNumb).block(faceNumb_,0,1,3) 
					<< tempMATR(0), tempMATR(2), tempMATR(3);
				faceNumb_ ++;
			}
			else if(tempMATR.rows() > 4){
				COORDINATE tempCoor;//try to add center point (not strictly)
				tempCoor(0) = 0.0;
				tempCoor(1) = 0.0;
				tempCoor(2) = 0.0;
				for(long tj = 0; tj < tempMATR.rows(); tj ++){
					map<long,COORDINATE>::iterator iter_2 = poinCoor.find(tempMATR(tj));
					tempCoor = tempCoor + iter_2->second;
				}
				tempCoor = tempCoor / tempMATR.rows();
				long centIden;
				map<COORDINATE,long>::iterator iter_3 = newNode.find(tempCoor);
				if(iter_3 == newNode.end()){
					centIden = nodeNumb;
					newNode.insert(map<COORDINATE,long>::value_type(tempCoor, nodeNumb));
					nodeNumb ++;
				}
				else{
					centIden = iter_3->second;
				}
				for(long tj = 0; tj < tempMATR.rows(); tj ++){
					if(tj < tempMATR.rows()-1){
						elemFace(elemNumb).block(faceNumb_,0,1,3) 
							<< centIden, tempMATR(tj), tempMATR(tj+1);
					}
					else{
						elemFace(elemNumb).block(faceNumb_,0,1,3) 
							<< centIden, tempMATR(tj), tempMATR(0);
					}
					faceNumb_ ++;
				}
			}
		}
		elemNumb ++;
	}
	nodeCoor.resize(nodeNumb, 3);
	for(map<long,COORDINATE>::iterator iter_0 = poinCoor.begin(); 
		iter_0 != poinCoor.end(); iter_0 ++){
		nodeCoor.block(iter_0->first,0,1,3) 
			<< (iter_0->second)(0), (iter_0->second)(1), (iter_0->second)(2);
	}
	for(map<COORDINATE,long>::iterator iter_0 = newNode.begin(); 
		iter_0 != newNode.end(); iter_0 ++){
		nodeCoor.block(iter_0->second,0,1,3) 
			<< (iter_0->first)(0), (iter_0->first)(1), (iter_0->first)(2);
	}
	//release
	poinCoor.clear();
	coorPoin.clear();
	linePoin.clear();
	poinLine.clear();
	faceLine.clear();
	lineFace.clear();
	lineUsed.clear();
	voluFace.clear();
	faceUsed.clear();
	lineCoun = 0;//for 0: +0=-0
	faceCoun = 0;
	voluCoun = 0;
	return true;
}

bool VEM::OUTPUT_ELEMENT(long fileIden){
	stringstream tempStre;
	tempStre << "resuNode_" << fileIden << ".txt";
	ofstream tempOfst(DIRECTORY(tempStre.str()), ios::out);
	tempStre.str("");
	tempStre.clear();
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	for(long ti = 0; ti < nodeCoor.rows(); ti ++){
		tempOfst<<setw(30)<<nodeCoor(ti,0)
			<<setw(30)<<nodeCoor(ti,1)
			<<setw(30)<<nodeCoor(ti,2)<<endl;
	}
	tempOfst.close();
	//
	tempStre << "resuElem_" << fileIden << ".txt";
	tempOfst.open(DIRECTORY(tempStre.str()), ios::out);
	tempStre.str("");
	tempStre.clear();
	for(long ti = 0; ti < elemFace.rows(); ti ++){
		for(long tj = 0; tj < elemFace(ti).rows(); tj ++){
			tempOfst<<setw(10)<<elemFace(ti)(tj,0)
				<<setw(10)<<elemFace(ti)(tj,1)
				<<setw(10)<<elemFace(ti)(tj,2)
				<<setw(10)<<ti<<endl;
		}
	}
	tempOfst.close();
	return true;
}

bool VEM::ELEMENT_STIFF(){
	OUTPUT_TIME("Element stiffness");
	double mateLame[2];
	double mateTrac;
	mateLame[0] = mateElas * matePois / (1.0 + matePois) / (1.0 - 2.0 * matePois);
	mateLame[1] = mateElas / 2.0 / (1.0 + matePois);
	mateTrac = ((1.0 - matePois) + 2.0 * (1.0 - 2.0 * matePois))
		* 3.0 * mateElas / (1.0 + matePois) / (1.0 - 2.0 * matePois);
	Matrix<Matrix<double,3,3>,12,1> epslBasi, sigmBasi;
	Matrix<double,12,12> matrG;
	for(long ti = 0; ti< 12; ti ++){
		epslBasi(ti) = MatrixXd::Zero(3, 3);
	}
	epslBasi(3)(0,0) = 1.0;
	epslBasi(4)(0,1) = 0.5;
	epslBasi(4)(1,0) = 0.5;
	epslBasi(5)(0,2) = 0.5;
	epslBasi(5)(2,0) = 0.5;
	epslBasi(6)(0,1) = 0.5;
	epslBasi(6)(1,0) = 0.5;
	epslBasi(7)(1,1) = 1.0;
	epslBasi(8)(1,2) = 0.5;
	epslBasi(8)(2,1) = 0.5;
	epslBasi(9)(0,2) = 0.5;
	epslBasi(9)(2,0) = 0.5;
	epslBasi(10)(1,2) = 0.5;
	epslBasi(10)(2,1) = 0.5;
	epslBasi(11)(2,2) = 1.0;
	for(long ti = 0; ti < 12; ti ++){
		sigmBasi(ti) = mateLame[0] * epslBasi(ti).trace() * MatrixXd::Identity(3,3)
			+ 2.0 * mateLame[1] * epslBasi(ti);
		for(long tj = 0; tj < 12; tj ++){
			matrG(ti,tj) = sigmBasi(ti).cwiseProduct(epslBasi(tj)).sum();
		}
	}
	for(long ti = 0; ti < elemFace.rows(); ti ++){
		double elemVolu = 0.0;
		map<long,long> nodeSequ;
		Matrix<Vector3d,Dynamic,Dynamic> faceNorm;
		faceNorm.resize(elemFace(ti).rows(), 1);
		Matrix<double,Dynamic,Dynamic> faceArea;
		faceArea.resize(elemFace(ti).rows(), 1);
		for(long tj = 0; tj < elemFace(ti).rows(); tj ++){
			Matrix<Vector3d,3,1> triaVert;
			for(long tk = 0; tk < 3; tk ++){
				triaVert(tk) = nodeCoor.block(elemFace(ti)(tj,tk),0,1,3).transpose();
				map<long,long>::iterator iter_0 = nodeSequ.find(elemFace(ti)(tj,tk));
				if(iter_0 == nodeSequ.end()){
					long tempSize = nodeSequ.size();
					nodeSequ.insert(map<long,long>::value_type(elemFace(ti)(tj,tk), tempSize));
				}
			}
			Vector3d n_F = (triaVert(1) - triaVert(0)).cross(triaVert(2) - triaVert(0));
			elemVolu += triaVert(0).dot(n_F);
			faceNorm(tj) = n_F.normalized();
			faceArea(tj) = n_F.norm() / 2.0;
		}
		elemVolu /= 6.0;
		Vector3d elemAver = MatrixXd::Zero(3,1);
		for(map<long,long>::iterator iter_0 = nodeSequ.begin(); 
			iter_0 != nodeSequ.end(); iter_0 ++){
			elemAver = elemAver + nodeCoor.block(iter_0->first,0,1,3).transpose();
		}
		elemAver = elemAver / nodeSequ.size();
		Matrix<Matrix<double,3,3>,Dynamic,Dynamic> gradBasi;
		gradBasi.resize(3*nodeSequ.size(), 1);
		for(long tj = 0; tj < elemFace(ti).rows(); tj ++){
			for(long tk = 0; tk < 3; tk ++){
				for(long tm = 0; tm < 3; tm ++){
					long idkm = 3 * nodeSequ[elemFace(ti)(tj,tk)] + tm;
					gradBasi(idkm) = MatrixXd::Zero(3,3);
				}
			}
		}
		for(long tj = 0; tj < elemFace(ti).rows(); tj ++){
			for(long tk = 0; tk < 3; tk ++){
				for(long tm = 0; tm < 3; tm ++){
					long idkm = 3 * nodeSequ[elemFace(ti)(tj,tk)] + tm;
					gradBasi(idkm).block(tm,0,1,3) = gradBasi(idkm).block(tm,0,1,3).eval() 
						+ faceArea(tj) / 3.0 * faceNorm(tj).transpose() / elemVolu;
				}
			}
		}
		Matrix<double,Dynamic,Dynamic> matrP = MatrixXd::Zero(12, 3 * nodeSequ.size());
		Matrix<double,Dynamic,Dynamic> matrD = MatrixXd::Zero(3 * nodeSequ.size(), 12);
		Matrix<double,Dynamic,Dynamic> matrN_C = MatrixXd::Zero(3 * nodeSequ.size(), 6);
		for(map<long,long>::iterator iter_0 = nodeSequ.begin(); 
			iter_0 != nodeSequ.end(); iter_0 ++){
			matrD(3 * (iter_0->second) + 0, 0) = 1.0;
			matrD(3 * (iter_0->second) + 1, 1) = 1.0;
			matrD(3 * (iter_0->second) + 2, 2) = 1.0;
			matrD(3 * (iter_0->second) + 0, 3) = nodeCoor(iter_0->first, 0) - elemAver(0);
			matrD(3 * (iter_0->second) + 1, 4) = nodeCoor(iter_0->first, 0) - elemAver(0);
			matrD(3 * (iter_0->second) + 2, 5) = nodeCoor(iter_0->first, 0) - elemAver(0);
			matrD(3 * (iter_0->second) + 0, 6) = nodeCoor(iter_0->first, 1) - elemAver(1);
			matrD(3 * (iter_0->second) + 1, 7) = nodeCoor(iter_0->first, 1) - elemAver(1);
			matrD(3 * (iter_0->second) + 2, 8) = nodeCoor(iter_0->first, 1) - elemAver(1);
			matrD(3 * (iter_0->second) + 0, 9) = nodeCoor(iter_0->first, 2) - elemAver(2);
			matrD(3 * (iter_0->second) + 1, 10) = nodeCoor(iter_0->first, 2) - elemAver(2);
			matrD(3 * (iter_0->second) + 2, 11) = nodeCoor(iter_0->first, 2) - elemAver(2);
			matrN_C(3 * (iter_0->second) + 0, 0) = nodeCoor(iter_0->first, 0) - elemAver(0);
			matrN_C(3 * (iter_0->second) + 1, 1) = nodeCoor(iter_0->first, 1) - elemAver(1);
			matrN_C(3 * (iter_0->second) + 2, 2) = nodeCoor(iter_0->first, 2) - elemAver(2);
			matrN_C(3 * (iter_0->second) + 0, 3) = nodeCoor(iter_0->first, 1) - elemAver(1);
			matrN_C(3 * (iter_0->second) + 1, 3) = nodeCoor(iter_0->first, 0) - elemAver(0);
			matrN_C(3 * (iter_0->second) + 1, 4) = nodeCoor(iter_0->first, 2) - elemAver(2);
			matrN_C(3 * (iter_0->second) + 2, 4) = nodeCoor(iter_0->first, 1) - elemAver(1);
			matrN_C(3 * (iter_0->second) + 0, 5) = nodeCoor(iter_0->first, 2) - elemAver(2);
			matrN_C(3 * (iter_0->second) + 2, 5) = nodeCoor(iter_0->first, 0) - elemAver(0);
		}
		for(map<long,long>::iterator iter_0 = nodeSequ.begin(); 
			iter_0 != nodeSequ.end(); iter_0 ++){
			for(long tj = 0; tj < 3; tj ++){
				matrP(tj, 3*(iter_0->second) + tj) = 1.0/ nodeSequ.size();
				matrP(3, 3 * (iter_0->second) + tj) = gradBasi(3 * (iter_0->second) + tj)(0,0);
				matrP(6, 3 * (iter_0->second) + tj) = gradBasi(3 * (iter_0->second) + tj)(0,1);
				matrP(9, 3 * (iter_0->second) + tj) = gradBasi(3 * (iter_0->second) + tj)(0,2);
				matrP(4, 3 * (iter_0->second) + tj) = gradBasi(3 * (iter_0->second) + tj)(1,0);
				matrP(7, 3 * (iter_0->second) + tj) = gradBasi(3 * (iter_0->second) + tj)(1,1);
				matrP(10, 3 * (iter_0->second) + tj) = gradBasi(3 * (iter_0->second) + tj)(1,2);
				matrP(5, 3 * (iter_0->second) + tj) = gradBasi(3 * (iter_0->second) + tj)(2,0);
				matrP(8, 3 * (iter_0->second) + tj) = gradBasi(3 * (iter_0->second) + tj)(2,1);
				matrP(11, 3 * (iter_0->second) + tj) = gradBasi(3 * (iter_0->second) + tj)(2,2);
			}
		}
		Matrix<double,Dynamic,Dynamic> matrS = matrD * matrP;
		double alph_E = 1.0 * elemVolu * mateTrac / ((matrN_C.transpose() * matrN_C).trace());
		Matrix<double,Dynamic,Dynamic> matrH = 
			MatrixXd::Identity(3 * nodeSequ.size(), 3 * nodeSequ.size()) - matrS;
		Matrix<double,Dynamic,Dynamic> K_h_E = matrP.transpose() * (elemVolu * matrG) * matrP 
			+ alph_E * matrH.transpose() * matrH;
		Matrix<double,Dynamic,Dynamic> tranK = 
			MatrixXd::Zero(3*nodeSequ.size(), 3*nodeSequ.size());
		for(map<long,long>::iterator iter_0 = nodeSequ.begin(); 
			iter_0 != nodeSequ.end(); iter_0 ++){
			tranK.block(3*(iter_0->second),3*(iter_0->second),3,3) = 
				freeTran(iter_0->first);
		}
		K_h_E = tranK.transpose() * K_h_E * tranK;
		for(map<long,long>::iterator iter_0 = nodeSequ.begin(); 
			iter_0 != nodeSequ.end(); iter_0 ++){
			for(long tj = 0; tj < 3; tj ++){
				if(consFlag(3 * iter_0->first + tj) != 0){
					long id_j = consDOF(3 * iter_0->first + tj);
					for(map<long,long>::iterator iter_1 = nodeSequ.begin(); 
						iter_1 != nodeSequ.end(); iter_1 ++){
						for(long tk = 0; tk < 3; tk ++){
							if(consFlag(3 * iter_1->first + tk) != 0){
								long id_k = consDOF(3 * iter_1->first + tk);
								stifList.emplace_back(id_j, id_k, 
									K_h_E(3*(iter_0->second)+tj , 3*(iter_1->second)+tk));
							}
						}
					}
				}
			}
		}
	}
	cout<<"There are "<<stifList.size()<<" nonzero components in stifList";
	OUTPUT_TIME("");
	return true;
}

bool VEM::OUTPUT_DISPLACEMENT(Matrix<double,Dynamic,Dynamic> tempDisp, long fileIden){
	stringstream tempStre;
	tempStre << "resuDisp_" << fileIden << ".txt";
	ofstream tempOfst(DIRECTORY(tempStre.str()), ios::out);
	tempStre.str("");
	tempStre.clear();
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	for(long ti = 0; ti < nodeCoor.rows(); ti ++){
		Vector3d nodeDisp;
		for(long tj = 0; tj < 3; tj ++){
			if(consFlag(3 * ti + tj)==0){
				nodeDisp(tj) = 0.0;
			}
			else{
				nodeDisp(tj) = tempDisp(consDOF(3*ti+tj));
			}
		}
		nodeDisp = freeTran(ti) * nodeDisp;
		tempOfst<<setw(30)<<nodeDisp(0)
			<<setw(30)<<nodeDisp(1)
			<<setw(30)<<nodeDisp(2)<<endl;
	}
	tempOfst.close();
	return true;
}

bool VEM::CONSTRAINT(const vector<long> & consFree){
	OUTPUT_TIME("Displacement constraint");
	long tempNumb = 0;
	consFlag.resize(3*nodeCoor.rows(),1);
	for(long ti = 0; ti < 3*nodeCoor.rows(); ti ++){
		consFlag(ti) = 1;
	}
	for(vector<long>::const_iterator iter_0 = consFree.begin(); 
		iter_0 != consFree.end(); iter_0 ++){
		consFlag(*iter_0) = 0;
		tempNumb ++;
	}
	freeNumb = 3 * nodeCoor.rows() - tempNumb;
	consDOF.resize(3 * nodeCoor.rows(),1);
	for(long ti = 0; ti < 3 * nodeCoor.rows(); ti ++){
		consDOF(ti) = 0;
	}
	for(long ti = 1; ti < 3 * nodeCoor.rows(); ti ++){
		consDOF(ti) = consDOF(ti-1) + consFlag(ti-1);
	}
	for(long ti = 1; ti < 3 * nodeCoor.rows(); ti ++){
		if(consFlag(ti) == 0){
			consDOF(ti) = -1;
		}
	}
	return true;
}

bool VEM::COUPLE(set<long> coupFree){
	long coupFlag = -1;
	long tempNumb = 0;
	for(long ti = 0; ti < 3 * nodeCoor.rows(); ti ++){
		set<long>::iterator iter_0 = coupFree.find(ti);
		if(iter_0 != coupFree.end()){
			if(consFlag(ti) == 0){
				cout<<"Error in VEM::COUPLE."<<endl;
			}
			if(coupFlag == -1){
				coupFlag = consDOF(ti);
			}
			else{
				tempNumb ++;
				consDOF(ti) = coupFlag;
			}
			consFlag(ti) = -1;
		}
		else{
			if(consFlag(ti) == 1){
				consDOF(ti) = consDOF(ti) - tempNumb;
			}
		}
	}
	freeNumb = freeNumb - tempNumb;
	coupFlag = -1;
	return true;
}

#endif
