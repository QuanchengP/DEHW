
#include "CONTACT.h"

class BLOCK: public CONTACT{
public:
	//
	//constraint and load settings
	//0 - case I with sparse mesh, 1 - case I with refined mesh
	//2 - case II with sparse mesh, 3 - case II with refined mesh
	//4 - case III with sparse mesh, 5 - case III with refined mesh
	long clse;
	Matrix<double,2,1> leng;//geometric length
	Matrix<long,2,1> diviNumb;//number of elements along different directions
	Matrix<long,2,1> buckNumb;//number for bucket search
	//
	BLOCK(long tempClse);
	bool MESH();//generate mesh
	//load on the non-overlapping part of lower block's upper surface
	bool LOAD_SUPPLEMENT(double pres, vector<pair<long,double>> & loadList);
	bool SOLVE();//solve for displacement and contact pressure
};

BLOCK::BLOCK(long tempClse){
	clse = tempClse;
	//
	mkdir("Block", 0755);
	outpDire = "Block/";
	//
	leng << 0.03, 0.025;
	if(clse % 2 == 0){
		diviNumb << 3, 2;
		buckNumb << 1, 1;//smaller than diviNumb
	}
	else{
		diviNumb << 19, 18;
		buckNumb << 100, 100;//smaller than diviNumb after refinement
	}
	//
	MESH();
}

bool BLOCK::MESH(){
	for(long tv = 0; tv < 2; tv ++){
		OUTPUT_TIME("Grid");
		//points
		for(long ti = 0; ti <= diviNumb(tv); ti ++){
			if(ti % 10 == 0){
				cout<<"Node layer "<<ti<<"/"<<diviNumb(tv)<<endl;
			}
			for(long tj = 0; tj <= diviNumb(tv); tj ++){
				for(long tk = 0; tk <= diviNumb(tv); tk ++){
					COORDINATE tempCoor;
					tempCoor <<
						leng(tv) / 2.0 - leng(tv) / diviNumb(tv) * (double)tj,
						leng(tv) / 2.0 - leng(tv) / diviNumb(tv) * (double)tk,
						- leng(tv) + leng(tv) / diviNumb(tv) * (double)ti;
					vem(tv).TRY_ADD_POINT(tempCoor);
				}
			}
		}
		//volumes
		long layeNumb = (diviNumb(tv)+1)*(diviNumb(tv)+1);
		for(long ti = 0; ti < diviNumb(tv); ti ++){
			if(ti % 10 == 0){
				cout<<"Volume layer "<<ti<<"/"<<diviNumb(tv)<<endl;
			}
			long tempI = ti * layeNumb;
			for(long tj = 0; tj < diviNumb(tv); tj ++){
				long tempJ = tempI + tj * (diviNumb(tv) + 1);
				for(long tk = 0; tk < diviNumb(tv); tk ++){
					Matrix<long,8,1> tempCorn;
					tempCorn << tempJ + tk + 1,
						tempJ + tk,
						tempJ + tk + (diviNumb(tv) + 1),
						tempJ + tk + 1 + (diviNumb(tv) + 1),
						tempJ + tk + 1 + layeNumb,
						tempJ + tk + layeNumb,
						tempJ + tk + (diviNumb(tv) + 1) + layeNumb,
						tempJ + tk + 1 + (diviNumb(tv) + 1) + layeNumb;
					vem(tv).ADD_BLOCK_FROM_POINT(tempCorn);
				}
			}
		}
		//refinement
		if(clse % 2 == 1){
			OUTPUT_TIME("Refine");
			list<long> spliVolu;
			Matrix<long,Dynamic,Dynamic> spliFlag;
			spliFlag.resize(diviNumb(tv)*diviNumb(tv),1);
			long tempNumb = (diviNumb(tv)-1)*diviNumb(tv)*diviNumb(tv);
			for(long ti=0; ti<diviNumb(tv); ti++){
				for(long tj=0; tj<diviNumb(tv); tj++){
					spliVolu.push_back(1+tempNumb+ti*diviNumb(tv)+tj);
					spliFlag(ti*diviNumb(tv)+tj) = 1;
				}
			}
			map<MATRIX,COORDINATE> planSurf;
			vem(tv).REFINE(spliVolu,spliFlag,planSurf);
			spliFlag.resize(spliVolu.size(),1);
			for(long ti=0; ti<spliVolu.size(); ti++){
				spliFlag(ti) = 1;
			}
			vem(tv).REFINE(spliVolu,spliFlag,planSurf);
			spliFlag.resize(spliVolu.size(),1);
			for(long ti=0; ti<spliVolu.size(); ti++){
				spliFlag(ti) = 1;
			}
			vem(tv).REFINE(spliVolu,spliFlag,planSurf);
			spliFlag.resize(spliVolu.size(),1);
			for(long ti=0; ti<spliVolu.size(); ti++){
				spliFlag(ti) = 1;
			}
			vem(tv).REFINE(spliVolu,spliFlag,planSurf);
		}
		vem(tv).VOLUME_2_ELEMENT();
		if(tv == 1){
			for(long ti = 0; ti < vem(tv).nodeCoor.rows(); ti ++){
				vem(tv).nodeCoor(ti,1) = - vem(tv).nodeCoor(ti,1);
				vem(tv).nodeCoor(ti,2) = - vem(tv).nodeCoor(ti,2);
			}
		}
		vem(tv).OUTPUT_ELEMENT(tv);
	}
	return true;
}

bool BLOCK::LOAD_SUPPLEMENT(double pres, vector<pair<long,double>> & loadList){
	list<pair<long,long>> tempList_0;
	for(long ti = 0; ti < vem(0).elemFace.rows(); ti ++){
		for(long tj = 0; tj < vem(0).elemFace(ti).rows(); tj ++){
			bool flag_ij = true;
			for(long tk = 0; tk < 3; tk ++){
				if(vem(0).nodeCoor(vem(0).elemFace(ti)(tj,tk),2) <= -1.0E-10){
					flag_ij = false;
				}
			}
			if(flag_ij){
				tempList_0.push_back(make_pair(ti,tj));
			}
		}
	}
	//16 triangles: the non-overlapping part of lower block's upper surface
	Matrix<Vector3d,16,3> tempList_1;
	tempList_1(0,0) << - leng(0) / 2.0, - leng(0) / 2.0, 0.0;
	tempList_1(0,1) << - leng(1) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(0,2) << - leng(0) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(1,0) << - leng(0) / 2.0, - leng(0) / 2.0, 0.0;
	tempList_1(1,1) << - leng(1) / 2.0, - leng(0) / 2.0, 0.0;
	tempList_1(1,2) << - leng(1) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(2,0) << - leng(1) / 2.0, - leng(0) / 2.0, 0.0;
	tempList_1(2,1) << leng(1) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(2,2) << - leng(1) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(3,0) << - leng(1) / 2.0, - leng(0) / 2.0, 0.0;
	tempList_1(3,1) << leng(1) / 2.0, - leng(0) / 2.0, 0.0;
	tempList_1(3,2) << leng(1) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(4,0) << leng(1) / 2.0, - leng(0) / 2.0, 0.0;
	tempList_1(4,1) << leng(0) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(4,2) << leng(1) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(5,0) << leng(1) / 2.0, - leng(0) / 2.0, 0.0;
	tempList_1(5,1) << leng(0) / 2.0, - leng(0) / 2.0, 0.0;
	tempList_1(5,2) << leng(0) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(6,0) << - leng(0) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(6,1) << - leng(1) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(6,2) << - leng(0) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(7,0) << - leng(0) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(7,1) << - leng(1) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(7,2) << - leng(1) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(8,0) << leng(1) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(8,1) << leng(0) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(8,2) << leng(1) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(9,0) << leng(1) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(9,1) << leng(0) / 2.0, - leng(1) / 2.0, 0.0;
	tempList_1(9,2) << leng(0) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(10,0) << - leng(0) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(10,1) << - leng(1) / 2.0, leng(0) / 2.0, 0.0;
	tempList_1(10,2) << - leng(0) / 2.0, leng(0) / 2.0, 0.0;
	tempList_1(11,0) << - leng(0) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(11,1) << - leng(1) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(11,2) << - leng(1) / 2.0, leng(0) / 2.0, 0.0;
	tempList_1(12,0) << - leng(1) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(12,1) << leng(1) / 2.0, leng(0) / 2.0, 0.0;
	tempList_1(12,2) << - leng(1) / 2.0, leng(0) / 2.0, 0.0;
	tempList_1(13,0) << - leng(1) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(13,1) << leng(1) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(13,2) << leng(1) / 2.0, leng(0) / 2.0, 0.0;
	tempList_1(14,0) << leng(1) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(14,1) << leng(0) / 2.0, leng(0) / 2.0, 0.0;
	tempList_1(14,2) << leng(1) / 2.0, leng(0) / 2.0, 0.0;
	tempList_1(15,0) << leng(1) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(15,1) << leng(0) / 2.0, leng(1) / 2.0, 0.0;
	tempList_1(15,2) << leng(0) / 2.0, leng(0) / 2.0, 0.0;
	for(list<pair<long,long>>::iterator iter_0 = tempList_0.begin(); 
		iter_0 != tempList_0.end(); iter_0 ++){
		for(long ti = 0; ti < 16; ti++){
			Matrix<Vector3d,3,1> nonmTria, mortTria;
			for(long tj = 0; tj < 3; tj ++){
				long iden_tj = vem(0).elemFace(iter_0->first)(iter_0->second,tj);
				nonmTria(tj) = vem(0).nodeCoor.block(iden_tj,0,1,3).transpose();
				mortTria(tj) = tempList_1(ti,tj);
			}
			Vector3d normVect;
			double origArea;
			Matrix<Vector2d,3,1> nonmProj, mortProj;
			vector<Vector2d> tempQuad;
			vector<double> tempWeig;
			//intersection of two triangles
			TI_SUB(nonmTria, mortTria, normVect, 
				origArea, nonmProj, mortProj, tempQuad, tempWeig);
			double triaArea = TRIANGLE_AREA_2D(nonmProj(0), nonmProj(1), nonmProj(2));
			for(long tk = 0; tk < tempQuad.size(); tk ++){
				Matrix<double,3,1> tempShap;
				for(long tj = 0; tj < 3; tj ++){
					tempShap((tj + 2) % 3) = 
						TRIANGLE_AREA_2D(nonmProj(tj), 
							nonmProj((tj + 1) % 3), tempQuad[tk]) / triaArea;
				}
				for(long tj = 0; tj < 3; tj ++){
					long iden_tj = vem(0).elemFace(iter_0->first)(iter_0->second,tj);
					loadList.push_back(make_pair(3 * iden_tj + 2,
						pres * tempShap(tj) * tempWeig[tk]));
				}
			}
		}
	}
	return true;
}

bool BLOCK::SOLVE(){
	//
	for(long tv = 0; tv < 2; tv ++){
		vem(tv).freeTran.resize(vem(tv).nodeCoor.rows(),1);
		for(long ti = 0; ti < vem(tv).nodeCoor.rows(); ti ++){
			vem(tv).freeTran(ti) << 1.0, 0.0, 0.0,
				0.0, 1.0, 0.0,
				0.0, 0.0, 1.0;
		}
	}
	//
	for(long tv = 0; tv < 2; tv ++){
		vector<long> consFree;
		for(long ti = 0; ti < vem(tv).nodeCoor.rows(); ti ++){
			if(clse == 0 || clse == 1 || clse == 4 || clse == 5){
				if(tv == 0 && vem(tv).nodeCoor(ti,2) <= - leng(tv) + 1.0E-10){
					consFree.push_back(3 * ti + 2);
				}
				if(vem(tv).nodeCoor(ti,0) <= - leng(tv) / 2.0 + 1.0E-10){
					consFree.push_back(3 * ti + 0);
				}
				if(vem(tv).nodeCoor(ti,1) <= - leng(tv) / 2.0 + 1.0E-10){
					consFree.push_back(3 * ti + 1);
				}
			}
			else{
				if(tv == 0 && vem(tv).nodeCoor(ti,2) <= - leng(tv) + 1.0E-10){
					consFree.push_back(3 * ti + 0);
					consFree.push_back(3 * ti + 1);
					consFree.push_back(3 * ti + 2);
				}
				if(tv == 1 && vem(tv).nodeCoor(ti,2) >= leng(tv) - 1.0E-10){
					consFree.push_back(3*ti+0);
					consFree.push_back(3*ti+1);
				}
			}
		}
		vem(tv).CONSTRAINT(consFree);
	}
	//
	OUTPUT_TIME("Load");
	Matrix<vector<pair<long,double>>,2,1> loadList;
	double contPres = - 1.0E7;
	for(long ti = 0; ti < vem(1).elemFace.rows(); ti ++){
		for(long tj = 0; tj < vem(1).elemFace(ti).rows(); tj ++){
			bool flag_ij = true;//the three nodes of the element face are all on z = leng(1)
			Matrix<double,3,3> tempXYZ = MatrixXd::Zero(3,3);
			for(long tk = 0; tk < 3; tk ++){
				tempXYZ.block(tk,0,1,3) = 
					vem(1).nodeCoor.block(vem(1).elemFace(ti)(tj,tk),0,1,3);
				if(vem(1).nodeCoor(vem(1).elemFace(ti)(tj,tk),2) < leng(1)-1.0E-10){
					flag_ij = false;
				}
			}
			if(flag_ij){
				//calculate the area of triangle by Heron's theorem
				double tempA = sqrt(pow(tempXYZ(0,0)-tempXYZ(1,0), 2.0) 
					+ pow(tempXYZ(0,1)-tempXYZ(1,1), 2.0) 
					+ pow(tempXYZ(0,2)-tempXYZ(1,2), 2.0));
				double tempB = sqrt(pow(tempXYZ(1,0)-tempXYZ(2,0), 2.0) 
					+ pow(tempXYZ(1,1)-tempXYZ(2,1), 2.0) 
					+ pow(tempXYZ(1,2)-tempXYZ(2,2), 2.0));
				double tempC = sqrt(pow(tempXYZ(2,0)-tempXYZ(0,0), 2.0) 
					+ pow(tempXYZ(2,1)-tempXYZ(0,1), 2.0) 
					+ pow(tempXYZ(2,2)-tempXYZ(0,2), 2.0));
				double tempP = (tempA + tempB + tempC)/2.0;
				double tempArea = sqrt(
					tempP * (tempP - tempA) * (tempP - tempB) * (tempP - tempC));
				double tempLoad = tempArea * contPres;
				for(long tk = 0; tk < 3; tk ++){
					loadList(1).push_back(
						make_pair(3 * vem(1).elemFace(ti)(tj,tk) + 2, tempLoad / 3.0));
				}
			}
		}
	}
	if(clse == 0 || clse == 1 || clse == 2 || clse == 3){
		LOAD_SUPPLEMENT(contPres, loadList(0));
	}
	//
	OUTPUT_TIME("Non-mortar and mortar segment");
	list<pair<long,long>> tempList_0;
	for(long ti = 0; ti < vem(0).elemFace.rows(); ti ++){
		for(long tj = 0; tj < vem(0).elemFace(ti).rows(); tj ++){
			bool flag_ij = true;
			for(long tk = 0; tk < 3; tk ++){
				if(vem(0).nodeCoor(vem(0).elemFace(ti)(tj,tk),2) <= -1.0E-10){
					flag_ij = false;
				}
			}
			if(flag_ij){
				tempList_0.push_back(make_pair(ti,tj));
			}
		}
	}
	nonmSegm.resize(tempList_0.size(),3);
	long tempIden = 0;
	for(list<pair<long,long>>::iterator iter_0 = tempList_0.begin(); 
		iter_0 != tempList_0.end(); iter_0 ++){
		nonmSegm.block(tempIden,0,1,3) << vem(0).elemFace(iter_0->first)(iter_0->second,0),
			vem(0).elemFace(iter_0->first)(iter_0->second,1),
			vem(0).elemFace(iter_0->first)(iter_0->second,2);
		tempIden ++;
	}
	tempList_0.clear();
	for(long ti = 0; ti < vem(1).elemFace.rows(); ti ++){
		for(long tj = 0; tj < vem(1).elemFace(ti).rows(); tj ++){
			bool flag_ij = true;
			for(long tk = 0; tk < 3; tk ++){
				if(vem(1).nodeCoor(vem(1).elemFace(ti)(tj,tk),2) >= 1.0E-10){
					flag_ij = false;
				}
			}
			if(flag_ij){
				tempList_0.push_back(make_pair(ti,tj));
			}
		}
	}
	mortSegm.resize(tempList_0.size(),3);
	tempIden = 0;
	for(list<pair<long,long>>::iterator iter_0 = tempList_0.begin(); 
		iter_0 != tempList_0.end(); iter_0 ++){
		mortSegm.block(tempIden,0,1,3) << vem(1).elemFace(iter_0->first)(iter_0->second,0),
			vem(1).elemFace(iter_0->first)(iter_0->second,1),
			vem(1).elemFace(iter_0->first)(iter_0->second,2);
		tempIden ++;
	}
	OUTPUT_CONTACT();
	//
	OUTPUT_TIME("Spatial search and contact detection");
	Matrix<double,Dynamic,Dynamic> locaCoor;
	locaCoor.resize(nonmSegm.rows(),2);
	for(long ti = 0; ti < nonmSegm.rows(); ti ++){
		locaCoor.block(ti,0,1,2) = (vem(0).nodeCoor.block(nonmSegm(ti,0),0,1,2).eval()
			+ vem(0).nodeCoor.block(nonmSegm(ti,1),0,1,2).eval()
			+ vem(0).nodeCoor.block(nonmSegm(ti,2),0,1,2).eval()) / 3.0;
	}
	MatrixXd tempFeat;
	tempFeat.resize(2,1);
	tempFeat << leng(0) / diviNumb(0), leng(0) / diviNumb(0);
	BUCKET_SORT(locaCoor, tempFeat, buckNumb(0), buckNumb(1));
	locaCoor.resize(mortSegm.rows(), 2);
	for(long ti = 0; ti < mortSegm.rows(); ti ++){
		locaCoor.block(ti,0,1,2) = (vem(1).nodeCoor.block(mortSegm(ti,0),0,1,2).eval()
			+ vem(1).nodeCoor.block(mortSegm(ti,1),0,1,2).eval()
			+ vem(1).nodeCoor.block(mortSegm(ti,2),0,1,2).eval()) / 3.0;
	}
	CONTACT_SEARCH(locaCoor);
	OUTPUT_CONTACT_ELEMENT();
	//
	Matrix<Vector3d,Dynamic,Dynamic> fricDire;
	fricDire.resize(0,0);
	CONTACT_ANALYSIS(loadList, fricDire);
	return true;
}
