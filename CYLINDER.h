
#include "CONTACT.h"

class CYLINDER: public CONTACT{
public:
	//
	//geometric settings
	//0 - case I, same length
	//1 - case II, different length, without corner rounding
	//2 - case III, different length, with corner rounding
	long gese;
	//0 - lower cylinder, 1 - upper cylinder
	Matrix<double,2,1> radi;//radius of cylinder
	Matrix<double,2,1> leng;//length of cylinder
	Matrix<long,2,5> diviNumb;//number of elements along different directions
	Matrix<long,2,1> buckNumb;//number for bucket search
	double bandWidt;//predicted contact band width > real value
	double loadInte;//load intensity
	double rounCorn;//radius of rounding corner
	double longIncr;//length difference / 2
	double longMini;//effective cylindrical length of upper cylinder = longMaxi - longMini
	double longMaxi;
	//element length of upper cylinder, used to calculate load vector
	Matrix<double,2,1> longDelt;
	//
	CYLINDER(long tempGese);
	bool MESH();//generate mesh
	//Cartesian - curvilinear - Cartesian
	bool PLNAE_2_SURFACE(long tv, const list<long> & spliVolu, 
		map<MATRIX,COORDINATE> & planSurf);
	//adjust a point on cylindrical onto the rounding corner
	COORDINATE ROUNDING_CORNER(COORDINATE tempCoor);
	//adjust a simplfy refined (straignt line/plane) point onto the rounding corner
	COORDINATE ROUNDING_CORNER_(long flag, double tempAngl, double tempThet);
	//is a node locating on the cylindrical surface, not inside, not on plane
	bool ON_SURFACE(long tv, COORDINATE tempCoor);
	bool SOLVE();//solve for displacement and contact pressure
};

CYLINDER::CYLINDER(long tempGese){
	gese = tempGese;
	//
	mkdir("Cylinder", 0755);
	outpDire = "Cylinder/";
	//
	radi << 0.022, 0.02;
	if(gese == 0){
		leng << 0.02, 0.02;
		diviNumb << 11, 14, 8, 31, 0,
			11, 14, 8, 40, 0;
		rounCorn = 0.0;
	}
	else if(gese == 1){
		leng << 0.02, 0.015;
		diviNumb << 11, 14, 8, 31, 0,
			11, 14, 8, 40, 0;
		rounCorn = 0.0;
	}
	else{
		leng << 0.02, 0.015;
		diviNumb << 11, 14, 8, 31, 5,
			11, 14, 8, 40, 5;
		rounCorn = 0.002;
	}
	buckNumb << 12, 512;
	bandWidt = 126E-6;
	loadInte = - 50.0E3;
	//
	MESH();
}

bool CYLINDER::MESH(){
	for(long tv = 0; tv < 2; tv ++){
		OUTPUT_TIME("Grid");
		//boundary nodes of block
		Matrix<double,Dynamic,Dynamic> upNode_0, downNode_0, upNode_1, downNode_1, upNode_2;
		upNode_0.resize(2, diviNumb(tv,0) + 1);
		downNode_0.resize(2, diviNumb(tv,0) + 1);
		upNode_1.resize(2, diviNumb(tv,1) + 1);
		downNode_1.resize(2, diviNumb(tv,1) + 1);
		upNode_2.resize(2, diviNumb(tv,1) + 1);
		double diviAngl_0 = - 5.0 / 8.0 * PI;
		double diviAngl_1 = - 3.0 / 8.0 * PI;
		Matrix<double,2,3> diviPoin;
		diviPoin << - radi(tv) / 2.0, - radi(tv) / 5.0, radi(tv) / 5.0,
			0.0, - radi(tv) * 2.0 / 3.0, - radi(tv) * 2.0 / 3.0;
		//1
		for(long tj = 0; tj <= diviNumb(tv,0); tj ++){
			upNode_0.block(0,tj,2,1) = 
				(1.0 - (double)tj / diviNumb(tv,0)) * diviPoin.block(0,0,2,1)
				+ (double)tj / diviNumb(tv,0) * diviPoin.block(0,1,2,1);
			double tempAngl = - PI + (diviAngl_0 + PI) * (double)tj / diviNumb(tv,0);
			downNode_0(0,tj) = radi(tv) * cos(tempAngl);
			downNode_0(1,tj) = radi(tv) * sin(tempAngl);
		}
		//2
		for(long tj = 0; tj <= diviNumb(tv,1); tj ++){
			upNode_1.block(0,tj,2,1) = 
				(1.0 - (double)tj / diviNumb(tv,1)) * diviPoin.block(0,1,2,1)
				+ (double)tj / diviNumb(tv,1) * diviPoin.block(0,2,2,1);
			double tempAngl = diviAngl_0 
				+ (diviAngl_1 - diviAngl_0) * (double)tj / diviNumb(tv,1);
			downNode_1(0,tj) = radi(tv) * cos(tempAngl);
			downNode_1(1,tj) = radi(tv) * sin(tempAngl);
			upNode_2(0,tj) = (1.0 - (double)tj / diviNumb(tv,1)) * (-radi(tv) / 2.0)
				+ (double)tj / diviNumb(tv,1) * (radi(tv) / 2.0);
			upNode_2(1,tj) = 0.0;
		}
		longIncr = (leng(0) - leng(1)) / 2.0;
		longMini = longIncr + rounCorn;
		longMaxi = leng(0) - longIncr - rounCorn;
		longDelt(0) = rounCorn / diviNumb(1,4);
		longDelt(1) = (leng(1) - 2*rounCorn) / (diviNumb(1,3) - 2 * diviNumb(1,4));
		//points
		Matrix<Matrix<long,Dynamic,Dynamic>,4,1> blocPoin;
		blocPoin(0).resize(diviNumb(tv,0) + 1, diviNumb(tv,2) + 1);
		blocPoin(1).resize(diviNumb(tv,1) + 1, diviNumb(tv,2) + 1);
		blocPoin(2).resize(diviNumb(tv,0) + 1, diviNumb(tv,2) + 1);
		blocPoin(3).resize(diviNumb(tv,1) + 1, diviNumb(tv,0) + 1);
		for(long tj = 0; tj <= diviNumb(tv,3); tj ++){
			double coor_z;
			if(tv == 0){
				coor_z = (double)tj / diviNumb(tv,3) * leng(tv);
			}
			else{
				if(tj <= diviNumb(tv,4) - 1){
					coor_z = longIncr + rounCorn * (double)tj / diviNumb(tv,4);
				}
				else if(tj >= diviNumb(tv,3) - diviNumb(tv,4) + 1){
					coor_z = longMaxi + rounCorn * 
						(double)(tj - (diviNumb(tv,3) - diviNumb(tv,4))) / diviNumb(tv,4);
				}
				else{
					coor_z = longIncr + rounCorn 
						+ (double)(tj - diviNumb(tv,4)) / (diviNumb(tv,3) - 2 * diviNumb(tv,4)) 
						* (leng(tv) - 2 * rounCorn);
				}
			}
			//1-1
			for(long tk = 0; tk <= diviNumb(tv,2); tk ++){
				for(long tm = 0; tm <= diviNumb(tv,0); tm ++){
					COORDINATE tempCoor;
					tempCoor <<
						(double)tk / diviNumb(tv,2) * upNode_0(0,tm)
						+ (1.0 - (double)tk / diviNumb(tv,2)) * downNode_0(0,tm),
						(double)tk / diviNumb(tv,2) * upNode_0(1,tm)
						+ (1.0 - (double)tk / diviNumb(tv,2)) * downNode_0(1,tm),
						coor_z;
					if(tv == 1){
						tempCoor = ROUNDING_CORNER(tempCoor);
					}
					long tempIden = vem(tv).TRY_ADD_POINT(tempCoor);
					if(tj == 0){
						blocPoin(0)(tm,tk) = tempIden;
					}
				}
			}
			//2
			for(long tk = 0; tk <= diviNumb(tv,2); tk ++){
				for(long tm = 0; tm <= diviNumb(tv,1); tm ++){
					COORDINATE tempCoor;
					tempCoor <<
						(double)tk / diviNumb(tv,2) * upNode_1(0,tm)
						+ (1.0 - (double)tk / diviNumb(tv,2)) * downNode_1(0,tm),
						(double)tk / diviNumb(tv,2) * upNode_1(1,tm)
						+ (1.0 - (double)tk / diviNumb(tv,2)) * downNode_1(1,tm),
						coor_z;
					if(tv == 1){
						tempCoor = ROUNDING_CORNER(tempCoor);
					}
					long tempIden = vem(tv).TRY_ADD_POINT(tempCoor);
					if(tj == 0){
						blocPoin(1)(tm,tk) = tempIden;
					}
				}
			}
			//1-2
			for(long tk = 0; tk <= diviNumb(tv,2); tk ++){
				for(long tm = 0; tm <= diviNumb(tv,0); tm ++){
					COORDINATE tempCoor;
					tempCoor <<
						- ((double)tk / diviNumb(tv,2) * upNode_0(0, diviNumb(tv,0) - tm)
						+ (1.0 - (double)tk / diviNumb(tv,2)) 
						* downNode_0(0, diviNumb(tv,0) - tm)),
						(double)tk / diviNumb(tv,2) * upNode_0(1, diviNumb(tv,0) - tm)
						+ (1.0 - (double)tk / diviNumb(tv,2)) 
						* downNode_0(1, diviNumb(tv,0) - tm),
						coor_z;
					if(tv == 1){
						tempCoor = ROUNDING_CORNER(tempCoor);
					}
					long tempIden = vem(tv).TRY_ADD_POINT(tempCoor);
					if(tj == 0){
						blocPoin(2)(tm,tk) = tempIden;
					}
				}
			}
			//3
			for(long tk = 0; tk <= diviNumb(tv,0); tk ++){
				for(long tm = 0; tm <= diviNumb(tv,1); tm ++){
					COORDINATE tempCoor;
					tempCoor <<
						(double)tk / diviNumb(tv,0) * upNode_2(0,tm)
						+ (1.0 - (double)tk / diviNumb(tv,0)) * upNode_1(0,tm),
						(double)tk / diviNumb(tv,0) * upNode_2(1,tm)
						+ (1.0 - (double)tk / diviNumb(tv,0)) * upNode_1(1,tm),
						coor_z;
					if(tv == 1){
						tempCoor = ROUNDING_CORNER(tempCoor);
					}
					long tempIden = vem(tv).TRY_ADD_POINT(tempCoor);
					if(tj == 0){
						blocPoin(3)(tm,tk) = tempIden;
					}
				}
			}
		}
		//volumes
		long layeNumb = (1 + 2 * diviNumb(tv,0) + diviNumb(tv,1)) * (diviNumb(tv,2) + 1)
			+ (diviNumb(tv,1) - 1) * diviNumb(tv,0);
		for(long tj = 0; tj < diviNumb(tv,3); tj ++){
			long idTemp_j = layeNumb * tj;
			Matrix<long,8,1> tempCorn;
			//1-1
			for(long tk = 0; tk < diviNumb(tv,2); tk ++){
				for(long tm = 0; tm < diviNumb(tv,0); tm ++){
					tempCorn << idTemp_j + blocPoin(0)(tm, tk),
						idTemp_j + blocPoin(0)(tm + 1, tk),
						idTemp_j + blocPoin(0)(tm + 1, tk + 1),
						idTemp_j + blocPoin(0)(tm, tk + 1),
						idTemp_j + blocPoin(0)(tm, tk) + layeNumb,
						idTemp_j + blocPoin(0)(tm + 1, tk) + layeNumb,
						idTemp_j + blocPoin(0)(tm + 1, tk + 1) + layeNumb,
						idTemp_j + blocPoin(0)(tm, tk + 1) + layeNumb;
					vem(tv).ADD_BLOCK_FROM_POINT(tempCorn);
				}
			}
			//2
			for(long tk = 0; tk < diviNumb(tv,2); tk ++){
				for(long tm = 0; tm < diviNumb(tv,1); tm ++){
					tempCorn << idTemp_j + blocPoin(1)(tm, tk),
						idTemp_j + blocPoin(1)(tm + 1, tk),
						idTemp_j + blocPoin(1)(tm + 1, tk + 1),
						idTemp_j + blocPoin(1)(tm, tk + 1),
						idTemp_j + blocPoin(1)(tm, tk) + layeNumb,
						idTemp_j + blocPoin(1)(tm + 1, tk) + layeNumb,
						idTemp_j + blocPoin(1)(tm + 1, tk + 1) + layeNumb,
						idTemp_j + blocPoin(1)(tm, tk + 1) + layeNumb;
					vem(tv).ADD_BLOCK_FROM_POINT(tempCorn);
				}
			}
			//1-2
			for(long tk = 0; tk < diviNumb(tv,2); tk ++){
				for(long tm = 0; tm < diviNumb(tv,0); tm ++){
					tempCorn << idTemp_j + blocPoin(2)(tm, tk),
						idTemp_j + blocPoin(2)(tm + 1, tk),
						idTemp_j + blocPoin(2)(tm + 1, tk + 1),
						idTemp_j + blocPoin(2)(tm , tk + 1),
						idTemp_j + blocPoin(2)(tm, tk) + layeNumb,
						idTemp_j + blocPoin(2)(tm + 1, tk) + layeNumb,
						idTemp_j + blocPoin(2)(tm + 1, tk + 1) + layeNumb,
						idTemp_j + blocPoin(2)(tm, tk + 1) + layeNumb;
					vem(tv).ADD_BLOCK_FROM_POINT(tempCorn);
				}
			}
			//3
			for(long tk = 0; tk < diviNumb(tv,0); tk ++){
				for(long tm = 0; tm < diviNumb(tv,1); tm ++){
					tempCorn << idTemp_j + blocPoin(3)(tm, tk),
						idTemp_j + blocPoin(3)(tm + 1, tk),
						idTemp_j + blocPoin(3)(tm + 1, tk + 1),
						idTemp_j + blocPoin(3)(tm, tk + 1),
						idTemp_j + blocPoin(3)(tm, tk) + layeNumb,
						idTemp_j + blocPoin(3)(tm + 1, tk) + layeNumb,
						idTemp_j + blocPoin(3)(tm + 1, tk + 1) + layeNumb,
						idTemp_j + blocPoin(3)(tm, tk + 1) + layeNumb;
					vem(tv).ADD_BLOCK_FROM_POINT(tempCorn);
				}
			}
		}
		//refinement
		OUTPUT_TIME("Refine");
		list<long> spliVolu;
		Matrix<long,Dynamic,Dynamic> spliFlag;
		spliFlag.resize(2 * diviNumb(tv, 3), 1);
		layeNumb = (2 * diviNumb(tv,0) + diviNumb(tv,1)) * diviNumb(tv,2)
			+ diviNumb(tv,1) * diviNumb(tv,0);
		long tempNumb = 0;
		for(long tj = 0; tj < diviNumb(tv,3); tj ++){
			for(long tk = diviNumb(tv,1) / 2; tk <= diviNumb(tv,1) / 2 + 1; tk ++){
				spliVolu.push_back(tj * layeNumb + diviNumb(tv,0) * diviNumb(tv,2) + tk);
				spliFlag(tempNumb) = -1;
				if(tk == diviNumb(tv,1) / 2){
					spliFlag(tempNumb) = 7;
				}
				else if(tk == diviNumb(tv,1) / 2 + 1){
					spliFlag(tempNumb) = 8;
				}
				tempNumb ++;
			}
		}
		map<MATRIX,COORDINATE> planSurf;
		PLNAE_2_SURFACE(tv, spliVolu, planSurf);
		vem(tv).REFINE(spliVolu, spliFlag, planSurf);//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		spliFlag.resize(spliVolu.size(),1);
		for(long ti = 0; ti < spliVolu.size(); ti ++){
			if(ti % 4 <= 1){
				spliFlag(ti) = 7;
			}
			else{
				spliFlag(ti) = 8;
			}
		}
		planSurf.clear();
		PLNAE_2_SURFACE(tv, spliVolu, planSurf);
		vem(tv).REFINE(spliVolu, spliFlag, planSurf);//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		spliFlag.resize(spliVolu.size(),1);
		for(long ti = 0; ti < spliVolu.size(); ti ++){
			if(ti % 8 <= 3){
				spliFlag(ti) = 3;
			}
			else{
				spliFlag(ti) = 3;
			}
		}
		planSurf.clear();
		PLNAE_2_SURFACE(tv, spliVolu, planSurf);
		vem(tv).REFINE(spliVolu, spliFlag, planSurf);//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		spliFlag.resize(spliVolu.size(),1);
		for(long ti = 0; ti < spliVolu.size(); ti ++){
			if(ti % 32 <= 15){
				if((ti % 32) % 4 == 1 || (ti % 32) % 4 == 3){
					spliFlag(ti) = 3;
				}
				else{
					spliFlag(ti) = -1;
				}
			}
			else{
				if((ti % 32 - 16) % 4 ==0 || (ti % 32 - 16) % 4 == 2){
					spliFlag(ti) = 3;
				}
				else{
					spliFlag(ti) = -1;
				}
			}
		}
		planSurf.clear();
		PLNAE_2_SURFACE(tv, spliVolu, planSurf);
		vem(tv).REFINE(spliVolu, spliFlag, planSurf);//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		spliFlag.resize(spliVolu.size(),1);
		for(long ti = 0; ti < spliVolu.size(); ti ++){
			spliFlag(ti) = 3;
		}
		planSurf.clear();
		PLNAE_2_SURFACE(tv, spliVolu, planSurf);
		vem(tv).REFINE(spliVolu, spliFlag, planSurf);//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		spliFlag.resize(spliVolu.size(),1);
		for(long ti = 0; ti < spliVolu.size(); ti ++){
			if(ti % 256 <= 127){
				if((ti % 256) % 8 == 1 || (ti % 256) % 8 >= 3){
					spliFlag(ti) = 3;
				}
				else{
					spliFlag(ti) = 7;
				}
			}
			else{
				if((ti % 256 - 128) % 8 <= 4 || (ti % 256 - 128) % 8 == 6){
					spliFlag(ti) = 3;
				}
				else{
					spliFlag(ti) = 8;
				}
			}
		}
		planSurf.clear();
		PLNAE_2_SURFACE(tv, spliVolu, planSurf);
		vem(tv).REFINE(spliVolu, spliFlag, planSurf);//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		spliFlag.resize(spliVolu.size(),1);
		for(long ti = 0; ti < spliVolu.size(); ti ++){
			spliFlag(ti) = -1;
		}
		planSurf.clear();
		PLNAE_2_SURFACE(tv, spliVolu, planSurf);
		vem(tv).REFINE(spliVolu, spliFlag, planSurf);//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		vem(tv).VOLUME_2_ELEMENT();
		if(tv == 0){
			for(long tj = 0; tj < vem(tv).nodeCoor.rows(); tj ++){
				vem(tv).nodeCoor(tj,0) = - vem(tv).nodeCoor(tj,0);
				vem(tv).nodeCoor(tj,1) = - vem(tv).nodeCoor(tj,1) - radi(0)- radi(1);
			}
		}
		vem(tv).OUTPUT_ELEMENT(tv);
	}
	return true;
}

bool CYLINDER::PLNAE_2_SURFACE(long tv, const list<long> & spliVolu, 
	map<MATRIX,COORDINATE> & planSurf){
	for(list<long>::const_iterator iter_0 = spliVolu.begin(); 
		iter_0 != spliVolu.end(); iter_0 ++){
		long tempVolu = *iter_0;
		map<long,MATRIX>::iterator iter_V = vem(tv).voluFace.find(tempVolu);
		for(long ti = 0; ti < (iter_V->second).rows(); ti ++){
			long tempFace = abs((iter_V->second)(ti));
			map<long,MATRIX>::iterator iter_F = vem(tv).faceLine.find(tempFace);
			bool faceFlag = true;
			COORDINATE faceCoor;
			faceCoor(0) = faceCoor(1) = faceCoor(2) = 0.0;
			MATRIX faceMatr;
			faceMatr.resize(4,1);
			for(long tj = 0; tj < (iter_F->second).rows(); tj ++){
				long tempLine = abs((iter_F->second)(tj));
				map<long,MATRIX>::iterator iter_L = vem(tv).linePoin.find(tempLine);
				long facePoin;
				if((iter_F->second)(tj) > 0){
					facePoin = (iter_L->second)(0);
				}
				else{
					facePoin = (iter_L->second)(1);
				}
				faceMatr(tj) = facePoin; 
				map<long,COORDINATE>::iterator face_P = vem(tv).poinCoor.find(facePoin);
				if(ON_SURFACE(tv, face_P->second)==false){
					faceFlag = false;
				}
				faceCoor = faceCoor + face_P->second;
				bool lineFlag = true;
				COORDINATE lineCoor;
				lineCoor(0) = lineCoor(1) = lineCoor(2) = 0.0;
				MATRIX lineMatr;
				lineMatr.resize(2,1);
				for(long tk = 0; tk < 2; tk ++){
					long linePoin = (iter_L->second)(tk);
					lineMatr(tk) = linePoin;
					map<long,COORDINATE>::iterator line_P = vem(tv).poinCoor.find(linePoin);
					if(ON_SURFACE(tv, line_P->second) == false){
						lineFlag = false;
					}
					lineCoor = lineCoor + line_P->second;
				}
				if(lineFlag){
					lineCoor = lineCoor / 2.0;
					COORDINATE resuCoor;
					double tempAngl = atan2(lineCoor(1),lineCoor(0));
					if(tv == 0 || 
						(longMini - 1.0E-10 <= lineCoor(2) 
						&& lineCoor(2) <= longMaxi + 1.0E-10)){
						resuCoor(0) = radi(tv) * cos(tempAngl);
						resuCoor(1) = radi(tv) * sin(tempAngl);
						resuCoor(2) = lineCoor(2);
					}
					else if(lineCoor(2) < longMini-1.0E-10){
						double tempThet = atan2(longMini - lineCoor(2), 
							- (radi(1) - rounCorn) - lineCoor(1));
						resuCoor = ROUNDING_CORNER_(0, tempAngl, tempThet);
					}
					else{
						double tempThet = atan2(lineCoor(2) - longMaxi, 
							-(radi(1) - rounCorn) - lineCoor(1));
						resuCoor = ROUNDING_CORNER_(1, tempAngl, tempThet);
					}
					planSurf.insert(map<MATRIX,COORDINATE>::value_type(
						lineMatr, resuCoor));//automatically avoid duplication
				}
			}
			if(faceFlag){
				faceCoor = faceCoor / 4.0;
				COORDINATE resuCoor;
				double tempAngl = atan2(faceCoor(1), faceCoor(0));
				if(tv == 0 || 
					(longMini - 1.0E-10 <= faceCoor(2) && faceCoor(2) <= longMaxi + 1.0E-10)){
					resuCoor(0) = radi(tv) * cos(tempAngl);
					resuCoor(1) = radi(tv) * sin(tempAngl);
					resuCoor(2) = faceCoor(2);
				}
				else if(faceCoor(2) < longMini - 1.0E-10){
					double tempThet = atan2(longMini - faceCoor(2), 
						- (radi(1) - rounCorn) - faceCoor(1));
					resuCoor = ROUNDING_CORNER_(0, tempAngl, tempThet);
				}
				else{
					double tempThet = atan2(faceCoor(2) - longMaxi, 
						- (radi(1) - rounCorn) - faceCoor(1));
					resuCoor = ROUNDING_CORNER_(1, tempAngl, tempThet);
				}
				planSurf.insert(map<MATRIX,COORDINATE>::value_type(
					faceMatr, resuCoor));//automatically avoid duplication
			}
		}
	}
	return true;
}

COORDINATE CYLINDER::ROUNDING_CORNER(COORDINATE tempCoor){
	COORDINATE resuCoor = tempCoor;
	if(abs(rounCorn) > 1.0E-12){
		if(tempCoor(2) < longMini -1.0E-10){
			double targRadi = radi(1) - rounCorn 
				+ sqrt(pow(rounCorn, 2) - pow(longMini - tempCoor(2), 2));
			resuCoor(0) = tempCoor(0) * targRadi / radi(1);
			resuCoor(1) = tempCoor(1) * targRadi / radi(1);
		}
		else if(tempCoor(2) > longMaxi +1.0E-10){
			double targRadi = radi(1) - rounCorn 
				+ sqrt(pow(rounCorn, 2) - pow(tempCoor(2) - longMaxi, 2));
			resuCoor(0) = tempCoor(0) * targRadi / radi(1);
			resuCoor(1) = tempCoor(1) * targRadi / radi(1);
		}
	}
	return resuCoor;
}

COORDINATE CYLINDER::ROUNDING_CORNER_(long flag, double tempAngl, double tempThet){
	COORDINATE resuCoor;
	if(tempThet > PI/2.0){
		tempThet = PI/2.0;
	}
	if(flag == 0){
		double targRadi = radi(1) - rounCorn + rounCorn * cos(tempThet);
		resuCoor(0) = targRadi * cos(tempAngl);
		resuCoor(1) = targRadi * sin(tempAngl);
		resuCoor(2) = longMini - rounCorn * sin(tempThet);
	}
	else{
		double targRadi = radi(1) - rounCorn + rounCorn * cos(tempThet);
		resuCoor(0) = targRadi * cos(tempAngl);
		resuCoor(1) = targRadi * sin(tempAngl);
		resuCoor(2) = longMaxi + rounCorn * sin(tempThet);
	}
	return resuCoor;
}

bool CYLINDER::ON_SURFACE(long tv, COORDINATE tempCoor){
	double tempRadi = sqrt(tempCoor(0) * tempCoor(0) + tempCoor(1) * tempCoor(1));
	if(tv == 0 || (longMini - 1.0E-10 <= tempCoor(2) && tempCoor(2) <= longMaxi + 1.0E-10)){
		if(abs(tempRadi - radi(tv)) > 1.0E-10){
			return false;
		}
		else{
			return true;
		}
	}
	else if(tempCoor(2) < longMini-1.0E-10){
		double targRadi = radi(1) - rounCorn 
			+ sqrt(pow(rounCorn, 2) - pow(longMini - tempCoor(2), 2));
		if(abs(tempRadi - targRadi) > 1.0E-10){
			return false;
		}
		else{
			return true;
		}
	}
	else{
		double targRadi = radi(1) - rounCorn 
			+ sqrt(pow(rounCorn, 2) - pow(tempCoor(2) - longMaxi, 2));
		if(abs(tempRadi - targRadi) > 1.0E-10){
			return false;
		}
		else{
			return true;
		}
	}
}

bool CYLINDER::SOLVE(){
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
			if(tv == 0 && vem(tv).nodeCoor(ti,1) <= - radi(0) - radi(1) + 1.0E-10){
				consFree.push_back(3 * ti + 0);
				consFree.push_back(3 * ti + 1);
				consFree.push_back(3 * ti + 2);
			}
			else if(tv == 1 && vem(tv).nodeCoor(ti,1) >= -1.0E-10){
				consFree.push_back(3 * ti + 0);
				consFree.push_back(3 * ti + 2);
			}
		}
		vem(tv).CONSTRAINT(consFree);
	}
	//
	OUTPUT_TIME("Load");
	Matrix<vector<pair<long,double>>,2,1> loadList;
	for(long ti = 0; ti < vem(1).nodeCoor.rows(); ti ++){
		if(vem(1).nodeCoor(ti,1) >= -1.0E-10 && abs(vem(1).nodeCoor(ti,0)) <= 1.0E-10){
			if(abs(rounCorn) < 1.0E-12){
				if(vem(1).nodeCoor(ti,2) <= longMini+1.0E-10){
					loadList(1).push_back(make_pair(3 * ti + 1, loadInte * longDelt(1) / 2.0));
				}
				else if(vem(1).nodeCoor(ti,2) <= longMaxi-1.0E-10){
					loadList(1).push_back(make_pair(3 * ti + 1, loadInte * longDelt(1)));
				}
				else{
					loadList(1).push_back(make_pair(3 * ti + 1, loadInte * longDelt(1) / 2.0));
				}
			}
			else{
				if(vem(1).nodeCoor(ti,2) <= longIncr + 1.0E-10){
					loadList(1).push_back(make_pair(3 * ti + 1, loadInte * longDelt(0) / 2.0));
				}
				else if(vem(1).nodeCoor(ti,2) <= longMini - 1.0E-10){
					loadList(1).push_back(make_pair(3 * ti + 1, loadInte * longDelt(0)));
				}
				else if(vem(1).nodeCoor(ti,2) <= longMini + 1.0E-10){
					loadList(1).push_back(make_pair(3 * ti + 1, 
						loadInte * (longDelt(0) + longDelt(1)) / 2.0));
				}
				else if(vem(1).nodeCoor(ti,2) <= longMaxi - 1.0E-10){
					loadList(1).push_back(make_pair(3 * ti + 1, loadInte * longDelt(1)));
				}
				else if(vem(1).nodeCoor(ti,2) <= longMaxi + 1.0E-10){
					loadList(1).push_back(make_pair(3 * ti + 1, 
						loadInte * (longDelt(0) + longDelt(1)) / 2.0));
				}
				else if(vem(1).nodeCoor(ti,2) <= leng(0) - longIncr - 1.0E-10){
					loadList(1).push_back(make_pair(3 * ti + 1, loadInte * longDelt(0)));
				}
				else{
					loadList(1).push_back(make_pair(3 * ti + 1, loadInte * longDelt(0) / 2.0));
				}
			}
		}
	}
	//
	OUTPUT_TIME("Non-mortar and mortar segment");
	list<pair<long,long>> tempList_0;
	for(long ti = 0; ti < vem(0).elemFace.rows(); ti ++){
		for(long tj = 0; tj < vem(0).elemFace(ti).rows(); tj ++){
			bool flag_ij = true;
			for(long tk = 0; tk < 3; tk ++){
				double tempX = vem(0).nodeCoor(vem(0).elemFace(ti)(tj,tk),0);
				double tempY = vem(0).nodeCoor(vem(0).elemFace(ti)(tj,tk),1);
				double tempZ = vem(0).nodeCoor(vem(0).elemFace(ti)(tj,tk),2);
				if(abs(tempX) > bandWidt){
					flag_ij = false;
					break;
				}
				double tempRadi = sqrt(tempX * tempX
					+ (tempY + radi(0) + radi(1)) * (tempY + radi(0) + radi(1)));
				if(abs(tempRadi - radi(0)) > 1.0E-10){
					flag_ij = false;
					break;
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
				double tempX = vem(1).nodeCoor(vem(1).elemFace(ti)(tj,tk),0);
				double tempY = vem(1).nodeCoor(vem(1).elemFace(ti)(tj,tk),1);
				double tempZ = vem(1).nodeCoor(vem(1).elemFace(ti)(tj,tk),2);
				if(abs(tempX) > bandWidt){
					flag_ij = false;
					break;
				}
				double tempRadi = sqrt(tempX * tempX + tempY * tempY);
				if(longMini - 1.0E-10 <= tempZ && tempZ <= longMaxi + 1.0E-10){
					if(abs(tempRadi - radi(1)) > 1.0E-10){
						flag_ij = false;
						break;
					}
				}
				else if(tempZ < longMini - 1.0E-10){
					double targRadi = radi(1) - rounCorn 
						+ sqrt(pow(rounCorn, 2) - pow(longMini - tempZ, 2));
					if(abs(tempRadi - targRadi) > 1.0E-10 || radi(1) - tempRadi > 2.0E-6){
						flag_ij = false;
						break;
					}
				}
				else if(tempZ > longMaxi + 1.0E-10){
					double targRadi = radi(1) - rounCorn 
						+ sqrt(pow(rounCorn, 2) - pow(tempZ - longMaxi, 2));
					if(abs(tempRadi - targRadi) > 1.0E-10 || radi(1) - tempRadi > 2.0E-6){
						flag_ij = false;
						break;
					}
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
		locaCoor(ti,0) = (vem(0).nodeCoor(nonmSegm(ti,0),0)
			+ vem(0).nodeCoor(nonmSegm(ti,1),0)
			+ vem(0).nodeCoor(nonmSegm(ti,2),0)) / 3.0;
		locaCoor(ti,1) = (vem(0).nodeCoor(nonmSegm(ti,0),2)
			+ vem(0).nodeCoor(nonmSegm(ti,1),2)
			+ vem(0).nodeCoor(nonmSegm(ti,2),2)) / 3.0;
	}
	MatrixXd tempFeat;
	tempFeat.resize(2,1);
	tempFeat << bandWidt / buckNumb(0), leng(0) / buckNumb(1);//!!!@@@###
	BUCKET_SORT(locaCoor, tempFeat, buckNumb(0), buckNumb(1));
	locaCoor.resize(mortSegm.rows(),2);
	for(long ti = 0; ti < mortSegm.rows(); ti ++){
		locaCoor(ti,0) = (vem(1).nodeCoor(mortSegm(ti,0),0)
			+ vem(1).nodeCoor(mortSegm(ti,1),0)
			+ vem(1).nodeCoor(mortSegm(ti,2),0)) / 3.0;
		locaCoor(ti,1) = (vem(1).nodeCoor(mortSegm(ti,0),2)
			+ vem(1).nodeCoor(mortSegm(ti,1),2)
			+ vem(1).nodeCoor(mortSegm(ti,2),2)) / 3.0;
	}
	CONTACT_SEARCH(locaCoor);
	OUTPUT_CONTACT_ELEMENT();
	//
	Matrix<Vector3d,Dynamic,Dynamic> fricDire;
	fricDire.resize(0,0);
	CONTACT_ANALYSIS(loadList, fricDire);
	return true;
}
