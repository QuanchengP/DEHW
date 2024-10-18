
#include "EINM.h"
#include "CONTACT.h"

/*************************************************************************************************/
/*******************************************DECLARATION*******************************************/
/*************************************************************************************************/

class FUNC_DATA;

class DEHW: public CONTACT{
public:
	//****************************************************************************************
	DEHW();
	//****************************************************************************************
	Matrix<long,2,1> z;//teeth number
	double a_h2;//working center distance
	double modiTran;//modification of transmission ratio
	double modiCent;//modification of center distance
	double r_b2;//radius of base circle of worm wheel
	double beta_c;//inclination angle of generating plane
	//****************************************************************************************
	double z_k;//encircled teeth
	Matrix<double,2,1> d;//reference circle diameter at throat
	Matrix<double,2,1> R_a;//tip arc
	double offsR_a;//offset of tip arc center of worm wheel
	double wheeWidt;//width of worm wheel
	Matrix<double,2,1> inneRadi;//inner rdius of hub
	//****************************************************************************************
	double a_1c;//center distance of primary enveloping
	double i_1c;//transmission ratio of primary enveloping
	double i_c1;
	double i_h2;//transmission ratio of second enveloping
	double i_2h;
	double m_t;//transverse module
	Matrix<double,2,1> h_a_s;//addendum
	Matrix<double,2,1> h_f_s;//dedendum
	Matrix<double,2,1> h_a;//addendum
	Matrix<double,2,1> h_f;//dedendum
	Matrix<double,2,1> d_f;//root circle
	Matrix<double,2,1> d_a;//tip circle
	Matrix<double,2,1> R_f;//root arc
	Matrix<double,2,1> R_t;//transition arc
	double alph;//nominal pressure angle
	double leadAngl;//nominal lead angle
	double pitcAngl;//pitch angle
	Matrix<double,2,1> tootThicCoef;
	double halfAngl;//half working angle
	double starAngl;//starting angle
	double termAngl;//terminating angle
	Matrix<double,3,1> wormCurv;//curvilinear coordinate of worm
	double widtAngl;//face width angle of worm wheel
	double backlash;//backlash
	//axial tooth thickness of worm, transverse tooth thickness of worm wheel
	Matrix<double,2,1> tootThic;
	Matrix<double,2,1> tootThicAngl;//tooth thickness angle
	//simultaneous envelope of tooth surface and tooth back
	//the angle between the two rigidly connected coordinates
	Matrix<double,2,1> backAngl;
	long BASIC_PARAMETER();
	//****************************************************************************************
	MatrixXd thet_ci; 
	MatrixXd thet_hi;
	MatrixXd thet_cu;
	Matrix<bool,Dynamic,Dynamic> f_cu;
	long GEOMETRY();
	long PRESET();//calculate thet_ci, thet_hi, thet_cu, f_cu
	//ICL of the first envelope
	long ICLOFE(long V, double x_dL, double x_dH, 
		double thet_c, double thet_h, ofstream &tempOfst);
	//ICL of the second envelope
	long ICLOSE(double thet_h, ofstream &tempOfst);
	//meshing limit line of the first envelope
	long MLLOFE(long V, double thet_cL, double thet_cH, ofstream &tempOfst);
	//curvature interference limit line of the first envelope
	long CILLOFE(long V, double thet_cL, double thet_cH, 
		MatrixXd &thet_cr, MatrixXd &thet_hr, ofstream &tempOfst);
	//meshing limit line of the second envelope
	long MLLOSE(double thet_cL, double thet_cH, ofstream &tempOfst);
	//curvature interference limit line of the second envelope
	long CILLOSE(double thet_cL, double thet_cH, 
		const MatrixXd &thet_cr, const MatrixXd &thet_hr, ofstream &tempOfst);
	//singular thet_h to thet_c
	long SINGULAR_H2C(double thet_h, double &thet_cmini, double &thet_cmaxi);
	//singular thet_c to thet_h
	long SINGULAR_C2H(double thet_c, double &thet_hs, double &thet_hm);
	//first and second meshing equations
	long FSME(double thet_1, double thet_h, double &x_d, double &y_d);
	//partial derivative of meshing equations
	long PD_FSME(double thet_1, double thet_h, 
		double &x_d, double &y_d, Matrix<double,2,2> &Pxy_d);
	//x_d, y_d, thet_c to r_1
	long WORM_DC2R(double x_d, double y_d, double thet_c, Vector3d &r_1_1);
	//x_d, y_d, thet_1, thet_h to r_2
	long WHEE_1H2R(double x_d, double y_d, double thet_1, double thet_h, 
		Vector3d &r_2_2);
	//partial derivative of r_2_2 relative to thet_1, thet_h
	long PD_WHEE_1H2R(double x_d, double y_d, double thet_1, double thet_h, 
		Vector3d &r_2_2, Matrix<double,3,2> &Dr_2_2);
	//curvature interference limit function of first envelope
	long CILFOFE(double thet_1, double x_d, double y_d,
		double &Psi_1, double &kapp_cxd, double &kapp_cyd, double &tau_cxd);
	//partial derivative for -curvature interference limit function of first envelope
	long PD_CILFOFE(double thet_1, double x_d, double y_d, double Px_d, double Py_d, 
		double &Psi_1, double &kapp_cxd, double &kapp_cyd, double &tau_cxd, 
		double &PPsi_1, double &Pkapp_1xd, double &Pkapp_1yd, double &Ptau_1xd);
	//function data setup, curvature interference limit function of second envelope
	long FUDA_SETUP(const double thet_c, FUNC_DATA &fuda);
	//curvature interference limit function of second envelope, non-interval version
	double CILFOSE_NI(double thet_1, double thet_h, double &kapp_h2N);
	//derivative for curvature interference limit function of second envelope, 
	//non-interval version
	double PD_CILFOSE_NI(double thet_1, double thet_h);
	//on worm wheel tooth surface
	long ON_WHEEL_TS(double x_d, double y_d, double thet_1, Vector3d r_2_2);
	//total length of ICL, must call WHEE_TS_GRID() in advance
	double TICL_LENGTH(double thet_h, const Matrix<double,Dynamic,Dynamic> &sd12);
	//analytical contact pressure, must call WHEE_TS_GRID() in advance
	long AC_PRESSURE(double thet_1, double thet_h, double &sigm_H, double &b_H);
	//analytical contact analysis
	long ANALYTICAL_CONTACT();
	//****************************************************************************************
	long refiLeve;//refinement level
	long reliSwit;//0 - no tooth flank relif, 1 - tooth flank relief
	Matrix<long,2,6> gridNumb;//number of grid division
	//curvilinear coordinate of tooth profile
	Matrix<Matrix<Vector2d,Dynamic,Dynamic>,2,1> curvCoor;
	//Cartesian coordinate of tooth profile
	Matrix<Matrix<Vector3d,Dynamic,Dynamic>,2,1> cartCoor;
	//flag of worm wheel tooth surface
	//1 - left new contact zone, 2 - right new contact zone, 3 - former contact zone
	//4 - head transition zone,  5 - rear transition zone,   0 - fail to solve
	Matrix<int,Dynamic,Dynamic> fpha;
	long TS_GRID();//grid discretization of tooth surface
	long WORM_TS_GRID();//grid discretization of worm tooth surface
	long WHEE_TS_GRID();//grid discretization of worm wheel tooth surface
	//worm tooth surface: curvilinear coordinate to Cartesian coordinate
	long WORM_CURV_2_CART(double xi_11, double xi_12, Vector3d &r_1_1, double &thet_c);
	//global coordinate r_2_2 to local coordinate of worm wheel
	long WHEE_G2L(Vector3d r_2_2, double &angl_f, double &radi_f, 
		double &R_fmini, double &R_fmaxi);
	//worm wheel tooth surface: curvilinear coordinate to Cartesian coordinate
	//thet_c, thet_h are initial values, new contact zone
	long WHEE_CURV_2_CART_1(double xi_21, double xi_22, Vector3d &r_2_2, 
		double &thet_c, double &thet_h, long f_lr, double &x_d, double &y_d);
	//x_d, y_d are initial values, former contact zone
	long WHEE_CURV_2_CART_2(double xi_21, double xi_22, Vector3d &r_c_c, 
		double &thet_c, double &x_d, double &y_d);
	//transition zone
	long WHEE_CURV_2_CART_3(double xi_21, double xi_22, Vector3d &r_2_2, 
		double &thet_c, double &thet_h, double xi_11);
	//transition zone of worm wheel, xi_11 - head transition zone, rear transition zone
	long WHEE_TRAN(double thet_c, double thet_h, double xi_11, 
		Vector3d &r_2_2, Matrix<double,3,2> &Dr_2_2);
	//phase analysis of worm wheel tooth surface
	long WHEE_PHAS(long ti, long tj, long f_ij, Vector3d r_2_2);
	//new contact zone, 1 - left, 2 - right
	long NEW_CONT_ZONE(long f_lr);
	long FORMER_CONT_ZONE();//former contact zone
	long TRANSITION_ZONE(long f_hr);//transition zone, 1 - head, 2 - rear
	//tooth flank relief of worm
	long WORM_RELI(Vector3d &tempXYZ, long ti, long tj);
	//tooth flank relief of worm wheel
	long WHEE_RELI(Vector3d &tempXYZ, long ti, long tj);
	//****************************************************************************************
	//discretization point on profile of worm, tooth surface and tooth back
	Matrix<Matrix<double,3,2>,Dynamic,Dynamic> wormProf;
	//discretization point on profile of worm wheel, tooth surface and tooth back
	Matrix<Matrix<double,3,2>,Dynamic,Dynamic> wheeProf;
	Matrix<double,Dynamic,Dynamic> profCurv;//-xi_11
	//four node forms one segment, the segment is which volume's face
	Matrix<vector<pair<MATRIX,long>>,2,1> segmVolu;
	//point coordinate to the numbering in curvCoor/cartCoor
	Matrix<map<COORDINATE,Matrix<long,2,3>>,2,1> poinNumb;
	long WORM_GRID();//generate grid of worm
	long WHEE_GRID();//generate grid of worm wheel
	//calculate the radius of root transition arc
	long WORM_ROOT_RADIUS(long flag, Matrix<double,2,3> tempPoin, 
		Matrix<double,2,1> &tempCent, double &tempRadi, Matrix<double,2,1> &tempAngl);
	//calculate nodes on root transition arc
	long WORM_ROOT(long id, long flag, Matrix<double,Dynamic,Dynamic> &rootProf);
	//transfer from "in unfolded cone surface" to r_2_2
	Matrix<double,3,1> WHEE_CONE(Matrix<double,2,1> tempXY, double tempAlph_3);
	//calculate nodes on root transition arc "in unfolded cone surface"
	long WHEE_ROOT(Matrix<double,Dynamic,Dynamic> &rootProf, 
		Matrix<Matrix<double,2,1>,Dynamic,Dynamic> tempProf, 
		long id, double r_f, double tempPitc);
	//****************************************************************************************
	double centErro;//center distance error
	Matrix<double,6,1> distCrit;//critical gap used for adaptive mesh refinement
	double fricCoef;//frictional coefficient
	double inpuTorq;//input torque
	long analNumb;//contact analysis at analNumb meshing positions
	Matrix<double,2,1> analAngl;//rotating angle of worm and worm wheel
	//bucket search
	Matrix<Matrix<double,2,3>,2,1> buckLoca_R;//local coordinate
	Matrix<long,2,2> buckNumb_R;//number of buckets
	Matrix<Matrix<list<long>,Dynamic,Dynamic>,2,1> buck_R;//bucket
	//use tempCrit as critical gap for the adaptive mesh refinement
	//higher refinement level, more buckets
	bool ADAPTIVE_REFINE(double tempCrit, long buckFact);
	//put non-mortar and mortar segments into different buckets
	bool BUCKET_SORT_R(long buckFact);
	//distance from point to volume
	double DISTANCE(COORDINATE tempCoor, long voluIden);
	//distance from point to triangle
	double TRIANGLE_DISTANCE(COORDINATE tempCoor, COORDINATE tempVect_0, 
		COORDINATE tempVect_1, COORDINATE tempVect_2);
	//Cartesian - curvilinear - Cartesian
	COORDINATE PLNAE_2_SURFACE(long tv, MATRIX tempMatr, 
		Matrix<map<COORDINATE,Matrix<long,2,3>>,2,1> & poinNumb_0);
	long SOLVE();//solve for displacement and contact pressure
};

//function data
class FUNC_DATA{
public:
	double thet_1;
	double C_11;
	double C_12;
	double C_13;
	double C_m11;
	double C_m12;
	double C_m13;
	double C_a11;
	double thet_a1;
	double C_m21;
	double C_m22;
	double C_m23;
	double C_a21;
	double thet_a2;
	double C_m31;
	double C_m32;
	double C_m33;
	double C_a31;
	double thet_a3;
	double C_f21;
	double C_f22;
	double C_f23;
	double C_a41;
	double C_a42;
	double C_a51;
	double C_a52;
	double thet_a5;
	double C_a61;
	double C_a62;
	double C_a63;
	double C_b11;
	double C_b12;
	double thet_b1;
	double C_b21;
	double C_b22;
	double C_b31;
	double C_b32;
	double thet_b3;
	double C_b41;
	double C_b42;
	double thet_b4;
	double C_b51;
	double C_b52;
	double thet_b5;
	double C_b61;
	double C_b62;
	double C_b71;
	double C_b72;
	double C_b81;
	double thet_b81;
	double thet_b82;
	double C_b83;
	double C_c11;
	double C_c12;
	double C_c21;
	double C_c22;
	double thet_c2;
	double C_c31;
	double C_c32;
	double C_c41;
	double C_c42;
	double thet_c4;
	double C_c43;
	double C_c51;
	double C_c52;
	double thet_c5;
	double C_c53;
	double C_c54;
	double C_c61;
	double C_s21;
	double C_s22;
	double C_s23;
	double C_c62;
	double thet_c6;
	//
	double thet_h;
	double C_s11;
	double C_s12;
	double C_s13;
};

class INDE_INIT{
public:
	long ti;
	long tj;
	double init_1;
	double init_2;
	INDE_INIT(long inde_1, long inde_2, double valu_1, double valu_2)
		: ti(inde_1), tj(inde_2), init_1(valu_1), init_2(valu_2){}
};

//curvature interference limit function of second envelope, interval version
FINTERVAL CILFOSE(const FINTERVAL &THET_H, const FUNC_DATA &fuda);
//refinement version
FINTERVAL CILFOSE_(const FINTERVAL &THET_H, const FUNC_DATA &fuda);

//thet_h is known, X is unknwon
FINTERVAL CILFOSE_X(const FINTERVAL &X, const FUNC_DATA &fuda);

//derivative for curvature interference limit function of second envelope, interval version
FINTERVAL PD_CILFOSE(const FINTERVAL &THET_H, const FUNC_DATA &fuda);
//refinement version
FINTERVAL PD_CILFOSE_(const FINTERVAL &THET_H, const FUNC_DATA &fuda);

//thet_h is known, X is unknwon
FINTERVAL PD_CILFOSE_X(const FINTERVAL &X, const FUNC_DATA &fuda);

/*************************************************************************************************/
/*****************************************IMPLEMENTATION******************************************/
/*************************************************************************************************/

//
DEHW::DEHW(){
	//
	mkdir("Dehw", 0755);
	outpDire = "Dehw/";
	FINTERVAL::precision(16); //Set interval output precision to 16 in Filib++
	//
	z << 1, 40;
	a_h2 = 0.25;
	modiTran = 0.2;
	modiCent = 0.0014;
	r_b2 = 0.158/2.0;
	beta_c = 11.0 * PI / 180.0;
	//
	z_k = 4.2;
	d(0) = 0.082;
	h_a_s << 0.6, 0.7;//smaller, unused
	h_f_s << 0.95, 1.05;//larger, larger
	R_a(1) = 0.0385;
	offsR_a = 0.003;
	wheeWidt = 0.06;
	inneRadi << 0.018, 0.15;
	inpuTorq = 180.0;//N*m
	fricCoef = 0.08;
	//
	reliSwit = 0;
	centErro = 0.0E-6;//30.0E-6/0.0
	distCrit << 34.0E-6, 31.0E-6, 28.0E-6, 25.0E-6, 22.0E-6, 1.0E-6;
	//
	refiLeve = 0;
	vem(1).mateElas = 110.0E9;
	BASIC_PARAMETER();
}

long DEHW::BASIC_PARAMETER(){
	//
	a_1c = a_h2 + modiCent;
	i_h2 = (double)(z(1)) / z(0);
	i_1c = i_h2 + modiTran;
	i_c1 = 1.0 / i_1c;
	i_2h = 1.0 / i_h2;
	//
	d(1) = 2.0 * a_h2 - d(0);
	m_t = d(1) / z(1);
	h_a(0) = h_a_s(0) * m_t;
	h_a(1) = h_a_s(1) * m_t;
	h_f(0) = h_f_s(0) * m_t;
	h_f(1) = h_f_s(1) * m_t;
	d_f(0) = d(0) - 2.0 * h_f(0);
	d_f(1) = d(1) - 2.0 * h_f(1);
	d_a(0) = d(0) + 2.0 * h_a(0);
	d_a(1) = d(1) + 2.0 * h_a(1);
	R_a(0) = a_h2 - 0.5 * d_a(0);
	R_f(0) = a_h2 - 0.5 * d_f(0);
	R_f(1) = a_h2 - 0.5 * d_f(1);
	R_t(0) = a_h2 - 0.5 * d(0) + 0.8 * m_t;//larger
	R_t(1) = a_h2 - 0.5 * d(1) + 0.9 * m_t;//larger
	//
	alph = asin(2.0 * r_b2 / d(1));
	leadAngl = atan(d(1) / i_h2 / d(0));
	pitcAngl = 2.0 * PI / z(1);
	tootThicCoef << 0.45, 0.55;
	halfAngl = 0.5 * (z_k - tootThicCoef(0)) * pitcAngl;
	starAngl = alph - halfAngl;
	termAngl = starAngl + z_k * pitcAngl;
	wormCurv(0) = i_h2 * starAngl;
	wormCurv(2) = i_h2 * termAngl;
	wormCurv(1) = (wormCurv(0) + wormCurv(2)) / 2.0;
	while(wormCurv(1) - 2.0 * PI >= wormCurv(0)){
		wormCurv(1) = wormCurv(1) - 2.0 * PI;
	}
	//
	widtAngl = asin(wheeWidt / 2.0 / R_f(1));
	backlash = 0.0;
	tootThic(0) = tootThicCoef(0) * PI * m_t - backlash;
	tootThic(1) = tootThicCoef(1) * PI * m_t;
	tootThicAngl(0) = tootThic(0) / (d(1) / 2.0);
	tootThicAngl(1) = tootThic(1) / (d(1) / 2.0);
	backAngl(0) = 2.0 * alph + tootThicAngl(0);
	backAngl(1) = 2.0 * alph - tootThicAngl(1);
	return 1;
}

long DEHW::GEOMETRY(){
	//
	PRESET();
	
	//ICL of the second envelope
	ofstream tempOfst;
	tempOfst.open(DIRECTORY("resuICLOSE.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	long tiMaxi = i_h2 / 6 * z(0);
	for(long ti = 1; ti <= tiMaxi; ti ++){
		double thet_h = ti * 2.0 * PI / z(0);
		if(thet_h < PI){
			continue;
		}
		ICLOSE(thet_h, tempOfst);
	}
	for(long ti = 0; ti < thet_ci.rows(); ti ++){
		ICLOSE(thet_hi(ti), tempOfst);
	}
	tempOfst.close();
	
	//meshing limit line of the first envelope
	tempOfst.open(DIRECTORY("resuMLLOFE.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	// double thet_cL = i_c1 * 1.1 * PI;
	// double thet_cH = i_c1 * (i_h2 / 3) * PI;
	double thet_cL = 0.01 * PI;
	double thet_cH = 0.49 * PI;
	MLLOFE(1, thet_cL, thet_cH, tempOfst);
	tempOfst.close();
	
	//curvature interference limit line of the first envelope
	tempOfst.open(DIRECTORY("resuCILLOFE.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	MatrixXd thet_cr, thet_hr;
	CILLOFE(1, thet_cL, thet_cH, thet_cr, thet_hr, tempOfst);
	tempOfst.close();
	
	//meshing limit line of the second envelope
	tempOfst.open(DIRECTORY("resuMLLOSE.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	MLLOSE(thet_cL, thet_cH, tempOfst);
	tempOfst.close();
	
	//curvature interference limit line of the second envelope
	tempOfst.open(DIRECTORY("resuCILLOSE.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	CILLOSE(thet_cL, thet_cH, thet_cr, thet_hr, tempOfst);
	tempOfst.close();
	//curvature interference limit function of the second envelope
	/*tempOfst.open(DIRECTORY("resuCILFOSE.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	long tj = 120;
	cout<<"thet_c="<<thet_cr(tj)<<"thet_hr="<<thet_hr(tj,4)<<endl;
	for(long tk = 0; tk <= 99999; tk ++){
		double thet_c = thet_cr(tj);
		double thet_h = thet_hr(tj, 0)+1E-9 
			+ (thet_hr(tj, 4)-1E-9 - thet_hr(tj, 0)-1E-9) / 99999 * tk;
		double kapp_h2N;
		tempOfst<<setw(30)<<thet_h
			<<setw(30)<<CILFOSE_NI(i_1c * thet_c, thet_h, kapp_h2N)<<endl;
	}
	tempOfst.close();*/
	
	return 1;
}

long DEHW::PRESET(){
	if(modiTran != 0.0 || modiCent != 0.0){
		double C_i1 = (i_c1 * i_c1 * i_2h * i_2h * (a_h2 * a_h2 - a_1c * a_1c) 
			- (i_2h * a_h2 - i_c1 * a_1c) * (i_2h * a_h2 - i_c1 * a_1c)) 
			* cos(beta_c) * cos(beta_c);
		double C_i2 = 2.0 * i_c1 * i_2h * (i_2h * a_h2 - i_c1 * a_1c)
			* a_h2 * cos(beta_c) * sin(beta_c);
		double C_i3 = (i_2h * a_h2 - i_c1 * a_1c) * (i_2h * a_h2 - i_c1 * a_1c);
		double Delt = C_i2 * C_i2 - 4.0 * C_i1 * C_i3;
		MatrixXd solu;
		if(Delt == 0.0){
			solu.resize(1,1);
			solu(0) = - C_i2 / 2.0 / C_i1;
		}
		else if(Delt > 0.0){
			solu.resize(2,1);
			solu(0) = (- C_i2 - sqrt(Delt)) / 2.0 / C_i1;
			solu(1) = (- C_i2 + sqrt(Delt)) / 2.0 / C_i1;
		}
		MatrixXd thet_ci1, thet_hi1;
		thet_ci1.resize(solu.rows(), 1);
		thet_hi1.resize(solu.rows(), 1);
		long numb_c = 0;
		for(long ti = 0; ti < solu.rows(); ti ++){
			if(solu(ti) <= 0.0 || solu(ti) >= 1.0){
				continue;
			}
			double thet_c = acos(solu(ti));
			double thet_1 = i_1c * thet_c;
			double thet_h, thet_hs, thet_hm;
			thet_h = thet_1 + atan2((i_2h * a_h2 - i_c1 * a_1c) * sin(thet_c), 
				(i_2h * a_h2 - i_c1 * a_1c) * tan(beta_c) + i_c1 * i_2h * a_h2 * cos(thet_c));
			SINGULAR_C2H(thet_c, thet_hs, thet_hm);
			while(thet_h < thet_hs){
				thet_h = thet_h + 2.0 * PI;
			}
			while(thet_h >= thet_hs + 2.0 * PI){
				thet_h = thet_h - 2.0 * PI;
			}
			cout<<"ICL_INFINITE_"<<ti<<": thet_c="<<setprecision(20)<<thet_c
				<<", thet_h="<<thet_h<<endl;
			cout<<"                thet_hs="<<setprecision(20)<<thet_hs
				<<", thet_hm="<<thet_hm<<endl;
			thet_ci1(numb_c) = thet_c;
			thet_hi1(numb_c) = thet_h;
			numb_c ++;
		}
		thet_ci.resize(numb_c, 1);
		thet_hi.resize(numb_c, 1);
		thet_ci.block(0, 0, numb_c, 1) = thet_ci1.block(0, 0, numb_c, 1);
		thet_hi.block(0, 0, numb_c, 1) = thet_hi1.block(0, 0, numb_c, 1);
	}
	else{
		thet_ci.resize(0, 0);
		thet_hi.resize(0, 0);
	}
	double C_u1 = cos(beta_c) * cos(beta_c);
	double C_u2 = - 3.0 * i_c1 * cos(beta_c) * sin(beta_c);
	double C_u3 = i_c1 * i_c1 * sin(beta_c) * sin(beta_c) 
		- 2.0 * i_c1 * i_c1 * cos(beta_c) * cos(beta_c) - 1.0;
	double C_u4 = i_c1 * i_c1 * i_c1 * cos(beta_c) * sin(beta_c);
	MatrixXd solu;
	CUBIC(C_u1, C_u2, C_u3, C_u4, solu);
	MatrixXd thet_cu1;
	thet_cu1.resize(solu.rows(), 1);
	long numb_c = 0;
	for(long ti = 0; ti < solu.rows(); ti ++){
		if(solu(ti) <= 0.0 || solu(ti) >= 1.0){
			continue;
		}
		double thet_c = acos(solu(ti));
		thet_cu1(numb_c) = thet_c;
		numb_c ++;
	}
	thet_cu.resize(numb_c, 1);
	thet_cu.block(0,0,numb_c,1) = thet_cu1.block(0,0,numb_c,1);
	f_cu.resize(numb_c, 1);
	for(long ti = 0; ti < thet_cu.rows(); ti ++){
		double thet_c = thet_cu(ti);
		Matrix<double,2,2> coef;
		coef << - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c),
			- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c),
			cos(beta_c) * cos(thet_c) * sin(thet_c),
			r_b2 * cos(beta_c) * sin(thet_c) * sin(thet_c)
			+ i_c1 * i_c1 * r_b2 * cos(beta_c) - a_1c * cos(beta_c) * sin(thet_c);
		if(abs(coef.determinant()) > 1.0E-12){
			f_cu(ti) = false;
			cout<<"CILLOFE_NONEXISTENCE: thet_c="<<setprecision(20)<<thet_c<<endl;
			continue;
		}
		f_cu(ti) = true;
		cout<<"CILLOFE_INFINITE: thet_c="<<setprecision(20)<<thet_c<<endl;
	}
	return 1;
}

long DEHW::ICLOFE(long V, double x_dL, double x_dH, 
	double thet_c, double thet_h, ofstream &tempOfst){
	long N = 9999;
	double deltX_d = (x_dH - x_dL) / (N + 1.0);
	for(long tj = 1; tj <= N; tj ++){
		double x_d = x_dL + deltX_d * (double)tj;
		double C_11 = - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c);
		double C_12 = sin(thet_c);
		double C_13 = - r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c);
		double y_d = (- C_13 - C_11 * x_d) / C_12;
		if(V == 1){
			Vector3d r_1_1;
			WORM_DC2R(x_d, y_d, thet_c, r_1_1);
			tempOfst<<setw(30)<<r_1_1(0)
				<<setw(30)<<r_1_1(1)<<setw(30)<<r_1_1(2)<<endl;
		}
		else{
			Vector3d r_2_2;
			WHEE_1H2R(x_d, y_d, i_1c * thet_c, thet_h, r_2_2);
			tempOfst<<setw(30)<<r_2_2(0)
				<<setw(30)<<r_2_2(1)<<setw(30)<<r_2_2(2)<<endl;
		}
	}
	return 1;
}

long DEHW::ICLOSE(double thet_h, ofstream &tempOfst){
	//
	if(thet_h <= PI){
		cout<<"ERROR in DEHW::ICLOSE!"<<endl;
		return -1;
	}
	//
	double thet_hsH, thet_hm;
	SINGULAR_C2H(PI / 2.0, thet_hsH, thet_hm);
	if(thet_h >= thet_hsH){
		cout<<"ERROR in DEHW::ICLOSE!"<<endl;
		return -1;
	}
	//
	if(modiTran == 0.0 && modiCent == 0.0){
		ICLOFE(2, - a_1c, a_1c, i_c1 * thet_h, thet_h, tempOfst);
		tempOfst<<"#"<<endl;
	}
	//
	for(long ti = 0; ti < thet_ci.rows(); ti ++){
		if(abs(thet_hi(ti) - thet_h) < 1.0E-14){
			ICLOFE(2, - a_1c, a_1c, thet_ci(ti), thet_hi(ti), tempOfst);
			tempOfst<<"#"<<endl;
		}
	}
	//
	double thet_cL, thet_cH;
	SINGULAR_H2C(thet_h, thet_cL, thet_cH);
	double thet_hmR;
	//
	long N = 9999;
	double deltThet_c = (thet_cH - thet_cL) / (N + 1.0);
	default_random_engine dren;
	uniform_real_distribution<double> urdd(0.25, 0.75);
	for(long tj = 1; tj <= N; tj ++){
		//
		double thet_c = thet_cL + deltThet_c * (double)tj;
		//
		double thet_hm;
		while(true){
			double thet_hstemp;
			SINGULAR_C2H(thet_c, thet_hstemp, thet_hm);
			if(abs(thet_h - thet_hm) <= 1.0E-14){
				thet_c = thet_cL + deltThet_c * (tj - urdd(dren));
			}
			else{
				break;
			}
		}
		//
		if(tj >= 2 && ((thet_hmR < thet_h && thet_hm > thet_h) 
			|| (thet_hmR > thet_h && thet_hm < thet_h))){
			tempOfst<<"#"<<endl;
		}
		thet_hmR = thet_hm;
		//
		double thet_1 = i_1c * thet_c;
		double x_d, y_d;
		FSME(thet_1, thet_h, x_d, y_d);
		Vector3d r_2_2;
		WHEE_1H2R(x_d, y_d, thet_1, thet_h, r_2_2);
		tempOfst<<setw(30)<<r_2_2(0)
			<<setw(30)<<r_2_2(1)<<setw(30)<<r_2_2(2)
			<<setw(30)<<thet_h<<setw(30)<<thet_c<<endl;
	}
	tempOfst<<"@"<<endl;
	return 1;
}

long DEHW::MLLOFE(long V, double thet_cL, double thet_cH, ofstream &tempOfst){
	long N = 9999;
	//
	double deltThet_c = (thet_cH - thet_cL) / (N + 1.0);
	for(long ti = 1; ti <= N; ti ++){
		double thet_c = thet_cL + deltThet_c * ti;
		Matrix<double,2,2> coeA;
		coeA << - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c), sin(thet_c),
			i_c1 * sin(beta_c) * sin(thet_c), i_c1 * cos(thet_c);
		Matrix<double,2,1> coeB;
		coeB << - r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c),
			- i_c1 * r_b2 * sin(beta_c) * cos(thet_c);
		Matrix<double,2,1> xy_d = coeA.fullPivLu().solve(-coeB);
		Vector3d r;
		if(V == 1){
			WORM_DC2R(xy_d(0), xy_d(1), thet_c, r);
		}
		else{
			WHEE_1H2R(xy_d(0), xy_d(1), i_1c * thet_c, i_1c * thet_c, r);
		}
		tempOfst<<setw(30)<<r(0)
			<<setw(30)<<r(1)<<setw(30)<<r(2)<<endl;
	}	
	return 1;
}

long DEHW::CILLOFE(long V, double thet_cL, double thet_cH, 
	MatrixXd &thet_cr, MatrixXd &thet_hr, ofstream &tempOfst){
	long N = 999;
	if(V == 1){
		thet_cr.resize(N, 1);
		thet_hr.resize(N, 5);
	}
	for(long ti = 0; ti < thet_cu.rows(); ti ++){
		if(f_cu(ti) == true){
			ICLOFE(V, - a_1c, a_1c, thet_cu(ti), i_1c * thet_cu(ti), tempOfst);
			if(V == 1){
				tempOfst<<"#"<<endl;
			}
		}
	}
	double deltThet_c = (thet_cH - thet_cL) / (N + 1);
	for(long tj = 1; tj <= N; tj ++){
		double thet_c = thet_cL + deltThet_c * tj;
		for(long ti = 0; ti < thet_cu.rows(); ti ++){
			if(abs(thet_cu(ti) - thet_c) < 1.0E-14){
				thet_c = thet_cL + deltThet_c * (tj - 0.5);
			}
		}
		if(V == 1 && tj >= 2){
			for(long ti = 0; ti< thet_cu.rows(); ti ++){
				if(thet_cr(tj - 2) < thet_cu(ti) && thet_cu(ti) < thet_c){
					tempOfst<<"#"<<endl;
					break;
				}
			}
		}
		Matrix<double,2,2> coeA;
		coeA << - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c), sin(thet_c),
			cos(beta_c) * cos(thet_c) * sin(thet_c),
			cos(beta_c) * sin(beta_c) * cos(thet_c) * cos(thet_c)
			+ 2.0 * i_c1 * cos(beta_c) * cos(beta_c) * cos(thet_c)
			- i_c1 * i_c1 * cos(beta_c) * sin(beta_c);
		Matrix<double,2,1> coeB;
		coeB << - r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c),
			r_b2 * cos(beta_c) * sin(thet_c) * sin(thet_c)
			+ i_c1 * i_c1 * r_b2 * cos(beta_c) - a_1c * cos(beta_c) * sin(thet_c);
		Matrix<double,2,1> xy_d = coeA.fullPivLu().solve(-coeB);
		Vector3d r;
		if(V == 1){
			WORM_DC2R(xy_d(0), xy_d(1), thet_c, r);
		}
		else{
			WHEE_1H2R(xy_d(0), xy_d(1), i_1c * thet_c, i_1c * thet_c, r);
		}
		tempOfst<<setw(30)<<r(0)
			<<setw(30)<<r(1)<<setw(30)<<r(2)<<endl;
		//
		if(V == 1){
			thet_cr(tj - 1) = thet_c;
			double thet_1 = i_1c * thet_c;
			double C_31 = i_2h * xy_d(0) * cos(beta_c) 
				- i_2h * a_1c * cos(beta_c) * cos(thet_c);
			double C_32 = - i_2h * xy_d(0) * sin(beta_c) * sin(thet_c)
				- i_2h * xy_d(1) * cos(thet_c)
				+ i_2h * r_b2 * sin(beta_c) * cos(thet_c);
			double C_33 = xy_d(0) * sin(beta_c) * cos(thet_c) - xy_d(1) * sin(thet_c)
				+ r_b2 * sin(beta_c) * sin(thet_c)
				- a_1c * sin(beta_c) + i_2h * a_h2 * cos(beta_c) * cos(thet_c);
			double thet_hs, thet_hm;
			SINGULAR_C2H(thet_c, thet_hs, thet_hm);
			thet_hr(tj - 1, 0) = thet_hs;
			thet_hr(tj - 1, 1) = thet_hm;
			thet_hr(tj - 1, 4) = thet_hs + 2.0 * PI;
			if(C_33 * C_33 > C_31 * C_31 + C_32 * C_32){
				thet_hr(tj - 1, 2) = thet_hs;
				thet_hr(tj - 1, 3) = thet_hs;
			}
			else{
				double thet_h0 = thet_1 - atan2(C_31, C_32) 
					- asin(C_33 / sqrt(C_31 * C_31 + C_32 * C_32));
				while(thet_h0 < thet_hs){
					thet_h0 += 2.0 * PI;
				}
				while(thet_h0 >= thet_hs + 2.0 * PI){
					thet_h0 -= 2.0 * PI;
				}
				double thet_h1 = thet_1 - atan2(C_31, C_32) 
					+ PI + asin(C_33 / sqrt(C_31 * C_31 + C_32 * C_32));
				while(thet_h1 < thet_hs){
					thet_h1 += 2.0 * PI;
				}
				while(thet_h1 >= thet_hs + 2.0 * PI){
					thet_h1 -= 2.0 * PI;
				}
				thet_hr(tj - 1, 2) = thet_h0;
				thet_hr(tj - 1, 3) = thet_h1;
			}
		}
	}
	if(V == 1){
		thet_hr = SORT(thet_hr, 0, true);
	}
	return 1;
}

long DEHW::MLLOSE(double thet_cL, double thet_cH, ofstream &tempOfst){
	if(modiTran == 0.0 && modiCent == 0.0){
		MLLOFE(2, thet_cL, thet_cH, tempOfst);
	}
	else{
		for(long ti = 0; ti< thet_ci.rows(); ti ++){
			double thet_c = thet_ci(ti);
			double thet_1 = i_1c * thet_c;
			//singular limit point at D_3 locates on the mllose
			double slml = (i_c1 * a_1c - i_2h * a_h2) 
				* (i_2h * a_h2 - i_c1 * a_1c 
				+ i_c1 * i_2h * a_h2 * cos(beta_c) * sin(beta_c) * cos(thet_ci(ti)));
			cout<<"singular limit point at D_3 locates on the mllose="
				<<setprecision(20)<<slml<<endl;
			double thet_h = thet_hi(ti);
			Matrix<double,2,2> coeA;
			coeA << - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c), sin(thet_c),
				- i_2h * cos(beta_c) * sin(thet_h - thet_1)
				- i_2h * sin(beta_c) * sin(thet_c) * cos(thet_h - thet_1),
				- i_2h * cos(thet_c) * cos(thet_h - thet_1);
			if(abs(coeA.determinant()) < 1.0E-12){
				cout<<"MLLOSE_NONEXISTENCE: thet_c="<<setprecision(20)<<thet_c<<endl;
				continue;
			}
			Matrix<double,2,1> coeB;
			coeB << - r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c),
				i_2h * a_1c * cos(beta_c) * cos(thet_c) * sin(thet_h - thet_1)
				+ i_2h * r_b2 * sin(beta_c) * cos(thet_c) * cos(thet_h - thet_1);
			Matrix<double,2,1> xy_d;
			xy_d = coeA.fullPivLu().solve(-coeB);
			Vector3d r_2_2;
			WHEE_1H2R(xy_d(0), xy_d(1), thet_1, thet_h, r_2_2);
			tempOfst<<setw(30)<<r_2_2(0)
				<<setw(30)<<r_2_2(1)<<setw(30)<<r_2_2(2)
				<<setw(30)<<r_2_2(0)
				<<setw(30)<<r_2_2(1)<<setw(30)<<r_2_2(2)<<endl;
		}
		default_random_engine dren;
		uniform_real_distribution<double> urdd(0.25, 0.75);
		long F_1 = 0;
		long F_2 = 0;
		long N = 9999;
		double deltThet_c = (thet_cH - thet_cL) / (N + 1.0);
		for(long tj = 1; tj <= N; tj ++){
			double thet_c = thet_cL + deltThet_c * tj;
			while(true){
				long flag_3 = 0;
				for(long ti = 0; ti < thet_ci.rows(); ti ++){
					if(abs(thet_c - thet_ci(ti)) < 1.0E-14){
						flag_3 = 1;
					}
				}
				if(flag_3 == 0){
					break;
				}
				else{
					thet_c = thet_cL + deltThet_c * (tj - urdd(dren));
				}
			}
			double thet_1 = i_1c * thet_c;
			double C_i1 = (i_c1 * i_c1 * i_2h * i_2h * (a_h2 * a_h2 - a_1c * a_1c) 
				- (i_2h * a_h2 - i_c1 * a_1c) * (i_2h * a_h2 - i_c1 * a_1c)) 
				* cos(beta_c) * cos(beta_c);
			double C_i2 = 2.0 * i_c1 * i_2h * (i_2h * a_h2 - i_c1 * a_1c)
				* a_h2 * cos(beta_c) * sin(beta_c);
			double C_i3 = (i_2h * a_h2 - i_c1 * a_1c) * (i_2h * a_h2 - i_c1 * a_1c);
			if(C_i1 * cos(thet_c) * cos(thet_c) + C_i2 * cos(thet_c) + C_i3 < 0.0){
				F_2 = 0;
				if(tj % 100 == 0){
					cout<<"MLLOSE_NONEXISTENCE: thet_c="<<setprecision(20)<<thet_c<<endl;
				}
				continue;
			}
			double thet_hs, thet_hm;
			SINGULAR_C2H(thet_c, thet_hs, thet_hm);
			double C_s21 = (- i_2h * a_h2 + i_c1 * a_1c) * sin(beta_c)
				- i_c1 * i_2h * a_h2 * cos(beta_c) * cos(thet_c);
			double C_s22 = (- i_2h * a_h2 + i_c1 * a_1c) * cos(beta_c) * sin(thet_c);
			double C_s23 = i_c1 * i_2h * a_1c * cos(beta_c) * cos(thet_c);
			double a2CC = atan2(C_s21, C_s22);
			Matrix<double,2,1> thet_h;
			thet_h(0) = thet_1 - a2CC 
				- asin(C_s23 / sqrt(C_s21 * C_s21 + C_s22 * C_s22));
			thet_h(1) = thet_1 - a2CC + PI 
				+ asin(C_s23 / sqrt(C_s21 * C_s21 + C_s22 * C_s22));
			while(thet_h(0) < thet_hs){
				thet_h(0) += 2.0 * PI;
			}
			while(thet_h(0) >= thet_hs + 2.0 * PI){
				thet_h(0) -= 2.0 * PI;
			}
			while(thet_h(1) < thet_hs){
				thet_h(1) += 2.0 * PI;
			}
			while(thet_h(1) >= thet_hs + 2.0 * PI){
				thet_h(1) -= 2.0 * PI;
			}
			Matrix<double,2,2> xy_d;
			FSME(thet_1, thet_h(0), xy_d(0,0), xy_d(1,0));
			FSME(thet_1, thet_h(1), xy_d(0,1), xy_d(1,1));
			if(tj >= 2 && F_1 == 1 && F_2 == 0){
				tempOfst<<"#"<<endl;
			}
			F_1 = F_2 = 1;
			Vector3d r_2_2;
			WHEE_1H2R(xy_d(0,0), xy_d(1,0), thet_1, thet_h(0), r_2_2);
			tempOfst<<setw(30)<<r_2_2(0)
				<<setw(30)<<r_2_2(1)<<setw(30)<<r_2_2(2);
			WHEE_1H2R(xy_d(0,1), xy_d(1,1), thet_1, thet_h(1), r_2_2);
			tempOfst<<setw(30)<<r_2_2(0)
				<<setw(30)<<r_2_2(1)<<setw(30)<<r_2_2(2)<<endl;
		}
	}
	return 1;
}

long DEHW::CILLOSE(double thet_cL, double thet_cH, 
	const MatrixXd &thet_cr, const MatrixXd &thet_hr, ofstream &tempOfst){
	OUTPUT_TIME("CILLOSE");
	double epsl_x = 1.0E-6;
	double epsl_t = 1.0E-4;
	if(modiTran == 0.0 && modiCent == 0.0){
		MLLOFE(2, thet_cL, thet_cH, tempOfst);
	}
	for(long ti = 0; ti < thet_ci.rows(); ti ++){
		long F_cicu = -1;
		for(long tj = 0; tj < thet_cu.rows(); tj ++){
			if(abs(thet_ci(ti) - thet_cu(tj)) <= 1.0E-14){
				F_cicu = tj;
				break;
			}
		}
		if(F_cicu != -1 && f_cu(F_cicu)){
			continue;
		}
		FUNC_DATA fuda;
		FUDA_SETUP(thet_ci(ti), fuda);
		fuda.thet_h = thet_hi(ti);
		fuda.C_s11 = - i_2h * cos(beta_c) * sin(fuda.thet_h - fuda.thet_1) 
			- i_2h * sin(beta_c) * sin(thet_ci(ti)) * cos(fuda.thet_h - fuda.thet_1);
		fuda.C_s12 = - i_2h * cos(thet_ci(ti)) * cos(fuda.thet_h - fuda.thet_1);
		fuda.C_s13 = i_2h * a_1c * cos(beta_c) * cos(thet_ci(ti)) * sin(fuda.thet_h - fuda.thet_1) 
			+ i_2h * r_b2 * sin(beta_c) * cos(thet_ci(ti)) * cos(fuda.thet_h - fuda.thet_1);
		long numbX;
		Matrix<double,Dynamic,Dynamic> limiX;
		if(F_cicu == -1){
			numbX = 2;
			limiX = MatrixXd::Zero(2,2);
			double x_du = (fuda.C_f22 * fuda.C_13 - fuda.C_f23 * fuda.C_12)
				/ (fuda.C_f21 * fuda.C_12 - fuda.C_f22 * fuda.C_11);
			double x_dL = min(- 10.0 * a_1c, x_du - 10.0 * a_1c);
			double x_dH = max(10.0 * a_1c, x_du + 10.0 * a_1c);
			limiX << x_dL, x_du - epsl_x,
				x_du + epsl_x, x_dH;
				
		}
		else{
			numbX = 1;
			limiX = MatrixXd::Zero(1,2);
			limiX << - 10.0 * a_1c, 10.0 * a_1c;
		}
		for(long tj = 0; tj < numbX; tj ++){
			EINM<FUNC_DATA> einm;
			einm.F = CILFOSE_X;
			einm.DF = PD_CILFOSE_X;
			einm.fuda = fuda;
			einm.searRang = FINTERVAL(limiX(tj,0), limiX(tj,1));
			einm.epsl_rX = 1.0E-9;
			einm.epsl_X = 1.0E-9;
			einm.epsl_f = 1.0E-9;
			einm.SOLVE();
			//
			for(list<pair<long, FINTERVAL>>::const_iterator iter_d = einm.done.begin(); 
				iter_d != einm.done.end(); iter_d ++){
				double x_d = filib::mid((*iter_d).second);
				cout<<"("<<ti<<"):"<<(*iter_d).second<<","
					<<CILFOSE_X(FINTERVAL(x_d), fuda)<<" - "<<(*iter_d).first<<endl;
				double y_d = - (fuda.C_11 * x_d + fuda.C_13) / fuda.C_12;
				Vector3d r_2_2;
				WHEE_1H2R(x_d, y_d, i_1c * thet_ci(ti), thet_hi(ti), r_2_2);
				tempOfst<<setw(30)<<r_2_2(0)
					<<setw(30)<<r_2_2(1)<<setw(30)<<r_2_2(2)<<endl;
			}
		}
	}
	for(long tj = 0; tj < thet_cr.rows(); tj ++){
		FUNC_DATA fuda;
		FUDA_SETUP(thet_cr(tj), fuda);
		for(long tk = 0; tk < thet_hr.cols() - 1; tk ++){
			if(thet_hr(tj, tk + 1) - thet_hr(tj, tk) > 2.0 * epsl_t){
				EINM<FUNC_DATA> einm;
				einm.F = CILFOSE_;
				einm.DF = PD_CILFOSE_;
				einm.fuda = fuda;
				einm.searRang = FINTERVAL(
					thet_hr(tj, tk) + epsl_t, thet_hr(tj, tk + 1) - epsl_t);
				einm.epsl_rX = 1.0E-9;
				einm.epsl_X = 1.0E-9;
				einm.epsl_f = 1.0E-9;
				einm.SOLVE();
				//
				for(list<pair<long, FINTERVAL>>::const_iterator iter_d = einm.done.begin(); 
					iter_d != einm.done.end(); iter_d ++){
					double thet_h = filib::mid((*iter_d).second);
					if((*iter_d).first != 1 || tj % 20 == 0){//if(tj % 20 == 0){
						cout<<"("<<tj<<","<<tk<<"):"<<(*iter_d).second<<","
							<<CILFOSE(FINTERVAL(thet_h), fuda)<<" - "<<(*iter_d).first<<endl;
					}
					double thet_1 = i_1c * thet_cr(tj);
					double x_d, y_d;
					FSME(thet_1, thet_h, x_d, y_d);
					Vector3d r_2_2;
					WHEE_1H2R(x_d, y_d, thet_1, thet_h, r_2_2);
					tempOfst<<setw(30)<<r_2_2(0)
						<<setw(30)<<r_2_2(1)<<setw(30)<<r_2_2(2)<<endl;
				}
			}
		}
	}
	OUTPUT_TIME("Done");
	return 1;
}

long DEHW::SINGULAR_H2C(double thet_h, double &thet_cmini, double &thet_cmaxi){
	double thet_cL = 0.0;
	double thet_cH = PI / 2.0;
	while(thet_cH - thet_cL > 1.0E-12){
		thet_cmaxi = (thet_cH + thet_cL) / 2.0;
		double thet_hs, thet_hm;
		SINGULAR_C2H(thet_cmaxi, thet_hs, thet_hm);
		if(thet_hs < thet_h){
			thet_cL = thet_cmaxi;
		}
		else{
			thet_cH = thet_cmaxi;
		}
	}
	thet_cL = 0.0;
	thet_cH = PI / 2.0;
	while(thet_cH - thet_cL > 1.0E-12){
		thet_cmini = (thet_cH + thet_cL) / 2.0;
		double thet_hs, thet_hm;
		SINGULAR_C2H(thet_cmini, thet_hs, thet_hm);
		thet_hs = thet_hs + 2.0 * PI;
		if(thet_hs < thet_h){
			thet_cL = thet_cmini;
		}
		else{
			thet_cH = thet_cmini;
		}
	}
	
	return 1;
}

long DEHW::SINGULAR_C2H(double thet_c, double &thet_hs, double &thet_hm){
//0 < thet_c < PI / 2.0
	double thet_1 = i_1c * thet_c;
	double C_m11 = - i_2h * cos(beta_c) * sin(thet_c);
	double C_m12 = i_c1 * i_2h * cos(beta_c) * cos(thet_c)
		+ i_2h * sin(beta_c);
	double C_m13 = i_c1 * cos(beta_c) * sin(thet_c);
	double a2CC = atan2(C_m11, C_m12);
	if(C_m13 > sqrt(C_m11 * C_m11 + C_m12 * C_m12)){
		thet_hs = thet_1 - a2CC - PI / 2.0;
		thet_hm = thet_hs;
	}
	else{
		thet_hs = thet_1 - PI - a2CC + asin(C_m13 / sqrt(C_m11 * C_m11 + C_m12 * C_m12));
		thet_hm = thet_1 - a2CC - asin(C_m13 / sqrt(C_m11 * C_m11 + C_m12 * C_m12));
	}
	
	return 1;
}

long DEHW::FSME(double thet_1, double thet_h, double &x_d, double &y_d){
	//
	double thet_c = i_c1 * thet_1;
	//
	Matrix<double,2,2> coefA;
	coefA << - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c), sin(thet_c),
		sin(beta_c) * cos(thet_c)
		+ i_2h * cos(beta_c) * cos(thet_h - thet_1) 
		- i_2h * sin(beta_c) * sin(thet_c) * sin(thet_h - thet_1),
		-sin(thet_c) - i_2h * cos(thet_c) * sin(thet_h -thet_1);
	Matrix<double,2,1> coefB;
	coefB << - r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c),
		+ r_b2 * sin(beta_c) * sin(thet_c) 
		+ i_2h * r_b2 * sin(beta_c) * cos(thet_c) * sin(thet_h - thet_1)
		- a_1c * sin(beta_c)
		- i_2h * a_1c * cos(beta_c) * cos(thet_c) * cos(thet_h - thet_1)
		+ i_2h * a_h2 * cos(beta_c) * cos(thet_c);
	Matrix<double,2,1> xy_d = coefA.fullPivLu().solve(-coefB);
	x_d = xy_d(0);
	y_d = xy_d(1);
	
	return 1;
}

long DEHW::PD_FSME(double thet_1, double thet_h, 
	double &x_d, double &y_d, Matrix<double,2,2> &Pxy_d){
	//
	double thet_c = i_c1 * thet_1;
	//
	Matrix<double,2,2> coefA;
	coefA << - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c), sin(thet_c),
		sin(beta_c) * cos(thet_c)
		+ i_2h * cos(beta_c) * cos(thet_h - thet_1) 
		- i_2h * sin(beta_c) * sin(thet_c) * sin(thet_h - thet_1),
		-sin(thet_c) - i_2h * cos(thet_c) * sin(thet_h -thet_1);
	Matrix<double,2,1> coefB;
	coefB << - r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c),
		+ r_b2 * sin(beta_c) * sin(thet_c) 
		+ i_2h * r_b2 * sin(beta_c) * cos(thet_c) * sin(thet_h - thet_1)
		- a_1c * sin(beta_c)
		- i_2h * a_1c * cos(beta_c) * cos(thet_c) * cos(thet_h - thet_1)
		+ i_2h * a_h2 * cos(beta_c) * cos(thet_c);
	Matrix<double,2,1> xy_d = coefA.fullPivLu().solve(-coefB);
	x_d = xy_d(0);
	y_d = xy_d(1);
	//derivative relative to thet_1
	Matrix<double,2,2> PcoefA;
	PcoefA << sin(beta_c) * sin(thet_c) * i_c1, cos(thet_c) * i_c1,
		- sin(beta_c) * sin(thet_c) * i_c1
		- i_2h * cos(beta_c) * sin(thet_h - thet_1) * -1.0
		- i_2h * sin(beta_c) * cos(thet_c) * i_c1 * sin(thet_h - thet_1)
		- i_2h * sin(beta_c) * sin(thet_c) * cos(thet_h - thet_1) * -1.0,
		- cos(thet_c) * i_c1 
		+ i_2h * sin(thet_c) * i_c1 * sin(thet_h -thet_1) 
		- i_2h * cos(thet_c) * cos(thet_h -thet_1) * -1.0;
	Matrix<double,2,1> PcoefB;
	PcoefB << - r_b2 * sin(beta_c) * cos(thet_c) * i_c1,
		r_b2 * sin(beta_c) * cos(thet_c) * i_c1
		- i_2h * r_b2 * sin(beta_c) * sin(thet_c) * i_c1 * sin(thet_h - thet_1)
		+ i_2h * r_b2 * sin(beta_c) * cos(thet_c) * cos(thet_h - thet_1) * -1.0
		+ i_2h * a_1c * cos(beta_c) * sin(thet_c) * i_c1 * cos(thet_h - thet_1)
		+ i_2h * a_1c * cos(beta_c) * cos(thet_c) * sin(thet_h - thet_1) * -1.0
		- i_2h * a_h2 * cos(beta_c) * sin(thet_c) * i_c1;
	Pxy_d.block(0,0,2,1) = coefA.fullPivLu().solve(- PcoefB - PcoefA * xy_d);
	//derivative relative to thet_h
	PcoefA << 0.0, 0.0,
		0.0 
		- i_2h * cos(beta_c) * sin(thet_h - thet_1) 
		- i_2h * sin(beta_c) * sin(thet_c) * cos(thet_h - thet_1),
		- 0.0 - i_2h * cos(thet_c) * cos(thet_h -thet_1);
	PcoefB << 0.0 + 0.0,
		+ 0.0 
		+ i_2h * r_b2 * sin(beta_c) * cos(thet_c) * cos(thet_h - thet_1)
		- 0.0
		+ i_2h * a_1c * cos(beta_c) * cos(thet_c) * sin(thet_h - thet_1)
		+ 0.0;
	Pxy_d.block(0,1,2,1) = coefA.fullPivLu().solve(- PcoefB - PcoefA * xy_d);
	return 1;
}

long DEHW::WORM_DC2R(double x_d, double y_d, double thet_c, Vector3d &r_1_1){
	double thet_1 = i_1c * thet_c;
	Vector3d xyz;
	xyz << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
	Matrix<double,3,3> rota;
	//R_oc,c
	rota << cos(thet_c), -sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	xyz = rota * xyz;
	//R_o1,oc
	rota << 1.0, 0.0, 0.0,
		0.0, 0.0, -1.0,
		0.0, 1.0, 0.0;
	xyz = rota * xyz;
	//T_o1,oc
	xyz(0) = xyz(0) + a_1c;
	//R_1,o1
	rota << cos(thet_1), sin(thet_1), 0.0,
		-sin(thet_1), cos(thet_1), 0.0,
		0.0, 0.0, 1.0;
	r_1_1 = rota * xyz;
	
	return 1;
}

long DEHW::WHEE_1H2R(double x_d, double y_d, double thet_1, double thet_h, 
	Vector3d &r_2_2){
	//
	double thet_c = i_c1 * thet_1;
	double thet_2 = i_2h * thet_h;
	//
	Vector3d xyz;
	Matrix<double,3,3> rota;
	WORM_DC2R(x_d, y_d, thet_c, xyz);
	//R_oh,h
	rota << cos(thet_h), -sin(thet_h), 0.0,
		sin(thet_h), cos(thet_h), 0.0,
		0.0, 0.0, 1.0;
	xyz = rota * xyz;
	//R_o2,oh
	rota << 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, -1.0, 0.0;
	xyz = rota * xyz;
	//T_o2,oh
	xyz(0) = xyz(0) - a_h2;
	//R_2,o2
	rota << cos(thet_2), sin(thet_2), 0.0,
		-sin(thet_2), cos(thet_2), 0.0,
		0.0, 0.0, 1.0;
	r_2_2 = rota * xyz;
	return 1;
}

long DEHW::PD_WHEE_1H2R(double x_d, double y_d, double thet_1, double thet_h, 
	Vector3d &r_2_2, Matrix<double,3,2> &Dr_2_2){
	//
	double thet_c = i_c1 * thet_1;
	double thet_2 = i_2h * thet_h;
	//relative to thet_1
	Matrix<double,2,2> Dxy_d;
	PD_FSME(thet_1, thet_h, x_d, y_d, Dxy_d);
	double Dx_d = Dxy_d(0,0);
	double Dy_d = Dxy_d(1,0);
	Vector3d r_c_c, Dr_c_c;
	r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
	Dr_c_c << - Dx_d, - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
	Matrix<double,3,3> R_oc_c, DR_oc_c;
	R_oc_c << cos(thet_c), - sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	DR_oc_c << - i_c1 * sin(thet_c), - i_c1 * cos(thet_c), 0.0,
		i_c1 * cos(thet_c), - i_c1 * sin(thet_c), 0.0,
		0.0, 0.0, 0.0;
	Vector3d r_c_oc = R_oc_c * r_c_c;
	Vector3d Dr_c_oc = DR_oc_c * r_c_c + R_oc_c * Dr_c_c;
	Matrix<double,3,3> R_o1_oc;
	Vector3d T_o1_oc;
	R_o1_oc << 1.0, 0.0, 0.0,
		0.0, 0.0, -1.0,
		0.0, 1.0, 0.0;
	T_o1_oc << a_1c, 0.0, 0.0;
	Vector3d r_c_o1 = R_o1_oc * r_c_oc + T_o1_oc;
	Vector3d Dr_c_o1 = R_o1_oc * Dr_c_oc;
	Matrix<double,3,3> R_1_o1, DR_1_o1;
	R_1_o1 << cos(thet_1), sin(thet_1), 0.0,
		-sin(thet_1), cos(thet_1), 0.0,
		0.0, 0.0, 1.0;
	DR_1_o1 << -sin(thet_1), cos(thet_1), 0.0,
		-cos(thet_1), -sin(thet_1), 0.0,
		0.0, 0.0, 0.0;
	Vector3d r_1_1 = R_1_o1 * r_c_o1;
	Vector3d Dr_1_1 = DR_1_o1 * r_c_o1 + R_1_o1 * Dr_c_o1;
	Matrix<double,3,3> R_oh_h;
	R_oh_h << cos(thet_h), - sin(thet_h), 0.0,
		sin(thet_h), cos(thet_h), 0.0,
		0.0, 0.0, 1.0;
	Vector3d r_h_oh = R_oh_h * r_1_1;
	Vector3d Dr_h_oh = R_oh_h * Dr_1_1;
	Matrix<double,3,3> R_o2_oh;
	Vector3d T_o2_oh;
	R_o2_oh << 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, -1.0, 0.0;
	T_o2_oh << - a_h2, 0.0, 0.0;
	Vector3d r_h_o2 = R_o2_oh * r_h_oh + T_o2_oh;
	Vector3d Dr_h_o2 = R_o2_oh * Dr_h_oh;
	Matrix<double,3,3> R_2_o2;
	R_2_o2 << cos(thet_2), sin(thet_2), 0.0,
		-sin(thet_2), cos(thet_2), 0.0,
		0.0, 0.0, 1.0;
	r_2_2 = R_2_o2 * r_h_o2;
	Dr_2_2.block(0,0,3,1) = R_2_o2 * Dr_h_o2;
	//relative to thet_h
	Dx_d = Dxy_d(0,1);
	Dy_d = Dxy_d(1,1);
	Dr_c_c << - Dx_d, - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
	Dr_c_oc = R_oc_c * Dr_c_c;
	Dr_c_o1 = R_o1_oc * Dr_c_oc;
	Dr_1_1 = R_1_o1 * Dr_c_o1;
	Matrix<double,3,3> DR_oh_h;
	DR_oh_h << -sin(thet_h), -cos(thet_h), 0.0,
		cos(thet_h), -sin(thet_h), 0.0,
		0.0, 0.0, 0.0;
	Dr_h_oh = DR_oh_h * r_1_1 + R_oh_h * Dr_1_1;
	Dr_h_o2 = R_o2_oh * Dr_h_oh;
	Matrix<double,3,3> DR_2_o2;
	DR_2_o2 << -i_2h * sin(thet_2), i_2h * cos(thet_2), 0.0,
		-i_2h * cos(thet_2), -i_2h * sin(thet_2), 0.0,
		0.0, 0.0, 0.0;
	Dr_2_2.block(0,1,3,1) = DR_2_o2 * r_h_o2 + R_2_o2 * Dr_h_o2;
	return 1;
}

long DEHW::FUDA_SETUP(const double thet_c, FUNC_DATA &fuda){
	fuda.thet_1 = i_1c * thet_c;
	fuda.C_11 = - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c);
	fuda.C_12 = sin(thet_c);
	fuda.C_13 = - r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c);
	fuda.C_m11 = - i_2h * cos(beta_c) * sin(thet_c);
	fuda.C_m12 = i_c1 * i_2h * cos(beta_c) * cos(thet_c) 
		+ i_2h * sin(beta_c);
	fuda.C_m13 = i_c1 * cos(beta_c) * sin(thet_c);
	fuda.C_a11 = sqrt(fuda.C_m11 * fuda.C_m11 + fuda.C_m12 * fuda.C_m12);
	fuda.thet_a1 = - fuda.thet_1 + atan2(fuda.C_m11, fuda.C_m12);
	fuda.C_m21 = i_2h * a_1c * cos(beta_c) * cos(thet_c) * sin(thet_c);
	fuda.C_m22 = - i_2h * a_1c * sin(beta_c) * cos(thet_c);
	fuda.C_m23 = - i_2h * a_h2 * cos(beta_c) * cos(thet_c) * sin(thet_c);
	fuda.C_a21 = sqrt(fuda.C_m21 * fuda.C_m21 + fuda.C_m22 * fuda.C_m22);
	fuda.thet_a2 = - fuda.thet_1 + atan2(fuda.C_m21, fuda.C_m22);
	fuda.C_m31 = 
		- i_2h * a_1c * cos(beta_c) * sin(beta_c) * sin(thet_c) * sin(thet_c)
		+ i_2h * r_b2 * cos(beta_c) * sin(beta_c) * sin(thet_c)
		+ i_c1 * i_2h * a_1c * cos(beta_c) * cos(beta_c) * cos(thet_c);
	fuda.C_m32 = - i_2h * r_b2 * sin(beta_c) * sin(beta_c)
		+ i_2h * a_1c * sin(beta_c) * sin(beta_c) * sin(thet_c)
		- i_c1 * i_2h * r_b2 * cos(beta_c) * sin(beta_c) * cos(thet_c);
	fuda.C_m33 = 
		- i_2h * a_h2 * cos(beta_c) * sin(beta_c) * cos(thet_c) * cos(thet_c)
		- i_c1 * r_b2 * cos(beta_c) * sin(beta_c) * sin(thet_c)
		+ i_c1 * a_1c * cos(beta_c) * sin(beta_c)
		- i_c1 * i_2h * a_h2 * cos(beta_c) * cos(beta_c) * cos(thet_c);
	fuda.C_a31 = sqrt(fuda.C_m31 * fuda.C_m31 + fuda.C_m32 * fuda.C_m32);
	fuda.thet_a3 = - fuda.thet_1 + atan2(fuda.C_m31, fuda.C_m32);
	fuda.C_f21 = cos(beta_c) * cos(thet_c) * sin(thet_c);
	fuda.C_f22 = cos(beta_c) * sin(beta_c) * cos(thet_c) * cos(thet_c)
		+ 2.0 * i_c1 * cos(beta_c) * cos(beta_c) * cos(thet_c)
		- i_c1 * i_c1 * cos(beta_c) * sin(beta_c);
	fuda.C_f23 = r_b2 * cos(beta_c) * sin(thet_c) * sin(thet_c)
		+ i_c1 * i_c1 * r_b2 * cos(beta_c)
		- a_1c * cos(beta_c) * sin(thet_c);
	fuda.C_a41 = - fuda.C_f21 * fuda.C_m21 - fuda.C_f22 * fuda.C_m31 + fuda.C_f23 * fuda.C_m11;
	fuda.C_a42 = - fuda.C_f21 * fuda.C_m22 - fuda.C_f22 * fuda.C_m32 + fuda.C_f23 * fuda.C_m12;
	fuda.C_a51 = sqrt(fuda.C_a41 * fuda.C_a41 + fuda.C_a42 * fuda.C_a42);
	fuda.C_a52 = - fuda.C_f21 * fuda.C_m23 - fuda.C_f22 * fuda.C_m33 + fuda.C_f23 * fuda.C_m13;
	fuda.thet_a5 = - fuda.thet_1 + atan2(fuda.C_a41, fuda.C_a42);
	fuda.C_a61 = - (sin(beta_c) * cos(thet_c) + i_c1 * cos(beta_c)) 
		* (sin(beta_c) * cos(thet_c) + i_c1 * cos(beta_c));
	fuda.C_a62 = - sin(thet_c) * sin(thet_c);
	fuda.C_a63 = (sin(beta_c) * cos(thet_c) 
		+ i_c1 * cos(beta_c)) * sin(thet_c);
	fuda.C_b11 = sqrt(i_2h * i_2h * sin(beta_c) * sin(beta_c)
		+ i_2h * i_2h * cos(beta_c) * cos(beta_c) * sin(thet_c) * sin(thet_c));
	fuda.C_b12 = - cos(beta_c) * cos(thet_c);
	fuda.thet_b1 = - fuda.thet_1 + atan2(i_2h * sin(beta_c), 
		i_2h * cos(beta_c) * sin(thet_c));
	fuda.C_b21 = - i_2h * r_b2 + i_2h * a_1c * sin(thet_c);
	fuda.C_b22 = - i_2h * a_h2 * sin(thet_c);
	fuda.C_b31 = sqrt(i_2h * i_2h * sin(beta_c) * sin(beta_c)
		+ i_2h * i_2h * cos(beta_c) * cos(beta_c) * sin(thet_c) * sin(thet_c));
	fuda.C_b32 = cos(beta_c) * cos(thet_c);
	fuda.thet_b3 = - fuda.thet_1 + atan2(- i_2h * sin(beta_c), 
		- i_2h * cos(beta_c) * sin(thet_c));
	fuda.C_b41 = sqrt(
		i_2h * i_2h * a_1c * a_1c 
		* sin(beta_c) * sin(beta_c) * cos(thet_c) * cos(thet_c)
		+ i_2h * i_2h * r_b2 * r_b2 
		* cos(beta_c) * cos(beta_c) * cos(thet_c) * cos(thet_c));
	fuda.C_b42 = r_b2 * cos(beta_c) * sin(thet_c) - a_1c * cos(beta_c) 
		- i_2h * a_h2 * sin(beta_c) * cos(thet_c);
	fuda.thet_b4 = - fuda.thet_1 + atan2(i_2h * a_1c * sin(beta_c) * cos(thet_c), 
		i_2h * r_b2 * cos(beta_c) * cos(thet_c));
	fuda.C_b51 = sqrt(i_2h * i_2h * cos(beta_c) * cos(beta_c)
		+ i_2h * i_2h * sin(beta_c) * sin(beta_c) * sin(thet_c) * sin(thet_c));
	fuda.C_b52 = - sin(beta_c) * cos(thet_c);
	fuda.thet_b5 = - fuda.thet_1 + atan2(- i_2h * cos(beta_c), 
		i_2h * sin(beta_c) * sin(thet_c));
	fuda.C_b61 = - i_2h * cos(thet_c);
	fuda.C_b62 = - sin(thet_c);
	fuda.C_b71 = - cos(thet_c);
	fuda.C_b72 = - sin(thet_c);
	fuda.C_b81 = sqrt(sin(beta_c) * sin(beta_c) * sin(thet_c) * sin(thet_c)
		+ cos(beta_c) * cos(beta_c));
	fuda.thet_b81 = - fuda.thet_1 + atan2(sin(beta_c) * sin(thet_c), cos(beta_c));
	fuda.thet_b82 = - fuda.thet_1 + atan2(- cos(beta_c), sin(beta_c) * sin(thet_c));
	fuda.C_b83 = - sin(beta_c) * cos(thet_c);
	fuda.C_c11 = cos(thet_c);
	fuda.C_c12 = - i_2h * sin(thet_c);
	fuda.C_c21 = sqrt(cos(beta_c) * cos(beta_c) 
		+ sin(beta_c) * sin(beta_c) * sin(thet_c) * sin(thet_c));
	fuda.C_c22 = - i_2h * sin(beta_c) * cos(thet_c);
	fuda.thet_c2 = - fuda.thet_1 + atan2(cos(beta_c), - sin(beta_c) * sin(thet_c));
	fuda.C_c31 = r_b2 * sin(thet_c) - a_1c;
	fuda.C_c32 = i_2h * r_b2 * cos(thet_c);
	fuda.C_c41 = - cos(thet_c);
	fuda.C_c42 = sqrt(sin(beta_c) * sin(beta_c) * sin(thet_c) * sin(thet_c)
		+ cos(beta_c) * cos(beta_c));
	fuda.thet_c4 = - fuda.thet_1 + atan2(sin(beta_c) * sin(thet_c), cos(beta_c));
	fuda.C_c43 = - r_b2 * sin(thet_c) + a_1c;
	fuda.C_c51 = i_2h * cos(thet_c);
	fuda.C_c52 = sqrt(i_2h * i_2h 
		* sin(beta_c) * sin(beta_c) * sin(thet_c) * sin(thet_c)
		+ i_2h * i_2h * cos(beta_c) * cos(beta_c));
	fuda.thet_c5 = - fuda.thet_1 + atan2(- i_2h * sin(beta_c) * sin(thet_c),
		- i_2h * cos(beta_c));
	fuda.C_c53 = i_2h * r_b2 * sin(thet_c) - i_2h * a_1c;
	fuda.C_c54 = i_2h * a_h2;
	fuda.C_c61 = i_2h * cos(beta_c) * cos(thet_c);
	fuda.C_s21 = (- i_2h * a_h2 + i_c1 * a_1c) * sin(beta_c)
		- i_c1 * i_2h * a_h2 * cos(beta_c) * cos(thet_c);
	fuda.C_s22 = (- i_2h * a_h2 + i_c1 * a_1c) 
		* cos(beta_c) * sin(thet_c);
	fuda.C_s23 = i_c1 * i_2h * a_1c * cos(beta_c) * cos(thet_c);
	fuda.C_c62 = sqrt(fuda.C_s21 * fuda.C_s21 + fuda.C_s22 * fuda.C_s22);
	fuda.thet_c6 = - fuda.thet_1 + atan2(fuda.C_s21, fuda.C_s22);
	return 1;
}

FINTERVAL CILFOSE_(const FINTERVAL &THET_H, const FUNC_DATA &fuda){
	long segmNumb = 10;
	FINTERVAL resu = FINTERVAL::EMPTY();
	for(long ti = 0; ti < segmNumb; ti ++){
		resu = hull(resu, CILFOSE(FINTERVAL(
			THET_H.inf() + (THET_H.sup() - THET_H.inf()) / segmNumb * ti, 
			THET_H.inf() + (THET_H.sup() - THET_H.inf()) / segmNumb * (ti + 1)), 
			fuda));
	}
	return resu;
}

FINTERVAL CILFOSE(const FINTERVAL &THET_H, const FUNC_DATA &fuda){
	FINTERVAL T_1 = fuda.C_a11 * sin(THET_H + fuda.thet_a1) + fuda.C_m13;//
	FINTERVAL T_2 = - fuda.C_a21 * sin(THET_H + fuda.thet_a2) - fuda.C_m23;
	FINTERVAL T_3 = T_2 / T_1;
	FINTERVAL T_4 = - fuda.C_a31 * sin(THET_H + fuda.thet_a3) - fuda.C_m33;
	FINTERVAL T_5 = T_4 / T_1;
	FINTERVAL T_6 = fuda.C_a51 * sin(THET_H + fuda.thet_a5) + fuda.C_a52;//
	FINTERVAL T_7 = T_6 / T_1;
	FINTERVAL T_8 = fuda.C_a61 / T_7;
	FINTERVAL T_9 = fuda.C_a62 / T_7;
	FINTERVAL T_10 = fuda.C_a63 / T_7;
	FINTERVAL T_11 = fuda.C_b11 * sin(THET_H + fuda.thet_b1) + fuda.C_b12;
	FINTERVAL T_12 = fuda.C_b21 * cos(THET_H - fuda.thet_1) + fuda.C_b22;
	FINTERVAL T_13 = T_11 * T_5 + T_12;
	FINTERVAL T_14 = fuda.C_b31 * sin(THET_H + fuda.thet_b3) + fuda.C_b32;
	FINTERVAL T_15 = fuda.C_b41 * sin(THET_H + fuda.thet_b4) + fuda.C_b42;
	FINTERVAL T_16 = T_14 * T_3 + T_15;
	FINTERVAL T_17 = fuda.C_b51 * sin(THET_H + fuda.thet_b5) + fuda.C_b52;
	FINTERVAL T_18 = fuda.C_b61 * sin(THET_H - fuda.thet_1) + fuda.C_b62;
	FINTERVAL T_19 = T_8 * T_13 + T_10 * T_16 + T_17;
	FINTERVAL T_20 = T_10 * T_13 + T_9 * T_16 - T_18;
	FINTERVAL T_21 = fuda.C_b71 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_22 = fuda.C_b71 * sin(THET_H - fuda.thet_1);
	FINTERVAL T_23 = fuda.C_b81 * sin(THET_H + fuda.thet_b81);
	FINTERVAL T_24 = fuda.C_b81 * sin(THET_H + fuda.thet_b82);
	FINTERVAL T_25 = T_19 * T_21 + T_20 * T_23;
	FINTERVAL T_26 = T_19 * T_22 + T_20 * T_24;
	FINTERVAL T_27 = fuda.C_b72 * T_19 + fuda.C_b83 * T_20;
	FINTERVAL T_28 = fuda.C_c11 * sin(THET_H - fuda.thet_1) + fuda.C_c12;
	FINTERVAL T_29 = fuda.C_c21 * sin(THET_H + fuda.thet_c2) + fuda.C_c22;
	FINTERVAL T_30 = fuda.C_c31 * sin(THET_H - fuda.thet_1) + fuda.C_c32;
	FINTERVAL T_31 = T_28 * T_3 + T_29 * T_5 + T_30;
	FINTERVAL T_32 = fuda.C_c41 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_33 = fuda.C_c42 * sin(THET_H + fuda.thet_c4);
	FINTERVAL T_34 = fuda.C_c43 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_35 = T_32 * T_3 + T_33 * T_5 + T_34;
	FINTERVAL T_36 = fuda.C_c51 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_37 = fuda.C_c52 * sin(THET_H + fuda.thet_c5);
	FINTERVAL T_38 = fuda.C_c53 * cos(THET_H - fuda.thet_1) + fuda.C_c54;
	FINTERVAL T_39 = T_36 * T_3 + T_37 * T_5 + T_38;
	FINTERVAL T_40 = fuda.C_c61 * (fuda.C_c62 * sin(THET_H + fuda.thet_c6) + fuda.C_s23);
	FINTERVAL T_41 = T_40 / T_1;
	FINTERVAL T_42 = T_25 * T_31 + T_26 * T_35 + T_27 * T_39 + T_41;
	return T_42;
}

FINTERVAL PD_CILFOSE_(const FINTERVAL &THET_H, const FUNC_DATA &fuda){
	long segmNumb = 10;
	FINTERVAL resu = FINTERVAL::EMPTY();
	for(long ti = 0; ti < segmNumb; ti ++){
		resu = hull(resu, PD_CILFOSE(FINTERVAL(
			THET_H.inf() + (THET_H.sup() - THET_H.inf()) / segmNumb * ti, 
			THET_H.inf() + (THET_H.sup() - THET_H.inf()) / segmNumb * (ti + 1)), 
			fuda));
	}
	return resu;
}

FINTERVAL PD_CILFOSE(const FINTERVAL &THET_H, const FUNC_DATA &fuda){
	FINTERVAL T_1 = fuda.C_a11 * sin(THET_H + fuda.thet_a1) + fuda.C_m13;
	FINTERVAL T_2 = - fuda.C_a21 * sin(THET_H + fuda.thet_a2) - fuda.C_m23;
	FINTERVAL T_3 = T_2 / T_1;
	FINTERVAL T_4 = - fuda.C_a21 * cos(THET_H + fuda.thet_a2);
	FINTERVAL T_5 = fuda.C_a11 * cos(THET_H + fuda.thet_a1);
	FINTERVAL T_6 = power(T_1, 2);
	FINTERVAL T_7 = (T_4 * T_1 - T_2 * T_5) / T_6;
	FINTERVAL T_8 = - fuda.C_a31 * sin(THET_H + fuda.thet_a3) - fuda.C_m33;
	FINTERVAL T_9 = T_8 / T_1;
	FINTERVAL T_10 = - fuda.C_a31 * cos(THET_H + fuda.thet_a3);
	FINTERVAL T_11 = (T_10 * T_1 - T_8 * T_5) / T_6;
	FINTERVAL T_12 = fuda.C_a51 * sin(THET_H + fuda.thet_a5) + fuda.C_a52;
	FINTERVAL T_13 = T_12 / T_1;
	FINTERVAL T_14 = fuda.C_a51 * cos(THET_H + fuda.thet_a5);
	FINTERVAL T_15 = (T_14 * T_1 - T_12 * T_5) / T_6;
	FINTERVAL T_16 = fuda.C_a61 / T_13;
	FINTERVAL T_17 = power(T_13, 2);
	FINTERVAL T_18 = - fuda.C_a61 * T_15 / T_17;
	FINTERVAL T_19 = fuda.C_a62 / T_13;
	FINTERVAL T_20 = - fuda.C_a62 * T_15 / T_17;
	FINTERVAL T_21 = fuda.C_a63 / T_13;
	FINTERVAL T_22 = - fuda.C_a63 * T_15 / T_17;
	FINTERVAL T_23 = fuda.C_b11 * sin(THET_H + fuda.thet_b1) + fuda.C_b12;
	FINTERVAL T_24 = fuda.C_b21 * cos(THET_H - fuda.thet_1) + fuda.C_b22;
	FINTERVAL T_25 = T_23 * T_9 + T_24;
	FINTERVAL T_26 = fuda.C_b11 * cos(THET_H + fuda.thet_b1);
	FINTERVAL T_27 = - fuda.C_b21 * sin(THET_H - fuda.thet_1);
	FINTERVAL T_28 = T_26 * T_9 + T_23 * T_11 + T_27;
	FINTERVAL T_29 = fuda.C_b31 * sin(THET_H + fuda.thet_b3) + fuda.C_b32;
	FINTERVAL T_30 = fuda.C_b41 * sin(THET_H + fuda.thet_b4) + fuda.C_b42;
	FINTERVAL T_31 = T_29 * T_3 + T_30;
	FINTERVAL T_32 = fuda.C_b31 * cos(THET_H + fuda.thet_b3);
	FINTERVAL T_33 = fuda.C_b41 * cos(THET_H + fuda.thet_b4);
	FINTERVAL T_34 = T_32 * T_3 + T_29 * T_7 + T_33;
	FINTERVAL T_35 = fuda.C_b51 * sin(THET_H + fuda.thet_b5) + fuda.C_b52;
	FINTERVAL T_36 = fuda.C_b51 * cos(THET_H + fuda.thet_b5);
	FINTERVAL T_37 = fuda.C_b61 * sin(THET_H - fuda.thet_1) + fuda.C_b62;
	FINTERVAL T_38 = fuda.C_b61 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_39 = fuda.C_b71 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_40 = fuda.C_b71 * sin(THET_H - fuda.thet_1);
	FINTERVAL T_41 = fuda.C_b81 * sin(THET_H + fuda.thet_b81);
	FINTERVAL T_42 = fuda.C_b81 * sin(THET_H + fuda.thet_b82);
	FINTERVAL T_43 = fuda.C_b81 * cos(THET_H + fuda.thet_b81);
	FINTERVAL T_44 = fuda.C_b81 * cos(THET_H + fuda.thet_b82);
	FINTERVAL T_45 = T_16 * T_25 + T_21 * T_31 + T_35;
	FINTERVAL T_46 = T_21 * T_25 + T_19 * T_31 - T_37;
	FINTERVAL T_47 = T_18 * T_25 + T_16 * T_28 + T_22 * T_31 + T_21 * T_34 + T_36;
	FINTERVAL T_48 = T_22 * T_25 + T_21 * T_28 + T_20 * T_31 + T_19 * T_34 - T_38;
	FINTERVAL T_49 = T_45 * T_39 + T_46 * T_41;
	FINTERVAL T_50 = T_45 * T_40 + T_46 * T_42;
	FINTERVAL T_51 = fuda.C_b72 * T_45 + fuda.C_b83 * T_46;
	FINTERVAL T_52 = T_47 * T_39 - T_45 * T_40 + T_48 * T_41 + T_46 * T_43;
	FINTERVAL T_53 = T_47 * T_40 + T_45 * T_39 + T_48 * T_42 + T_46 * T_44;
	FINTERVAL T_54 = fuda.C_b72 * T_47 + fuda.C_b83 * T_48;
	FINTERVAL T_55 = fuda.C_c11 * sin(THET_H - fuda.thet_1) + fuda.C_c12;
	FINTERVAL T_56 = fuda.C_c21 * sin(THET_H + fuda.thet_c2) + fuda.C_c22;
	FINTERVAL T_57 = fuda.C_c31 * sin(THET_H - fuda.thet_1) + fuda.C_c32;
	FINTERVAL T_58 = T_55 * T_3 + T_56 * T_9 + T_57;
	FINTERVAL T_59 = fuda.C_c41 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_60 = fuda.C_c42 * sin(THET_H + fuda.thet_c4);
	FINTERVAL T_61 = fuda.C_c43 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_62 = T_59 * T_3 + T_60 * T_9 + T_61;
	FINTERVAL T_63 = fuda.C_c51 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_64 = fuda.C_c52 * sin(THET_H + fuda.thet_c5);
	FINTERVAL T_65 = fuda.C_c53 * cos(THET_H - fuda.thet_1) + fuda.C_c54;
	FINTERVAL T_66 = T_63 * T_3 + T_64 * T_9 + T_65;
	FINTERVAL T_67 = fuda.C_c11 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_68 = fuda.C_c21 * cos(THET_H + fuda.thet_c2);
	FINTERVAL T_69 = fuda.C_c31 * cos(THET_H - fuda.thet_1);
	FINTERVAL T_70 = T_67 * T_3 + T_55 * T_7 + T_68 * T_9 + T_56 * T_11 + T_69;
	FINTERVAL T_71 = - fuda.C_c41 * sin(THET_H - fuda.thet_1);
	FINTERVAL T_72 = fuda.C_c42 * cos(THET_H + fuda.thet_c4);
	FINTERVAL T_73 = - fuda.C_c43 * sin(THET_H - fuda.thet_1);
	FINTERVAL T_74 = T_71 * T_3 + T_59 * T_7 + T_72 * T_9 + T_60 * T_11 + T_73;
	FINTERVAL T_75 = - fuda.C_c51 * sin(THET_H - fuda.thet_1);
	FINTERVAL T_76 = fuda.C_c52 * cos(THET_H + fuda.thet_c5);
	FINTERVAL T_77 = - fuda.C_c53 * sin(THET_H - fuda.thet_1);
	FINTERVAL T_78 = T_75 * T_3 + T_63 * T_7 + T_76 * T_9 + T_64 * T_11 + T_77;
	FINTERVAL T_79 = fuda.C_c61 * (fuda.C_c62 * sin(THET_H + fuda.thet_c6) + fuda.C_s23);
	FINTERVAL T_80 = T_79 / T_1;
	FINTERVAL T_81 = fuda.C_c61 * fuda.C_c62 * cos(THET_H + fuda.thet_c6);
	FINTERVAL T_82 = (T_81 * T_1 - T_79 * T_5) / T_6;
	FINTERVAL T_83 = T_49 * T_58 + T_50 * T_62 + T_51 * T_66 + T_80;
	FINTERVAL T_84 = T_52 * T_58 + T_53 * T_62 + T_54 * T_66 
		+ T_49 * T_70 + T_50 * T_74 + T_51 * T_78 + T_82;
	return T_84;
}

FINTERVAL CILFOSE_X(const FINTERVAL &X, const FUNC_DATA &fuda){
	FINTERVAL T_3 = X;
	FINTERVAL T_5 = - (fuda.C_11 * T_3 + fuda.C_13) / fuda.C_12;
	FINTERVAL T_7 = fuda.C_f21 * T_3 + fuda.C_f22 * T_5 + fuda.C_f23;
	FINTERVAL T_8 = fuda.C_a61 / T_7;
	FINTERVAL T_9 = fuda.C_a62 / T_7;
	FINTERVAL T_10 = fuda.C_a63 / T_7;
	double t_11 = fuda.C_b11 * sin(fuda.thet_h + fuda.thet_b1) + fuda.C_b12;
	double t_12 = fuda.C_b21 * cos(fuda.thet_h - fuda.thet_1) + fuda.C_b22;
	FINTERVAL T_13 = t_11 * T_5 + t_12;
	double t_14 = fuda.C_b31 * sin(fuda.thet_h + fuda.thet_b3) + fuda.C_b32;
	double t_15 = fuda.C_b41 * sin(fuda.thet_h + fuda.thet_b4) + fuda.C_b42;
	FINTERVAL T_16 = t_14 * T_3 + t_15;
	double t_17 = fuda.C_b51 * sin(fuda.thet_h + fuda.thet_b5) + fuda.C_b52;
	double t_18 = fuda.C_b61 * sin(fuda.thet_h - fuda.thet_1) + fuda.C_b62;
	FINTERVAL T_19 = T_8 * T_13 + T_10 * T_16 + t_17;
	FINTERVAL T_20 = T_10 * T_13 + T_9 * T_16 - t_18;
	double t_21 = fuda.C_b71 * cos(fuda.thet_h - fuda.thet_1);
	double t_22 = fuda.C_b71 * sin(fuda.thet_h - fuda.thet_1);
	double t_23 = fuda.C_b81 * sin(fuda.thet_h + fuda.thet_b81);
	double t_24 = fuda.C_b81 * sin(fuda.thet_h + fuda.thet_b82);
	FINTERVAL T_25 = T_19 * t_21 + T_20 * t_23;
	FINTERVAL T_26 = T_19 * t_22 + T_20 * t_24;
	FINTERVAL T_27 = fuda.C_b72 * T_19 + fuda.C_b83 * T_20;
	double t_28 = fuda.C_c11 * sin(fuda.thet_h - fuda.thet_1) + fuda.C_c12;
	double t_29 = fuda.C_c21 * sin(fuda.thet_h + fuda.thet_c2) + fuda.C_c22;
	double t_30 = fuda.C_c31 * sin(fuda.thet_h - fuda.thet_1) + fuda.C_c32;
	FINTERVAL T_31 = t_28 * T_3 + t_29 * T_5 + t_30;
	double t_32 = fuda.C_c41 * cos(fuda.thet_h - fuda.thet_1);
	double t_33 = fuda.C_c42 * sin(fuda.thet_h + fuda.thet_c4);
	double t_34 = fuda.C_c43 * cos(fuda.thet_h - fuda.thet_1);
	FINTERVAL T_35 = t_32 * T_3 + t_33 * T_5 + t_34;
	double t_36 = fuda.C_c51 * cos(fuda.thet_h - fuda.thet_1);
	double t_37 = fuda.C_c52 * sin(fuda.thet_h + fuda.thet_c5);
	double t_38 = fuda.C_c53 * cos(fuda.thet_h - fuda.thet_1) + fuda.C_c54;
	FINTERVAL T_39 = t_36 * T_3 + t_37 * T_5 + t_38;
	FINTERVAL T_41 = fuda.C_s11 * T_3 + fuda.C_s12 * T_5 + fuda.C_s13;
	FINTERVAL T_42 = T_25 * T_31 + T_26 * T_35 + T_27 * T_39 + T_41;
	return T_42;
}

FINTERVAL PD_CILFOSE_X(const FINTERVAL &X, const FUNC_DATA &fuda){
	FINTERVAL T_3 = X;
	double t_7 = 1.0;
	FINTERVAL T_9 = - (fuda.C_11 * T_3 + fuda.C_13) / fuda.C_12;
	double t_11 = - fuda.C_11 / fuda.C_12;
	FINTERVAL T_13 = fuda.C_f21 * T_3 + fuda.C_f22 * T_9 + fuda.C_f23;
	double t_15 = fuda.C_f21 + fuda.C_f22 * t_11;
	FINTERVAL T_16 = fuda.C_a61 / T_13;
	FINTERVAL T_17 = power(T_13, 2);
	FINTERVAL T_18 = - fuda.C_a61 * t_15 / T_17;
	FINTERVAL T_19 = fuda.C_a62 / T_13;
	FINTERVAL T_20 = - fuda.C_a62 * t_15 / T_17;
	FINTERVAL T_21 = fuda.C_a63 / T_13;
	FINTERVAL T_22 = - fuda.C_a63 * t_15 / T_17;
	double t_23 = fuda.C_b11 * sin(fuda.thet_h + fuda.thet_b1) + fuda.C_b12;
	double t_24 = fuda.C_b21 * cos(fuda.thet_h - fuda.thet_1) + fuda.C_b22;
	FINTERVAL T_25 = t_23 * T_9 + t_24;
	double t_28 = t_23 * t_11;
	double t_29 = fuda.C_b31 * sin(fuda.thet_h + fuda.thet_b3) + fuda.C_b32;
	double t_30 = fuda.C_b41 * sin(fuda.thet_h + fuda.thet_b4) + fuda.C_b42;
	FINTERVAL T_31 = t_29 * T_3 + t_30;
	double t_34 = t_29 * t_7;
	double t_35 = fuda.C_b51 * sin(fuda.thet_h + fuda.thet_b5) + fuda.C_b52;
	double t_36 = 0.0;
	double t_37 = fuda.C_b61 * sin(fuda.thet_h - fuda.thet_1) + fuda.C_b62;
	double t_38 = 0.0;
	double t_39 = fuda.C_b71 * cos(fuda.thet_h - fuda.thet_1);
	double t_40 = fuda.C_b71 * sin(fuda.thet_h - fuda.thet_1);
	double t_41 = fuda.C_b81 * sin(fuda.thet_h + fuda.thet_b81);
	double t_42 = fuda.C_b81 * sin(fuda.thet_h + fuda.thet_b82);
	double t_43 = 0.0;
	double t_44 = 0.0;
	FINTERVAL T_45 = T_16 * T_25 + T_21 * T_31 + t_35;
	FINTERVAL T_46 = T_21 * T_25 + T_19 * T_31 - t_37;
	FINTERVAL T_47 = T_18 * T_25 + T_16 * t_28 + T_22 * T_31 + T_21 * t_34;
	FINTERVAL T_48 = T_22 * T_25 + T_21 * t_28 + T_20 * T_31 + T_19 * t_34;
	FINTERVAL T_49 = T_45 * t_39 + T_46 * t_41;
	FINTERVAL T_50 = T_45 * t_40 + T_46 * t_42;
	FINTERVAL T_51 = fuda.C_b72 * T_45 + fuda.C_b83 * T_46;
	FINTERVAL T_52 = T_47 * t_39 + T_48 * t_41;
	FINTERVAL T_53 = T_47 * t_40 + T_48 * t_42;
	FINTERVAL T_54 = fuda.C_b72 * T_47 + fuda.C_b83 * T_48;
	double t_55 = fuda.C_c11 * sin(fuda.thet_h - fuda.thet_1) + fuda.C_c12;
	double t_56 = fuda.C_c21 * sin(fuda.thet_h + fuda.thet_c2) + fuda.C_c22;
	double t_57 = fuda.C_c31 * sin(fuda.thet_h - fuda.thet_1) + fuda.C_c32;
	FINTERVAL T_58 = t_55 * T_3 + t_56 * T_9 + t_57;
	double t_59 = fuda.C_c41 * cos(fuda.thet_h - fuda.thet_1);
	double t_60 = fuda.C_c42 * sin(fuda.thet_h + fuda.thet_c4);
	double t_61 = fuda.C_c43 * cos(fuda.thet_h - fuda.thet_1);
	FINTERVAL T_62 = t_59 * T_3 + t_60 * T_9 + t_61;
	double t_63 = fuda.C_c51 * cos(fuda.thet_h - fuda.thet_1);
	double t_64 = fuda.C_c52 * sin(fuda.thet_h + fuda.thet_c5);
	double t_65 = fuda.C_c53 * cos(fuda.thet_h - fuda.thet_1) + fuda.C_c54;
	FINTERVAL T_66 = t_63 * T_3 + t_64 * T_9 + t_65;
	double t_70 = t_55 * t_7 + t_56 * t_11;
	double t_74 = t_59 * t_7 + t_60 * t_11;
	double t_78 = t_63 * t_7 + t_64 * t_11;
	FINTERVAL T_80 = fuda.C_s11 * T_3 + fuda.C_s12 * T_9 + fuda.C_s13;
	double t_82 = fuda.C_s11 + fuda.C_s12 * t_11;
	FINTERVAL T_83 = T_49 * T_58 + T_50 * T_62 + T_51 * T_66 + T_80;
	FINTERVAL T_84 = T_52 * T_58 + T_53 * T_62 + T_54 * T_66 
		+ T_49 * t_70 + T_50 * t_74 + T_51 * t_78 + t_82;
	return T_84;
}

long DEHW::CILFOFE(double thet_1, double x_d, double y_d, 
	double &Psi_1, double &kapp_1xd, double &kapp_1yd, double &tau_1xd){
	//
	double thet_c = thet_1 / i_1c;
	//
	double kapp_cxd = 0.0;
	double kapp_cyd = 0.0;
	double tau_cxd = 0.0;
	Vector3d i_d_c, j_d_c, i_d_oc, j_d_oc, omeg_c1_oc, v_c1_oc;
	i_d_c << -1.0, 0.0, 0.0;
	j_d_c << 0.0, - sin(beta_c), cos(beta_c);
	//R_oc,c
	Matrix<double,3,3> rota;
	rota << cos(thet_c), -sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	i_d_oc = rota * i_d_c;
	j_d_oc = rota * j_d_c;
	omeg_c1_oc << 0.0, -1.0, i_c1;
	v_c1_oc << - y_d * cos(beta_c) 
		- i_c1 * (- x_d * sin(thet_c) + cos(thet_c) * (r_b2 - y_d * sin(beta_c))),
		i_c1 * (- x_d * cos(thet_c) - sin(thet_c) * (r_b2 - y_d * sin(beta_c))),
		- x_d * cos(thet_c) - sin(thet_c) * (r_b2 - y_d * sin(beta_c)) + a_1c;
	double N_1xd = kapp_cxd * v_c1_oc.dot(i_d_oc) + tau_cxd * v_c1_oc.dot(j_d_oc)
		+ omeg_c1_oc.dot(j_d_oc);
	double N_1yd = tau_cxd * v_c1_oc.dot(i_d_oc) + kapp_cyd * v_c1_oc.dot(j_d_oc)
		- omeg_c1_oc.dot(i_d_oc);
	Vector3d N_1_oc = N_1xd * i_d_oc + N_1yd * j_d_oc;
	double PPhi_1Pthet_1 = x_d * sin(beta_c) * sin(thet_c) / i_1c
		+ y_d * cos(thet_c) / i_1c
		- r_b2 * sin(beta_c) * cos(thet_c) / i_1c;
	Psi_1 = N_1_oc.dot(v_c1_oc) + PPhi_1Pthet_1;
	double kapp_c1xd = N_1xd * N_1xd / Psi_1;
	double kapp_c1yd = N_1yd * N_1yd / Psi_1;	
	double tau_c1xd = N_1xd * N_1yd / Psi_1;
	kapp_1xd = kapp_cxd - kapp_c1xd;
	kapp_1yd = kapp_cyd - kapp_c1yd;
	tau_1xd = tau_cxd - tau_c1xd;
	return 1;
}

long DEHW::PD_CILFOFE(double thet_1, double x_d, double y_d, double Px_d, double Py_d, 
	double &Psi_1, double &kapp_1xd, double &kapp_1yd, double &tau_1xd,
	double &PPsi_1, double &Pkapp_1xd, double &Pkapp_1yd, double &Ptau_1xd){
	//
	double thet_c = thet_1 / i_1c;
	//
	double kapp_cxd = 0.0;
	double kapp_cyd = 0.0;
	double tau_cxd = 0.0;
	Vector3d i_d_c, j_d_c, i_d_oc, j_d_oc, omeg_c1_oc, v_c1_oc;
	Vector3d Pi_d_oc, Pj_d_oc, Pv_c1_oc;
	i_d_c << -1.0, 0.0, 0.0;
	j_d_c << 0.0, - sin(beta_c), cos(beta_c);
	//R_oc,c
	Matrix<double,3,3> rota, Prota;
	rota << cos(thet_c), -sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	Prota << - sin(thet_c) * i_c1, - cos(thet_c) * i_c1, 0.0,
		cos(thet_c) * i_c1, - sin(thet_c) * i_c1, 0.0,
		0.0, 0.0, 0.0;
	Pi_d_oc = Prota * i_d_c;
	Pj_d_oc = Prota * j_d_c;
	i_d_oc = rota * i_d_c;
	j_d_oc = rota * j_d_c;
	omeg_c1_oc << 0.0, -1.0, i_c1;
	v_c1_oc << - y_d * cos(beta_c) 
		- i_c1 * (- x_d * sin(thet_c) + cos(thet_c) * (r_b2 - y_d * sin(beta_c))),
		i_c1 * (- x_d * cos(thet_c) - sin(thet_c) * (r_b2 - y_d * sin(beta_c))),
		- x_d * cos(thet_c) - sin(thet_c) * (r_b2 - y_d * sin(beta_c)) + a_1c;
	Pv_c1_oc << - Py_d * cos(beta_c) 
		- i_c1 * (- Px_d * sin(thet_c) - x_d * cos(thet_c) * i_c1 
		- sin(thet_c) * i_c1 * (r_b2 - y_d * sin(beta_c))
		+ cos(thet_c) * (- Py_d * sin(beta_c))),
		i_c1 * (- Px_d * cos(thet_c) + x_d * sin(thet_c) * i_c1
		- cos(thet_c) * i_c1 * (r_b2 - y_d * sin(beta_c))
		- sin(thet_c) * (- Py_d * sin(beta_c))),
		- Px_d * cos(thet_c) + x_d * sin(thet_c) * i_c1 
		- cos(thet_c) * i_c1 * (r_b2 - y_d * sin(beta_c))
		- sin(thet_c) * (- Py_d * sin(beta_c));
	double N_1xd = kapp_cxd * v_c1_oc.dot(i_d_oc) + tau_cxd * v_c1_oc.dot(j_d_oc)
		+ omeg_c1_oc.dot(j_d_oc);
	double N_1yd = tau_cxd * v_c1_oc.dot(i_d_oc) + kapp_cyd * v_c1_oc.dot(j_d_oc)
		- omeg_c1_oc.dot(i_d_oc);
	double PN_1xd = omeg_c1_oc.dot(Pj_d_oc);
	double PN_1yd = - omeg_c1_oc.dot(Pi_d_oc);
	Vector3d N_1_oc = N_1xd * i_d_oc + N_1yd * j_d_oc;
	Vector3d PN_1_oc = PN_1xd * i_d_oc + PN_1yd * j_d_oc 
		+ N_1xd * Pi_d_oc + N_1yd * Pj_d_oc;
	double PPhi_1Pthet_1 = x_d * sin(beta_c) * sin(thet_c) * i_c1
		+ y_d * cos(thet_c) * i_c1
		- r_b2 * sin(beta_c) * cos(thet_c) * i_c1;
	double PPPhi_1Pthet_1 = Px_d * sin(beta_c) * sin(thet_c) * i_c1
		+ x_d * sin(beta_c) * cos(thet_c) * i_c1 * i_c1
		+ Py_d * cos(thet_c) * i_c1
		- y_d * sin(thet_c) * i_c1 * i_c1
		+ r_b2 * sin(beta_c) * sin(thet_c) * i_c1 * i_c1;
	Psi_1 = N_1_oc.dot(v_c1_oc) + PPhi_1Pthet_1;
	PPsi_1 = PN_1_oc.dot(v_c1_oc) + N_1_oc.dot(Pv_c1_oc) + PPPhi_1Pthet_1;
	double kapp_c1xd = N_1xd * N_1xd / Psi_1;
	double kapp_c1yd = N_1yd * N_1yd / Psi_1;	
	double tau_c1xd = N_1xd * N_1yd / Psi_1;
	double Pkapp_c1xd = (2.0 * N_1xd * PN_1xd * Psi_1 - N_1xd * N_1xd * PPsi_1) / (Psi_1 * Psi_1);
	double Pkapp_c1yd = (2.0 * N_1yd * PN_1yd * Psi_1 - N_1yd * N_1yd * PPsi_1) / (Psi_1 * Psi_1);	
	double Ptau_c1xd = (PN_1xd * N_1yd * Psi_1 + N_1xd * PN_1yd * Psi_1 
		- N_1xd * N_1yd * PPsi_1) / (Psi_1 * Psi_1);
	kapp_1xd = kapp_cxd - kapp_c1xd;
	kapp_1yd = kapp_cyd - kapp_c1yd;
	tau_1xd = tau_cxd - tau_c1xd;
	Pkapp_1xd = - Pkapp_c1xd;
	Pkapp_1yd = - Pkapp_c1yd;
	Ptau_1xd = - Ptau_c1xd;
	return 1;
}

double DEHW::CILFOSE_NI(double thet_1, double thet_h, double &kapp_h2N){
	//
	double thet_c = thet_1 / i_1c;
	double thet_2 = thet_h / i_h2;
	//
	double x_d, y_d;
	FSME(thet_1, thet_h, x_d, y_d);
	double Psi_1, kapp_1xd, kapp_1yd, tau_1xd;
	CILFOFE(thet_1, x_d, y_d, Psi_1, kapp_1xd, kapp_1yd, tau_1xd);
	double kapp_hxd = kapp_1xd;
	double kapp_hyd = kapp_1yd;
	double tau_hxd = tau_1xd;
	//
	Vector3d i_d_c, j_d_c, i_d_oh, j_d_oh;
	i_d_c << -1.0, 0.0, 0.0;
	j_d_c << 0.0, - sin(beta_c), cos(beta_c);
	Vector3d r_c_c, r_h_oh;
	r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
	//R_oc,c
	Matrix<double,3,3> rota;
	rota << cos(thet_c), -sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	i_d_oh = rota * i_d_c;
	j_d_oh = rota * j_d_c;
	r_h_oh = rota * r_c_c;
	//R_o1,oc
	rota << 1.0, 0.0, 0.0,
		0.0, 0.0, -1.0,
		0.0, 1.0, 0.0;
	i_d_oh = rota * i_d_oh;
	j_d_oh = rota * j_d_oh;
	r_h_oh = rota * r_h_oh;
	//T_o1,oc
	r_h_oh(0) = r_h_oh(0) + a_1c;
	//R_1,o1
	rota << cos(thet_1), sin(thet_1), 0.0,
		-sin(thet_1), cos(thet_1), 0.0,
		0.0, 0.0, 1.0;
	i_d_oh = rota * i_d_oh;
	j_d_oh = rota * j_d_oh;
	r_h_oh = rota * r_h_oh;
	//R_oh,h
	rota << cos(thet_h), -sin(thet_h), 0.0,
		sin(thet_h), cos(thet_h), 0.0,
		0.0, 0.0, 1.0;
	i_d_oh = rota * i_d_oh;
	j_d_oh = rota * j_d_oh;
	r_h_oh = rota * r_h_oh;
	//
	Vector3d omeg_h2_oh, omeg_2_oh, o_h2_oh, v_h2_oh;
	omeg_h2_oh << 0.0, i_2h, 1.0;
	omeg_2_oh << 0.0, - i_2h, 0.0;
	o_h2_oh << - a_h2, 0.0, 0.0;
	v_h2_oh = omeg_h2_oh.cross(r_h_oh) - omeg_2_oh.cross(o_h2_oh);
	double N_2xd = kapp_hxd * v_h2_oh.dot(i_d_oh) + tau_hxd * v_h2_oh.dot(j_d_oh)
		+omeg_h2_oh.dot(j_d_oh);
	double N_2yd = tau_hxd * v_h2_oh.dot(i_d_oh) + kapp_hyd * v_h2_oh.dot(j_d_oh)
		-omeg_h2_oh.dot(i_d_oh);
	Vector3d N_2_oh = N_2xd * i_d_oh + N_2yd * j_d_oh;
	double B_11 = i_2h * x_d * cos(beta_c) - i_2h * a_1c * cos(beta_c) * cos(thet_c);
	double B_12 = - i_2h * x_d * sin(beta_c) * sin(thet_c) 
		- i_2h * y_d * cos(thet_c) + i_2h * r_b2 * sin(beta_c) * cos(thet_c);
	double PPhi_2Pthet_h = - B_11 * sin(thet_h - thet_1) + B_12 * cos(thet_h - thet_1);
	double Psi_2 = N_2_oh.dot(v_h2_oh) + PPhi_2Pthet_h;
	kapp_h2N = (N_2xd * N_2xd + N_2yd * N_2yd) / Psi_2;
	return Psi_2;
}

double DEHW::PD_CILFOSE_NI(double thet_1, double thet_h){
	//
	double thet_c = thet_1 / i_1c;
	//
	double x_d, y_d, Px_d, Py_d;
	Matrix<double,2,2> Pxy_d;
	PD_FSME(thet_1, thet_h, x_d, y_d, Pxy_d);
	Px_d = Pxy_d(0, 0);
	Py_d = Pxy_d(1, 0);
	double Psi_1, kapp_1xd, kapp_1yd, tau_1xd, PPsi_1, Pkapp_1xd, Pkapp_1yd, Ptau_1xd;
	PD_CILFOFE(thet_1, x_d, y_d, Px_d, Py_d, 
		Psi_1, kapp_1xd, kapp_1yd, tau_1xd, 
		PPsi_1, Pkapp_1xd, Pkapp_1yd, Ptau_1xd);
	double kapp_hxd = kapp_1xd;
	double kapp_hyd = kapp_1yd;
	double tau_hxd = tau_1xd;
	double Pkapp_hxd = Pkapp_1xd;
	double Pkapp_hyd = Pkapp_1yd;
	double Ptau_hxd = Ptau_1xd;
	//
	Vector3d i_d_c, j_d_c, i_d_oh, j_d_oh;
	Vector3d Pi_d_oh, Pj_d_oh;
	i_d_c << -1.0, 0.0, 0.0;
	j_d_c << 0.0, - sin(beta_c), cos(beta_c);
	Vector3d r_c_c, r_h_oh, Pr_c_c, Pr_h_oh;
	r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
	Pr_c_c << - Px_d, - Py_d * sin(beta_c), Py_d * cos(beta_c);
	//R_oc,c
	Matrix<double,3,3> rota, Prota;
	rota << cos(thet_c), -sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	Prota << -sin(thet_c) * i_c1, -cos(thet_c) * i_c1, 0.0,
		cos(thet_c) * i_c1, -sin(thet_c) * i_c1, 0.0,
		0.0, 0.0, 0.0;
	Pi_d_oh = Prota * i_d_c;
	Pj_d_oh = Prota * j_d_c;
	Pr_h_oh = Prota * r_c_c + rota * Pr_c_c;
	i_d_oh = rota * i_d_c;
	j_d_oh = rota * j_d_c;
	r_h_oh = rota * r_c_c;
	//R_o1,oc
	rota << 1.0, 0.0, 0.0,
		0.0, 0.0, -1.0,
		0.0, 1.0, 0.0;
	Pi_d_oh = rota * Pi_d_oh;
	Pj_d_oh = rota * Pj_d_oh;
	Pr_h_oh = rota * Pr_h_oh;
	i_d_oh = rota * i_d_oh;
	j_d_oh = rota * j_d_oh;
	r_h_oh = rota * r_h_oh;
	//T_o1,oc
	r_h_oh(0) = r_h_oh(0) + a_1c;
	//R_1,o1
	rota << cos(thet_1), sin(thet_1), 0.0,
		-sin(thet_1), cos(thet_1), 0.0,
		0.0, 0.0, 1.0;
	Prota << -sin(thet_1), cos(thet_1), 0.0,
		-cos(thet_1), -sin(thet_1), 0.0,
		0.0, 0.0, 0.0;
	Pi_d_oh = Prota * i_d_oh + rota * Pi_d_oh;
	Pj_d_oh = Prota * j_d_oh + rota * Pj_d_oh;
	Pr_h_oh = Prota * r_h_oh + rota * Pr_h_oh;
	i_d_oh = rota * i_d_oh;
	j_d_oh = rota * j_d_oh;
	r_h_oh = rota * r_h_oh;
	//R_oh,h
	rota << cos(thet_h), -sin(thet_h), 0.0,
		sin(thet_h), cos(thet_h), 0.0,
		0.0, 0.0, 1.0;
	Pi_d_oh = rota * Pi_d_oh;
	Pj_d_oh = rota * Pj_d_oh;
	Pr_h_oh = rota * Pr_h_oh;
	i_d_oh = rota * i_d_oh;
	j_d_oh = rota * j_d_oh;
	r_h_oh = rota * r_h_oh;
	//
	Vector3d omeg_h2_oh, omeg_2_oh, o_h2_oh, v_h2_oh, Pv_h2_oh;
	omeg_h2_oh << 0.0, i_2h, 1.0;
	omeg_2_oh << 0.0, - i_2h, 0.0;
	o_h2_oh << - a_h2, 0.0, 0.0;
	v_h2_oh = omeg_h2_oh.cross(r_h_oh) - omeg_2_oh.cross(o_h2_oh);
	Pv_h2_oh = omeg_h2_oh.cross(Pr_h_oh);
	double N_2xd = kapp_hxd * v_h2_oh.dot(i_d_oh) + tau_hxd * v_h2_oh.dot(j_d_oh)
		+omeg_h2_oh.dot(j_d_oh);
	double N_2yd = tau_hxd * v_h2_oh.dot(i_d_oh) + kapp_hyd * v_h2_oh.dot(j_d_oh)
		-omeg_h2_oh.dot(i_d_oh);
	double PN_2xd = Pkapp_hxd * v_h2_oh.dot(i_d_oh) 
		+ kapp_hxd * Pv_h2_oh.dot(i_d_oh)
		+ kapp_hxd * v_h2_oh.dot(Pi_d_oh)
		+ Ptau_hxd * v_h2_oh.dot(j_d_oh)
		+ tau_hxd * Pv_h2_oh.dot(j_d_oh)
		+ tau_hxd * v_h2_oh.dot(Pj_d_oh)
		+ omeg_h2_oh.dot(Pj_d_oh);
	double PN_2yd = Ptau_hxd * v_h2_oh.dot(i_d_oh) 
		+ tau_hxd * Pv_h2_oh.dot(i_d_oh) 
		+ tau_hxd * v_h2_oh.dot(Pi_d_oh)
		+ Pkapp_hyd * v_h2_oh.dot(j_d_oh)
		+ kapp_hyd * Pv_h2_oh.dot(j_d_oh)
		+ kapp_hyd * v_h2_oh.dot(Pj_d_oh)
		- omeg_h2_oh.dot(Pi_d_oh);
	Vector3d N_2_oh = N_2xd * i_d_oh + N_2yd * j_d_oh;
	Vector3d PN_2_oh = PN_2xd * i_d_oh + PN_2yd * j_d_oh + N_2xd * Pi_d_oh + N_2yd * Pj_d_oh;
	double B_11 = i_2h * x_d * cos(beta_c) - i_2h * a_1c * cos(beta_c) * cos(thet_c);
	double B_12 = - i_2h * x_d * sin(beta_c) * sin(thet_c) 
		- i_2h * y_d * cos(thet_c) + i_2h * r_b2 * sin(beta_c) * cos(thet_c);
	double PB_11 = i_2h * Px_d * cos(beta_c) + i_2h * a_1c * cos(beta_c) * sin(thet_c) * i_c1;
	double PB_12 = - i_2h * Px_d * sin(beta_c) * sin(thet_c) 
		- i_2h * x_d * sin(beta_c) * cos(thet_c) * i_c1 
		- i_2h * Py_d * cos(thet_c) + i_2h * y_d * sin(thet_c) * i_c1 
		- i_2h * r_b2 * sin(beta_c) * sin(thet_c) * i_c1;
	double PPhi_2Pthet_h = - B_11 * sin(thet_h - thet_1) + B_12 * cos(thet_h - thet_1);
	double PPPhi_2Pthet_h = - PB_11 * sin(thet_h - thet_1) + PB_12 * cos(thet_h - thet_1)
		 + B_11 * cos(thet_h - thet_1) + B_12 * sin(thet_h - thet_1);
	double Psi_2 = N_2_oh.dot(v_h2_oh) + PPhi_2Pthet_h;
	double PPsi_2 = PN_2_oh.dot(v_h2_oh) + N_2_oh.dot(Pv_h2_oh) + PPPhi_2Pthet_h;
	return PPsi_2;
}

long DEHW::ON_WHEEL_TS(double x_d, double y_d, double thet_1, Vector3d r_2_2){
	//xi_21, xi_22
	double angl_fi, radi_fi, R_fmini, R_fmaxi;
	WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
	if(!(-widtAngl <= angl_fi && angl_fi <= widtAngl 
		&& R_fmini <= radi_fi && radi_fi <= R_fmaxi)){
		return -1;
	}
	//first and last cutting edge of hob
	double thet_c = i_c1 * thet_1;
	Vector3d r_c_c, T_o1_oc, r_1_o1;
	Matrix<double,3,3> R_oc_c, R_o1_oc, R_1_o1;
	r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
	R_oc_c << cos(thet_c), - sin(thet_c), 0.0,
		sin(thet_c), cos(thet_c), 0.0,
		0.0, 0.0, 1.0;
	R_o1_oc << 1.0, 0.0, 0.0,
		0.0, 0.0, -1.0,
		0.0, 1.0, 0.0;
	R_1_o1 << cos(thet_1), sin(thet_1), 0.0,
		- sin(thet_1), cos(thet_1), 0.0,
		0.0, 0.0, 1.0;
	T_o1_oc << a_1c, 0.0, 0.0;
	r_1_o1 = R_o1_oc * R_oc_c * r_c_c + T_o1_oc;
	double woxi_11 =  thet_1 - atan2(r_1_o1(1), r_1_o1(0));
	if(!(wormCurv(0) - 1.0E-14 <= woxi_11 && woxi_11 <= wormCurv(2) + 1.0E-14)){
		return -1;
	}
	//tooth tip of worm
	Vector3d r_1_1;
	WORM_DC2R(x_d, y_d, thet_c, r_1_1);
	double xi_12 = sqrt(r_1_1(2) * r_1_1(2) + pow(
		a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1)), 2.0));
	if(xi_12 < R_a(0)){
		return -1;
	}
	return 1;
}

double DEHW::TICL_LENGTH(double thet_h, const Matrix<double,Dynamic,Dynamic> &sd12){
	double totaLeng = 0.0;
	//no check of input parameter validity
	//
	Vector3d r_2_2Back;
	long flagBack;
	if(modiTran == 0.0 && modiCent == 0.0){
		long N = 99999;
		double x_dH = a_1c;
		double x_dL = - a_1c;
		double thet_c = i_c1 * thet_h;
		double deltX_d = (x_dH - x_dL) / (N + 1.0);
		for(long tj = 1; tj <= N; tj ++){
			double x_d = x_dL + deltX_d * (double)tj;
			double C_11 = - sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c);
			double C_12 = sin(thet_c);
			double C_13 = - r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c);
			double y_d = (- C_13 - C_11 * x_d) / C_12;
			Vector3d r_2_2;
			WHEE_1H2R(x_d, y_d, i_1c * thet_c, thet_h, r_2_2);
			long flagTemp = ON_WHEEL_TS(x_d, y_d, thet_h, r_2_2);
			if(tj >= 2 && flagBack == 1 && flagTemp == 1){
				totaLeng += (r_2_2 - r_2_2Back).norm();
			}
			flagBack = flagTemp;
			r_2_2Back = r_2_2;
		}
	}
	//
	double thet_cL, thet_cH;
	SINGULAR_H2C(thet_h, thet_cL, thet_cH);
	//
	long N = 99999;
	double deltThet_c = (thet_cH - thet_cL) / (N + 1.0);
	default_random_engine dren;
	uniform_real_distribution<double> urdd(0.25, 0.75);
	for(long tj = 1; tj <= N; tj ++){
		//
		double thet_c = thet_cL + deltThet_c * (double)tj;
		//
		double thet_hm;
		while(true){
			double thet_hstemp;
			SINGULAR_C2H(thet_c, thet_hstemp, thet_hm);
			if(abs(thet_h - thet_hm) <= 1.0E-14){
				thet_c = thet_cL + deltThet_c * (tj - urdd(dren));
			}
			else{
				break;
			}
		}
		//
		double thet_1 = i_1c * thet_c;
		double x_d, y_d;
		FSME(thet_1, thet_h, x_d, y_d);
		Vector3d r_2_2;
		WHEE_1H2R(x_d, y_d, thet_1, thet_h, r_2_2);
		long flagTemp = ON_WHEEL_TS(x_d, y_d, thet_1, r_2_2);
		if(modiTran != 0.0 || modiCent != 0.0){
			double angl_fi, radi_fi, R_fmini, R_fmaxi;
			WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
			for(long tk = 0; tk <= sd12.rows() - 2; tk ++){
				if(sd12(tk,0) <= radi_fi && radi_fi <= sd12(tk + 1,0)){
					double inteAngl = sd12(tk,1) + (sd12(tk + 1,1) - sd12(tk,1)) 
						/ (sd12(tk + 1,0) - sd12(tk,0)) * (radi_fi - sd12(tk,0));
					if((thet_h >= thet_hm && angl_fi > inteAngl) ||
						(thet_h < thet_hm && angl_fi < inteAngl)){
						flagTemp = -1;
					}
					break;
				}
			}
			if(thet_h >= thet_hm && ((radi_fi < sd12(0,0) && angl_fi > sd12(0,1)) ||
				(radi_fi > sd12(sd12.rows() - 1,0) && angl_fi > sd12(sd12.rows() - 1,1)))){
				flagTemp = -1;
			}
		}
		double kapp_h2N;
		if(CILFOSE_NI(thet_1, thet_h, kapp_h2N) <= 0.0){
			flagTemp = -1;
		}
		if(tj >= 2 && flagBack == 1 && flagTemp == 1){
			totaLeng += (r_2_2 - r_2_2Back).norm();
		}
		flagBack = flagTemp;
		r_2_2Back = r_2_2;
	}
	return totaLeng;
}

long DEHW::AC_PRESSURE(double thet_1, double thet_h, double &sigm_H, double &b_H){
	//induced normal curvature
	double kapp_h2N;
	CILFOSE_NI(thet_1, thet_h, kapp_h2N);
	//D1_1 intersects with D1_2
	Matrix<double,Dynamic,Dynamic> sd12 = MatrixXd::Zero(curvCoor(1).rows(), 2);
	for(long ti = 0; ti < curvCoor(1).cols(); ti ++){
		sd12(ti,0) = curvCoor(1)(curvCoor(1).rows() - 1, ti)(1);
		sd12(ti,1) = curvCoor(1)(curvCoor(1).rows() - 1, ti)(0);
		for(long tj = curvCoor(1).rows() - 1; tj > 0; tj --){
			if(fpha(tj,ti) == 2 && fpha(tj-1,ti) != 2){
				sd12(ti,0) = (curvCoor(1)(tj,ti)(1) + curvCoor(1)(tj-1,ti)(1)) / 2.0;
				sd12(ti,1) = (curvCoor(1)(tj,ti)(0) + curvCoor(1)(tj-1,ti)(0)) / 2.0;
				break;
			}
		}
	}
	double totaLeng = 0.0;
	for(long ti = -2; ti <= 5; ti ++){
		double thet_h_ti = ti * 2.0 * PI / z(0) + thet_h;
		if(thet_h_ti < PI){
			continue;
		}
		totaLeng += TICL_LENGTH(thet_h_ti, sd12);
	}
	double alph_n = atan(tan(alph) * cos(leadAngl));
	double outpTorq = inpuTorq * z(1) / (double)(z(0));
	double W_n = outpTorq / (d(1) / 2.0) / cos(alph_n) / cos(leadAngl) / totaLeng;
	//composite elastic modulus
	double E_s = 1.0 / ((1.0 - vem(0).matePois * vem(0).matePois) / vem(0).mateElas 
		+ (1.0 - vem(1).matePois * vem(1).matePois) / vem(1).mateElas);
	//maximum contact pressure
	sigm_H = sqrt(E_s * kapp_h2N * W_n / PI);
	//half contact band width
	b_H = sqrt(4.0 * W_n / PI / kapp_h2N / E_s);
	cout<<setprecision(20);
	cout<<"totaLeng="<<totaLeng<<",outpTorq="<<outpTorq<<endl
		<<"W_n="<<W_n<<",E_s="<<E_s<<endl
		<<"kapp_h2N="<<kapp_h2N<<",sigm_H="<<sigm_H<<",b_H="<<b_H<<endl;
	return 1;
}

long DEHW::ANALYTICAL_CONTACT(){
	OUTPUT_TIME("ANALYTICAL_CONTACT");
	//D1_1 intersects with D1_2
	Matrix<double,Dynamic,Dynamic> sd12 = MatrixXd::Zero(curvCoor(1).cols(), 2);
	for(long ti = 0; ti < curvCoor(1).cols(); ti ++){
		sd12(ti,0) = curvCoor(1)(curvCoor(1).rows() - 1, ti)(1);
		sd12(ti,1) = curvCoor(1)(curvCoor(1).rows() - 1, ti)(0);
		for(long tj = curvCoor(1).rows() - 1; tj > 0; tj --){
			if(fpha(tj,ti) == 2 && fpha(tj-1,ti) != 2){
				sd12(ti,0) = (curvCoor(1)(tj,ti)(1) + curvCoor(1)(tj-1,ti)(1)) / 2.0;
				sd12(ti,1) = (curvCoor(1)(tj,ti)(0) + curvCoor(1)(tj-1,ti)(0)) / 2.0;
				break;
			}
		}
	}
	//
	ofstream tempOfst;
	tempOfst.open(DIRECTORY("resuTILCLENG.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	analNumb = 40;
	for(long ta = 0; ta < analNumb; ta ++){
		double totaLeng = 0.0;
		for(long ti = -2; ti <= 5; ti ++){
			cout<<"ta,ti="<<ta<<","<<ti<<endl;
			double thet_h = ti * 2.0 * PI / z(0) + 2.0 * PI / z(0) / analNumb * ta;
			if(thet_h < PI){
				continue;
			}
			totaLeng += TICL_LENGTH(thet_h, sd12);
		}
		tempOfst<<setw(30)<<totaLeng<<endl;
	}
	tempOfst.close();
	return 1;
}

long DEHW::TS_GRID(){
	//
	WORM_TS_GRID();
	ofstream tempOfst;
	tempOfst.open(DIRECTORY("resuWOTSGRID.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	for(long ti = 0; ti < curvCoor(0).rows(); ti ++){
		for(long tj = 0; tj < curvCoor(0).cols(); tj ++){
			tempOfst << setw(30) << cartCoor(0)(ti,tj)(0) 
				<< setw(30) << cartCoor(0)(ti,tj)(1) 
				<< setw(30) << cartCoor(0)(ti,tj)(2) << endl;
		}
	}
	tempOfst.close();
	//
	WHEE_TS_GRID();
	tempOfst.open(DIRECTORY("resuWHTSGRID.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	for(long ti = 0; ti < curvCoor(1).rows(); ti ++){
		for(long tj = 0; tj < curvCoor(1).cols(); tj ++){
			tempOfst << setw(30) << cartCoor(1)(ti,tj)(0) 
				<< setw(30) << cartCoor(1)(ti,tj)(1) 
				<< setw(30) << cartCoor(1)(ti,tj)(2) << endl;
		}
	}
	tempOfst.close();
	tempOfst.open(DIRECTORY("resuWHTSPHAS.txt"), ios::out);
	tempOfst<<setiosflags(ios::scientific)<<setprecision(20);
	for(long ti = 0; ti < curvCoor(1).rows(); ti ++){
		for(long tj = 0; tj < curvCoor(1).cols(); tj ++){
			tempOfst << setw(12) << fpha(ti,tj) << endl;
		}
	}
	tempOfst.close();
	return 1;
}

long DEHW::WORM_TS_GRID(){
	OUTPUT_TIME("WORM_TS_GRID");
	//curvilinear coordinate
	long tangNumb = (long)(130.0 * (1.0 + (z(0) - 1.0) / 4.0));
	tangNumb = tangNumb - tangNumb % z(0);
	double deltTang = 2.0 * PI / tangNumb;
	long moduNumb = ceil((wormCurv(1) - wormCurv(0)) / deltTang);
	double realStar = wormCurv(1) - moduNumb * deltTang;
	gridNumb(0,0) = 10;
	gridNumb(0,1) = 5;
	gridNumb(0,2) = 4;
	gridNumb(0,3) = 8;
	gridNumb(0,4) = tangNumb;
	gridNumb(0,5) = tangNumb * floor((wormCurv(2) - wormCurv(1)) / 2.0 / PI) + 2 * moduNumb;
	curvCoor(0).resize(gridNumb(0,5) * pow2[refiLeve] + 1, gridNumb(0,3) * pow2[refiLeve] + 1);
	deltTang /= pow2[refiLeve];
	for(long ti = 0; ti < curvCoor(0).rows(); ti ++){
		for(long tj = 0; tj < curvCoor(0).cols(); tj ++){
			curvCoor(0)(ti,tj) << realStar + ti * deltTang,
				R_t(0) + (R_a(0) - R_t(0)) / (curvCoor(0).cols() - 1) * (double)tj;
		}
	}
	//Cartesian coordinate
	cartCoor(0).resize(curvCoor(0).rows(), curvCoor(0).cols());
	for(long ti = 0; ti < curvCoor(0).rows(); ti ++){
		if(ti % 1000 == 0){
			cout<<ti<<"/"<<curvCoor(0).rows()<<endl;
		}
		for(long tj = 0; tj < curvCoor(0).cols(); tj ++){
			double thet_c;
			WORM_CURV_2_CART(curvCoor(0)(ti,tj)(0), 
				curvCoor(0)(ti,tj)(1), cartCoor(0)(ti,tj), thet_c);
			if(reliSwit == 1){
				WORM_RELI(cartCoor(0)(ti,tj), ti, tj);
			}
		}
	}
	return 1;
}

long DEHW::WHEE_TS_GRID(){
	OUTPUT_TIME("WHEE_TS_GRID");
	//curvilinear coordinate
	gridNumb(1,0) = 10;
	gridNumb(1,1) = 8;
	gridNumb(1,2) = 4;
	gridNumb(1,3) = 8;
	gridNumb(1,4) = 32;
	gridNumb(1,5) = 8 + z(0);
	curvCoor(1).resize(gridNumb(1,4) * pow2[refiLeve] + 1, gridNumb(1,3) * pow2[refiLeve] + 1);
	for(long ti = 0; ti < curvCoor(1).rows(); ti ++){
		double angl_fi = widtAngl - 2.0 * widtAngl / (curvCoor(1).rows() - 1) * (double)ti;
		double angl_ai = angl_fi - asin(offsR_a * sin(angl_fi) / R_a(1));
		double R_fmini = (R_a(1) * cos(angl_ai) - offsR_a) / cos(angl_fi);
		double R_fmaxi = R_t(1);
		for(long tj = 0; tj < curvCoor(1).cols(); tj ++){
			curvCoor(1)(ti,tj) << angl_fi, 
				R_fmini + (R_fmaxi - R_fmini) / (curvCoor(1).cols() - 1) * (double)tj;
		}
	}
	//Cartesian coordinate	
	cartCoor(1).resize(curvCoor(1).rows(), curvCoor(1).cols());
	fpha = MatrixXi::Zero(curvCoor(1).rows(), curvCoor(1).cols());
	NEW_CONT_ZONE(1);//left new contact zone
	NEW_CONT_ZONE(2);//right new contact zone
	if(modiTran == 0.0 && modiCent == 0.0){
		FORMER_CONT_ZONE();//former contact zone
	}
	TRANSITION_ZONE(1);//head transition zone
	TRANSITION_ZONE(2);//rear transition zone
	//tooth flank relief
	if(reliSwit == 1){
		for(long ti = 0; ti < curvCoor(1).rows(); ti ++){
			for(long tj = 0; tj < curvCoor(1).cols(); tj ++){
				WHEE_RELI(cartCoor(1)(ti,tj), ti, tj);
			}
		}
	}
	return 1;
}

long DEHW::WORM_CURV_2_CART(double xi_11, double xi_12, Vector3d &r_1_1, double &thet_c){
	Matrix<double,2,1> x;//0 - thet_c, 1 - x_d
	x << i_c1 * xi_11, d(1) / 2.0;//
	while(true){
		Matrix<double,2,1> func, deltX;
		Matrix<double,2,2> Dfunc;
		//function
		thet_c = x(0);
		double thet_1 = i_1c * thet_c;
		double x_d = x(1);
		double y_d = - ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * x_d 
			+ (- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c))) / sin(thet_c);
		Vector3d r_c_c, T_o1_oc, r_1_o1;
		Matrix<double,3,3> R_oc_c, R_o1_oc, R_1_o1;
		r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
		R_oc_c << cos(thet_c), - sin(thet_c), 0.0,
			sin(thet_c), cos(thet_c), 0.0,
			0.0, 0.0, 1.0;
		R_o1_oc << 1.0, 0.0, 0.0,
			0.0, 0.0, -1.0,
			0.0, 1.0, 0.0;
		R_1_o1 << cos(thet_1), sin(thet_1), 0.0,
			- sin(thet_1), cos(thet_1), 0.0,
			0.0, 0.0, 1.0;
		T_o1_oc << a_1c, 0.0, 0.0;
		r_1_o1 = R_o1_oc * R_oc_c * r_c_c + T_o1_oc;
		r_1_1 = R_1_o1 * r_1_o1;
		func << thet_1 - atan2(r_1_o1(1), r_1_o1(0)) - xi_11,
			r_1_1(2) * r_1_1(2) 
			+ (a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1))) 
			* (a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1)))
			- xi_12 * xi_12;
		//derivative 1
		double Dy_d = - (((+ sin(beta_c) * sin(thet_c) - 0.0) * x_d 
			+ (- r_b2 * sin(beta_c) * cos(thet_c) + 0.0)) * sin(thet_c) 
			- ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * x_d 
			+ (- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c))) * cos(thet_c))
			/ (sin(thet_c) * sin(thet_c));
		Vector3d Dr_c_c, Dr_1_o1, Dr_1_1;
		Matrix<double,3,3> DR_oc_c, DR_1_o1;
		Dr_c_c << 0.0, 0.0 - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
		DR_oc_c << - sin(thet_c), - cos(thet_c), 0.0,
			cos(thet_c), - sin(thet_c), 0.0,
			0.0, 0.0, 0.0;
		DR_1_o1 << - i_1c * sin(thet_1), i_1c * cos(thet_1), 0.0,
			- i_1c * cos(thet_1), - i_1c * sin(thet_1), 0.0,
			0.0, 0.0, 0.0;
		Dr_1_o1 = R_o1_oc * DR_oc_c * r_c_c + R_o1_oc * R_oc_c * Dr_c_c;
		Dr_1_1 = DR_1_o1 * r_1_o1 + R_1_o1 * Dr_1_o1;
		Dfunc(0,0) = i_1c - (Dr_1_o1(1) * r_1_o1(0) - r_1_o1(1) * Dr_1_o1(0))
			/ (r_1_o1(1) * r_1_o1(1) + r_1_o1(0) * r_1_o1(0));
		atan2(r_1_o1(1), r_1_o1(0));
		Dfunc(1,0) = 2.0 * r_1_1(2) * Dr_1_1(2) 
			+ 2.0 * (a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1))) 
			* (0.0 - (r_1_1(0) * Dr_1_1(0) + r_1_1(1) * Dr_1_1(1)) 
			/ sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1)));
		//derivative 2
		Dy_d = - ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) 
			+ 0.0) / sin(thet_c);
		Dr_c_c << - 1.0, 0.0 - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
		Dr_1_o1 = R_o1_oc * R_oc_c * Dr_c_c;
		Dr_1_1 = R_1_o1 * Dr_1_o1;
		Dfunc(0,1) = - (Dr_1_o1(1) * r_1_o1(0) - r_1_o1(1) * Dr_1_o1(0))
			/ (r_1_o1(1) * r_1_o1(1) + r_1_o1(0) * r_1_o1(0));
		Dfunc(1,1) = 2.0 * r_1_1(2) * Dr_1_1(2) 
			+ 2.0 * (a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1))) 
			* (0.0 - (r_1_1(0) * Dr_1_1(0) + r_1_1(1) * Dr_1_1(1)) 
			/ sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1)));
		//iteration
		deltX = Dfunc.fullPivLu().solve(-func);
		if(deltX.norm() < 1.0E-12){
			if(a_h2 - sqrt(r_1_1(0) * r_1_1(0) + r_1_1(1) * r_1_1(1)) < 0.0){
				cout<<"ERROR in DEHW::WORM_CURV_2_CART!"<<endl;
			}
			break;
		}
		x += deltX;
	}
	return 1;
}

long DEHW::WHEE_G2L(Vector3d r_2_2, double &angl_f, double &radi_f, 
	double &R_fmini, double &R_fmaxi){
	double radi_xi = a_h2 - sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
	double coor_zi = r_2_2(2);
	angl_f = atan2(coor_zi, radi_xi);
	radi_f = sqrt(radi_xi * radi_xi + coor_zi * coor_zi);
	double angl_ai = angl_f - asin(offsR_a * sin(angl_f) / R_a(1));
	R_fmini = (R_a(1) * cos(angl_ai) - offsR_a) / cos(angl_f);
	R_fmaxi = R_t(1);
	return 1;
}

long DEHW::WHEE_CURV_2_CART_1(double xi_21, double xi_22, Vector3d &r_2_2, 
	double &thet_c, double &thet_h, long f_lr, double &x_d, double &y_d){
	Matrix<double,2,1> x;
	x << thet_c, thet_h;
	for(long ti = 0; ti < 1000; ti ++){
		Matrix<double,2,1> func, deltX;
		Matrix<double,2,2> Dfunc;
		thet_c = x(0);
		thet_h = x(1);
		//function
		double thet_1 = i_1c * thet_c;
		FSME(thet_1, thet_h, x_d, y_d);
		Matrix<double,3,2> Dr_2_2;
		PD_WHEE_1H2R(x_d, y_d, thet_1, thet_h, r_2_2, Dr_2_2);
		double coorX = a_h2 - sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		func << atan2(r_2_2(2), coorX) - xi_21,
			r_2_2(2) * r_2_2(2) + coorX * coorX - xi_22 * xi_22;
		//derivative 1
		double DcoorX = - (r_2_2(0) * Dr_2_2(0,0) + r_2_2(1) * Dr_2_2(1,0)) 
			/ sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		Dfunc(0,0) = (Dr_2_2(2,0) * coorX - r_2_2(2) * DcoorX) 
			/ (coorX * coorX + r_2_2(2) * r_2_2(2));
		Dfunc(1,0) = 2.0 * r_2_2(2) * Dr_2_2(2,0) + 2.0 * coorX * DcoorX;
		Dfunc.block(0,0,2,1) = i_1c * Dfunc.block(0,0,2,1);
		//derivative 2
		DcoorX = - (r_2_2(0) * Dr_2_2(0,1) + r_2_2(1) * Dr_2_2(1,1)) 
			/ sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		Dfunc(0,1) = (Dr_2_2(2,1) * coorX - r_2_2(2) * DcoorX) 
			/ (coorX * coorX + r_2_2(2) * r_2_2(2));
		Dfunc(1,1) = 2.0 * r_2_2(2) * Dr_2_2(2,1) + 2.0 * coorX * DcoorX;
		//iteration
		deltX = Dfunc.fullPivLu().solve(-func);
		if(deltX.norm() < 1.0E-12){
			if(coorX < 0.0){
				cout<<"ERROR in DEHW::WHEE_CURV_2_CART_1!"<<endl;
			}
			break;
		}
		double rfac = 2.0;
		long rfacFlag = 0;
		while(rfac > 1.0E-10){
			rfac /= 2.0;
			Matrix<double,2,1> x_test = x + rfac * deltX;
			if(x_test(0) < 0.01 * PI || x_test(0) > 0.49 * PI){
				continue;
			}
			double thet_hs, thet_hm;
			SINGULAR_C2H(x(0), thet_hs, thet_hm);
			if((f_lr == 1 && (x_test(1) <= thet_hs + 1.0E-14 || thet_hm - 1.0E-14 <= x_test(1)))
				|| (f_lr == 2 
				&& (x_test(1) <= thet_hm + 1.0E-14 
				|| thet_hs + 2.0 * PI - 1.0E-14 <= x_test(1)))){
				continue;
			}
			double thet_ct = x_test(0);
			double thet_ht = x_test(1);
			//function
			double thet_1t = i_1c * thet_ct;
			double x_dt, y_dt;
			FSME(thet_1t, thet_ht, x_dt, y_dt);
			Vector3d r_2_2t;
			Matrix<double,3,2> Dr_2_2t;
			PD_WHEE_1H2R(x_dt, y_dt, thet_1t, thet_ht, r_2_2t, Dr_2_2t);
			double cooX_t = a_h2 - sqrt(r_2_2t(0) * r_2_2t(0) + r_2_2t(1) * r_2_2t(1));
			Matrix<double,2,1> func_t;
			func_t << atan2(r_2_2t(2), cooX_t) - xi_21,
				r_2_2t(2) * r_2_2t(2) + cooX_t * cooX_t - xi_22 * xi_22;
			if(func_t.norm() < func.norm()){
				x = x_test;
				rfacFlag = 1;
				break;
			}
		}
		if(rfacFlag == 0){
			break;
		}
	}
	return 1;
}

long DEHW::WHEE_CURV_2_CART_2(double xi_21, double xi_22, Vector3d &r_c_c, 
		double &thet_c, double &x_d, double &y_d){
	Matrix<double,2,1> x;
	x << thet_c, x_d;
	for(long ti = 0; ti < 1000; ti ++){
		Matrix<double,2,1> func, deltX;
		Matrix<double,2,2> Dfunc;
		thet_c = x(0);
		x_d = x(1);
		//function
		y_d = - ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * x_d 
			- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c)) 
			/ sin(thet_c);
		r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
		double coorX = a_h2 - sqrt(r_c_c(0) * r_c_c(0) + r_c_c(1) * r_c_c(1));
		func << atan2(r_c_c(2), coorX) - xi_21,
			r_c_c(2) * r_c_c(2) + coorX * coorX - xi_22 * xi_22;
		//derivative 1
		double Dy_d = - (((+ sin(beta_c) * sin(thet_c) - 0.0) * x_d 
			- r_b2 * sin(beta_c) * cos(thet_c) + 0.0) * sin(thet_c) 
			- ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * x_d 
			- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c)) * cos(thet_c))
			/ (sin(thet_c) * sin(thet_c));
		Vector3d Dr_c_c;
		Dr_c_c << - 0.0, 0.0 - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
		double DcoorX = - (r_c_c(0) * Dr_c_c(0) + r_c_c(1) * Dr_c_c(1)) 
			/ sqrt(r_c_c(0) * r_c_c(0) + r_c_c(1) * r_c_c(1));
		Dfunc(0,0) = (Dr_c_c(2) * coorX - r_c_c(2) * DcoorX) 
			/ (coorX * coorX + r_c_c(2) * r_c_c(2));
		Dfunc(1,0) = 2.0 * r_c_c(2) * Dr_c_c(2) + 2.0 * coorX * DcoorX;
		//derivative 2
		Dy_d = - ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * 1.0 
			- 0.0 + 0.0) / sin(thet_c);
		Dr_c_c << - 1.0, - Dy_d * sin(beta_c), Dy_d * cos(beta_c);
		DcoorX = - (r_c_c(0) * Dr_c_c(0) + r_c_c(1) * Dr_c_c(1)) 
			/ sqrt(r_c_c(0) * r_c_c(0) + r_c_c(1) * r_c_c(1));
		Dfunc(0,1) = (Dr_c_c(2) * coorX - r_c_c(2) * DcoorX) 
			/ (coorX * coorX + r_c_c(2) * r_c_c(2));
		Dfunc(1,1) = 2.0 * r_c_c(2) * Dr_c_c(2) + 2.0 * coorX * DcoorX;
		//iteration
		deltX = Dfunc.fullPivLu().solve(-func);
		if(deltX.norm() < 1.0E-12){
			if(coorX < 0.0){
				cout<<"ERROR in DEHW::WHEE_CURV_2_CART_2!"<<endl;
			}
			break;
		}
		double rfac = 2.0;
		long rfacFlag = 0;
		while(rfac > 1.0E-10){
			rfac /= 2.0;
			Matrix<double,2,1> x_test = x + rfac * deltX;
			if(x_test(0) < 0.01 * PI || x_test(0) > 0.49 * PI){
				continue;
			}
			double thet_ct = x_test(0);
			double x_dt = x_test(1);
			//function
			double y_dt = - ((- sin(beta_c) * cos(thet_ct) - i_c1 * cos(beta_c)) * x_dt 
				- r_b2 * sin(beta_c) * sin(thet_ct) + a_1c * sin(beta_c)) 
				/ sin(thet_ct);
			Vector3d r_c_ct;
			r_c_ct << - x_dt, r_b2 - y_dt * sin(beta_c), y_dt * cos(beta_c);
			double cooX_t = a_h2 - sqrt(r_c_ct(0) * r_c_ct(0) + r_c_ct(1) * r_c_ct(1));
			Matrix<double,2,1> func_t;
			func_t << atan2(r_c_ct(2), cooX_t) - xi_21,
				r_c_ct(2) * r_c_ct(2) + cooX_t * cooX_t - xi_22 * xi_22;
			if(func_t.norm() < func.norm()){
				x = x_test;
				rfacFlag = 1;
				break;
			}
		}
		if(rfacFlag == 0){
			break;
		}
	}
	return 1;
}

long DEHW::WHEE_CURV_2_CART_3(double xi_21, double xi_22, Vector3d &r_2_2, 
	double &thet_c, double &thet_h, double xi_11){
	Matrix<double,2,1> x;
	x << thet_c, thet_h;
	for(long ti = 0; ti < 1000; ti ++){
		Matrix<double,2,1> func, deltX;
		Matrix<double,2,2> Dfunc;
		thet_c = x(0);
		thet_h = x(1);
		//function
		Matrix<double,3,2> Dr_2_2;
		WHEE_TRAN(thet_c, thet_h, xi_11, r_2_2, Dr_2_2);
		double coorX = a_h2 - sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		func << atan2(r_2_2(2), coorX) - xi_21,
			r_2_2(2) * r_2_2(2) + coorX * coorX - xi_22 * xi_22;
		//derivative 1
		double DcoorX = - (r_2_2(0) * Dr_2_2(0,0) + r_2_2(1) * Dr_2_2(1,0)) 
			/ sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		Dfunc(0,0) = (Dr_2_2(2,0) * coorX - r_2_2(2) * DcoorX) 
			/ (coorX * coorX + r_2_2(2) * r_2_2(2));
		Dfunc(1,0) = 2.0 * r_2_2(2) * Dr_2_2(2,0) + 2.0 * coorX * DcoorX;
		//derivative 2
		DcoorX = - (r_2_2(0) * Dr_2_2(0,1) + r_2_2(1) * Dr_2_2(1,1)) 
			/ sqrt(r_2_2(0) * r_2_2(0) + r_2_2(1) * r_2_2(1));
		Dfunc(0,1) = (Dr_2_2(2,1) * coorX - r_2_2(2) * DcoorX) 
			/ (coorX * coorX + r_2_2(2) * r_2_2(2));
		Dfunc(1,1) = 2.0 * r_2_2(2) * Dr_2_2(2,1) + 2.0 * coorX * DcoorX;
		//iteration
		deltX = Dfunc.fullPivLu().solve(-func);
		if(deltX.norm() < 1.0E-11){
			if(coorX < 0.0){
				cout<<"ERROR in DEHW::WHEE_CURV_2_CART_3!"<<endl;
			}
			break;
		}
		double rfac = 2.0;
		long rfacFlag = 0;
		while(rfac > 1.0E-10){
			rfac /= 2.0;
			Matrix<double,2,1> x_test = x + rfac * deltX;
			if(x_test(0) < 0.01 * PI || x_test(0) > 0.49 * PI){
				continue;
			}
			double thet_ct = x_test(0);
			double thet_ht = x_test(1);
			//function
			Vector3d r_2_2t;
			WHEE_TRAN(thet_ct, thet_ht, xi_11, r_2_2t, Dr_2_2);
			double cooX_t = a_h2 - sqrt(r_2_2t(0) * r_2_2t(0) + r_2_2t(1) * r_2_2t(1));
			Matrix<double,2,1> func_t;
			func_t << atan2(r_2_2t(2), cooX_t) - xi_21,
				r_2_2t(2) * r_2_2t(2) + cooX_t * cooX_t - xi_22 * xi_22;
			if(func_t.norm() < func.norm()){
				x = x_test;
				rfacFlag = 1;
				break;
			}
		}
		if(rfacFlag == 0){
			break;
		}
	}
	return 1;
}

long DEHW::WHEE_TRAN(double thet_c, double thet_h, double xi_11, 
	Vector3d &r_2_2, Matrix<double,3,2> &Dr_2_2){
	double thet_1 = i_1c * thet_c;
	double thet_2 = i_2h * thet_h;
	//
	double C_1 = (tan(beta_c) * cos(thet_c) + i_c1) * cos(thet_1 - xi_11) 
		+ i_c1 * tan(beta_c) * sin(thet_c) * sin(thet_1 - xi_11) 
		- cos(thet_c) * sin(thet_c) * sin(thet_1 - xi_11);
	double C_2 = i_c1 * r_b2 * sin(thet_c) - i_c1 * a_1c;
	double x_a = - C_2 / C_1;
	double z_a = ((tan(beta_c) * sin(thet_1 - xi_11) + sin(thet_c) * cos(thet_1 - xi_11)) 
		* x_a + r_b2 - a_1c * sin(thet_c)) / cos(thet_c);
	Vector3d r_1_1;
	r_1_1 << x_a * cos(xi_11), - x_a * sin(xi_11), z_a;
	Matrix<double,3,3> R_oh_h;
	R_oh_h << cos(thet_h), - sin(thet_h), 0.0,
		sin(thet_h), cos(thet_h), 0.0,
		0.0, 0.0, 1.0;
	Vector3d r_h_oh = R_oh_h * r_1_1;
	Matrix<double,3,3> R_o2_oh;
	R_o2_oh << 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, -1.0, 0.0;
	Vector3d T_o2_oh;
	T_o2_oh << - a_h2, 0.0, 0.0;
	Vector3d r_h_o2 = R_o2_oh * r_h_oh + T_o2_oh;
	Matrix<double,3,3> R_2_o2;
	R_2_o2 << cos(thet_2), sin(thet_2), 0.0,
		- sin(thet_2), cos(thet_2), 0.0,
		0.0, 0.0, 1.0;
	r_2_2 = R_2_o2 * r_h_o2;
	//derivative relative to thet_c
	double DC_1 = (- tan(beta_c) * sin(thet_c) + 0.0) * cos(thet_1 - xi_11) 
		- (tan(beta_c) * cos(thet_c) + i_c1) * i_1c * sin(thet_1 - xi_11) 
		+ i_c1 * tan(beta_c) * cos(thet_c) * sin(thet_1 - xi_11) 
		+ i_c1 * tan(beta_c) * sin(thet_c) * i_1c * cos(thet_1 - xi_11) 
		+ sin(thet_c) * sin(thet_c) * sin(thet_1 - xi_11) 
		- cos(thet_c) * cos(thet_c) * sin(thet_1 - xi_11) 
		- cos(thet_c) * sin(thet_c) * i_1c * cos(thet_1 - xi_11);
	double DC_2 = i_c1 * r_b2 * cos(thet_c) - 0.0;
	double Dx_a = - (DC_2 * C_1 - C_2 * DC_1) / (C_1 * C_1);
	double Dz_a = (((tan(beta_c) * i_1c * cos(thet_1 - xi_11) 
		+ cos(thet_c) * cos(thet_1 - xi_11) 
		- sin(thet_c) * i_1c * sin(thet_1 - xi_11)) * x_a
		+ (tan(beta_c) * sin(thet_1 - xi_11) + sin(thet_c) * cos(thet_1 - xi_11)) * Dx_a
		+ 0.0 - a_1c * cos(thet_c)) * cos(thet_c) 
		+ ((tan(beta_c) * sin(thet_1 - xi_11) + sin(thet_c) * cos(thet_1 - xi_11)) 
		* x_a + r_b2 - a_1c * sin(thet_c)) * sin(thet_c)) 
		/ (cos(thet_c) * cos(thet_c));
	Vector3d Dr_1_1;
	Dr_1_1 << Dx_a * cos(xi_11), - Dx_a * sin(xi_11), Dz_a;
	Vector3d Dr_h_oh = R_oh_h * Dr_1_1;
	Vector3d Dr_h_o2 = R_o2_oh * Dr_h_oh;
	Dr_2_2.block(0,0,3,1) = R_2_o2 * Dr_h_o2;
	//derivative relative to thet_h
	Matrix<double,3,3> DR_oh_h;
	DR_oh_h << - sin(thet_h), - cos(thet_h), 0.0,
		cos(thet_h), - sin(thet_h), 0.0,
		0.0, 0.0, 0.0;
	Dr_h_oh = DR_oh_h * r_1_1;
	Dr_h_o2 = R_o2_oh * Dr_h_oh;
	Matrix<double,3,3> DR_2_o2;
	DR_2_o2 << - i_2h * sin(thet_2), i_2h * cos(thet_2), 0.0,
		- i_2h * cos(thet_2), - i_2h * sin(thet_2), 0.0,
		0.0, 0.0, 0.0;
	Dr_2_2.block(0,1,3,1) = DR_2_o2 * r_h_o2 + R_2_o2 * Dr_h_o2;
	return 1;
}

long DEHW::WHEE_PHAS(long ti, long tj, long f_ij, Vector3d r_2_2){
	if(fpha(ti,tj) == 0){
		fpha(ti,tj) = f_ij;
		cartCoor(1)(ti,tj) = r_2_2;
	}
	else{
		double phas_1 = atan2(cartCoor(1)(ti,tj)(1), cartCoor(1)(ti,tj)(0));
		if(phas_1 < 0.0){
			phas_1 += 2.0 * PI;
		}
		double phas_2 = atan2(r_2_2(1), r_2_2(0));
		if(phas_2 < 0.0){
			phas_2 += 2.0 *PI;
		}
		if(phas_2 > phas_1){
			fpha(ti,tj) = f_ij;
			cartCoor(1)(ti,tj) = r_2_2;
		}
	}
	return 1;
}

long DEHW::NEW_CONT_ZONE(long f_lr){
	//initial value
	double thet_c, thet_h;
	Vector3d r_2_2;
	double angl_fi, radi_fi, R_fmini, R_fmaxi;
	long numb_c = 1000;
	long numb_h = 10000;
	double thet_cL = 0.01 * PI;
	double thet_cH = 0.49 * PI;
	double epsl_t = 1.0E-8;
	long F_init = 0;
	for(long ti = 0; ti <= numb_c && F_init == 0; ti ++){
		thet_c = thet_cL + (thet_cH - thet_cL) / (double)numb_c * ti;
		double thet_hs, thet_hm;
		SINGULAR_C2H(thet_c, thet_hs, thet_hm);
		double thet_hL, thet_hH;
		if(f_lr == 1){
			thet_hL = thet_hs + epsl_t;
			thet_hH = thet_hm - epsl_t;
		}
		else{
			thet_hL = thet_hm + epsl_t;
			thet_hH = thet_hs + 2.0 * PI - epsl_t;
		}
		if(thet_hL >= thet_hH){
			continue;
		}
		for(long tj = 0; tj <= numb_h && F_init == 0; tj ++){
			if(f_lr == 1){
				thet_h = thet_hH - (thet_hH - thet_hL) / (double)numb_h * tj;
			}
			else{
				thet_h = thet_hL + (thet_hH - thet_hL) / (double)numb_h * tj;
			}
			//????????????????????cillofe, !!!!!the choosing order of thet_h!!!!!
			double thet_1 = i_1c * thet_c;
			double x_d, y_d;
			FSME(thet_1, thet_h, x_d, y_d);
			WHEE_1H2R(x_d, y_d, thet_1, thet_h, r_2_2);
			double kapp_h2N;
			if(CILFOSE_NI(thet_1, thet_h, kapp_h2N) > 0.0){
				WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
				if(-widtAngl <= angl_fi && angl_fi <= widtAngl 
					&& R_fmini <= radi_fi && radi_fi <= R_fmaxi){
					F_init = 1;
				}
			}
		}
	}
	if(F_init == 0){
		cout<<"ERROR in DEHW::NEW_CONT_ZONE"<<f_lr<<"!"<<endl;
		return -1;
	}
	//closest point
	double miniDist = 1.0E20;
	list<INDE_INIT> breaList;
	for(long ti = 0; ti < curvCoor(1).rows() - 1; ti ++){
		for(long tj = 0; tj < curvCoor(1).cols() - 1; tj ++){
			double epsl_x = (curvCoor(1)(ti,tj)(0) - curvCoor(1)(ti+1,tj)(0)) / 4.0;
			double epsl_y = (curvCoor(1)(ti,tj+1)(1) - curvCoor(1)(ti,tj)(1)) / 4.0;
			if(curvCoor(1)(ti+1,tj)(0) - epsl_x <= angl_fi && 
				angl_fi <= curvCoor(1)(ti,tj)(0) + epsl_x && 
				curvCoor(1)(ti,tj)(1) - epsl_y <= radi_fi && 
				radi_fi <= curvCoor(1)(ti,tj+1)(1) + epsl_y){
				breaList.push_back(INDE_INIT(ti, tj, thet_c, thet_h));
				breaList.push_back(INDE_INIT(ti + 1, tj, thet_c, thet_h));
				breaList.push_back(INDE_INIT(ti, tj + 1, thet_c, thet_h));
				breaList.push_back(INDE_INIT(ti + 1, tj + 1, thet_c, thet_h));
			}
			double dist_ij = radi_fi * abs(angl_fi - curvCoor(1)(ti,tj)(0)) 
				+ abs(radi_fi - curvCoor(1)(ti,tj)(1));
			if(dist_ij < miniDist){
				miniDist = dist_ij;
			}
		}
	}
	cout<<"DEHW::NEW_CONT_ZONE"<<f_lr<<", miniDist="<<miniDist
		<<", initial number="<<breaList.size()<<endl;
	//breadth-first search
	double epsl_d = 1.0E-9;
	long coun_w = 0;
	MatrixXi F_sear = MatrixXi::Zero(curvCoor(1).rows(), curvCoor(1).cols());
	while(!breaList.empty()){
		INDE_INIT temp_ij = breaList.front();
		breaList.pop_front();
		if(F_sear(temp_ij.ti, temp_ij.tj) == 1){
			continue;
		}
		//
		double x_d, y_d;
		WHEE_CURV_2_CART_1(curvCoor(1)(temp_ij.ti, temp_ij.tj)(0), 
			curvCoor(1)(temp_ij.ti, temp_ij.tj)(1), 
			r_2_2, temp_ij.init_1, temp_ij.init_2, f_lr, x_d, y_d);
		//
		WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
		double dist_ij = 
			radi_fi * abs(angl_fi - curvCoor(1)(temp_ij.ti, temp_ij.tj)(0)) 
			+ abs(radi_fi - curvCoor(1)(temp_ij.ti, temp_ij.tj)(1));
		if(dist_ij < epsl_d){
			F_sear(temp_ij.ti, temp_ij.tj) = 1;
			//
			thet_c = temp_ij.init_1;
			double thet_1 = i_1c * thet_c;
			Vector3d r_c_c, T_o1_oc, r_1_o1;
			Matrix<double,3,3> R_oc_c, R_o1_oc, R_1_o1;
			r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
			R_oc_c << cos(thet_c), - sin(thet_c), 0.0,
				sin(thet_c), cos(thet_c), 0.0,
				0.0, 0.0, 1.0;
			R_o1_oc << 1.0, 0.0, 0.0,
				0.0, 0.0, -1.0,
				0.0, 1.0, 0.0;
			R_1_o1 << cos(thet_1), sin(thet_1), 0.0,
				- sin(thet_1), cos(thet_1), 0.0,
				0.0, 0.0, 1.0;
			T_o1_oc << a_1c, 0.0, 0.0;
			r_1_o1 = R_o1_oc * R_oc_c * r_c_c + T_o1_oc;
			double woxi_11 =  thet_1 - atan2(r_1_o1(1), r_1_o1(0));
			if(curvCoor(0)(0,0)(0) - 1.0E-12 <= woxi_11 
				&& woxi_11 <= curvCoor(0)(curvCoor(0).rows() - 1,0)(0) + 1.0E-12){
				WHEE_PHAS(temp_ij.ti, temp_ij.tj, f_lr, r_2_2);
			}
			//
			Matrix<long,8,2> inde;
			inde << temp_ij.ti - 1, temp_ij.tj - 1, 
				temp_ij.ti, temp_ij.tj - 1, 
				temp_ij.ti + 1, temp_ij.tj - 1, 
				temp_ij.ti - 1, temp_ij.tj, 
				temp_ij.ti + 1, temp_ij.tj, 
				temp_ij.ti - 1, temp_ij.tj + 1, 
				temp_ij.ti, temp_ij.tj + 1, 
				temp_ij.ti + 1, temp_ij.tj + 1;
			for(long ti = 0; ti < inde.rows(); ti ++){
				if(inde(ti,0) >= 0 && inde(ti,1) >= 0 
					&& inde(ti,0) < curvCoor(1).rows() && inde(ti,1) < curvCoor(1).cols() 
					&& F_sear(inde(ti,0), inde(ti,1)) == 0){
					breaList.push_back(INDE_INIT(inde(ti,0), inde(ti,1), 
						temp_ij.init_1, temp_ij.init_2));
				}
			}
		}
		coun_w ++;
		if(coun_w % 4000 == 0){
			cout<<coun_w<<endl;
		}
	}
	return 1;
}

long DEHW::FORMER_CONT_ZONE(){
	//initial value
	double thet_c, x_d;
	Vector3d r_c_c;
	double angl_fi, radi_fi, R_fmini, R_fmaxi;
	long numb_c = 1000;
	double thet_cL = 0.01 * PI;
	double thet_cH = 0.49 * PI;
	long numb_d = 10000;
	double x_dL = - 10.0 * a_1c;
	double x_dH = 10.0 * a_1c;
	long F_init = 0;
	for(long ti = 0; ti <= numb_c && F_init == 0; ti ++){
		thet_c = thet_cL + (thet_cH - thet_cL) / (double)numb_c * ti;
		for(long tj = 1; tj < numb_d && F_init == 0; tj ++){
			x_d = x_dL + (x_dH - x_dL) / (double)numb_d * tj;
			double y_d = - ((- sin(beta_c) * cos(thet_c) - i_c1 * cos(beta_c)) * x_d 
				- r_b2 * sin(beta_c) * sin(thet_c) + a_1c * sin(beta_c)) / sin(thet_c);
			r_c_c << - x_d, r_b2 - y_d * sin(beta_c), y_d * cos(beta_c);
			WHEE_G2L(r_c_c, angl_fi, radi_fi, R_fmini, R_fmaxi);
			if(-widtAngl <= angl_fi && angl_fi <= widtAngl 
				&& R_fmini <= radi_fi && radi_fi <= R_fmaxi){
				F_init = 1;
			}
		}
	}
	if(F_init == 0){
		cout<<"ERROR in DEHW::FORMER_CONT_ZONE!"<<endl;
		return -1;
	}
	//closest point
	double miniDist = 1.0E20;
	list<INDE_INIT> breaList;
	for(long ti = 0; ti < curvCoor(1).rows() - 1; ti ++){
		for(long tj = 0; tj < curvCoor(1).cols() - 1; tj ++){
			double epsl_x = (curvCoor(1)(ti,tj)(0) - curvCoor(1)(ti+1,tj)(0)) / 4.0;
			double epsl_y = (curvCoor(1)(ti,tj+1)(1) - curvCoor(1)(ti,tj)(1)) / 4.0;
			if(curvCoor(1)(ti+1,tj)(0) - epsl_x <= angl_fi && 
				angl_fi <= curvCoor(1)(ti,tj)(0) + epsl_x && 
				curvCoor(1)(ti,tj)(1) - epsl_y <= radi_fi && 
				radi_fi <= curvCoor(1)(ti,tj+1)(1) + epsl_y){
				breaList.push_back(INDE_INIT(ti, tj, thet_c, x_d));
				breaList.push_back(INDE_INIT(ti + 1, tj, thet_c, x_d));
				breaList.push_back(INDE_INIT(ti, tj + 1, thet_c, x_d));
				breaList.push_back(INDE_INIT(ti + 1, tj + 1, thet_c, x_d));
			}
			double dist_ij = radi_fi * abs(angl_fi - curvCoor(1)(ti,tj)(0)) 
				+ abs(radi_fi - curvCoor(1)(ti,tj)(1));
			if(dist_ij < miniDist){
				miniDist = dist_ij;
			}
		}
	}
	cout<<"DEHW::FORMER_CONT_ZONE, miniDist="<<miniDist
		<<", initial number="<<breaList.size()<<endl;
	//breadth-first search
	double epsl_d = 1.0E-9;
	long coun_w = 0;
	MatrixXi F_sear = MatrixXi::Zero(curvCoor(1).rows(), curvCoor(1).cols());
	while(!breaList.empty()){
		INDE_INIT temp_ij = breaList.front();
		breaList.pop_front();
		if(F_sear(temp_ij.ti, temp_ij.tj) == 1){
			continue;
		}
		//
		double y_d;
		WHEE_CURV_2_CART_2(curvCoor(1)(temp_ij.ti, temp_ij.tj)(0), 
			curvCoor(1)(temp_ij.ti, temp_ij.tj)(1), 
			r_c_c, temp_ij.init_1, temp_ij.init_2, y_d);
		//
		WHEE_G2L(r_c_c, angl_fi, radi_fi, R_fmini, R_fmaxi);
		double dist_ij = 
			radi_fi * abs(angl_fi - curvCoor(1)(temp_ij.ti, temp_ij.tj)(0)) 
			+ abs(radi_fi - curvCoor(1)(temp_ij.ti, temp_ij.tj)(1));
		if(dist_ij < epsl_d){
			F_sear(temp_ij.ti, temp_ij.tj) = 1;
			//
			thet_c = temp_ij.init_1;
			double thet_1 = i_1c * thet_c;
			Vector3d T_o1_oc, r_1_o1;
			Matrix<double,3,3> R_oc_c, R_o1_oc, R_1_o1;
			R_oc_c << cos(thet_c), - sin(thet_c), 0.0,
				sin(thet_c), cos(thet_c), 0.0,
				0.0, 0.0, 1.0;
			R_o1_oc << 1.0, 0.0, 0.0,
				0.0, 0.0, -1.0,
				0.0, 1.0, 0.0;
			R_1_o1 << cos(thet_1), sin(thet_1), 0.0,
				- sin(thet_1), cos(thet_1), 0.0,
				0.0, 0.0, 1.0;
			T_o1_oc << a_1c, 0.0, 0.0;
			r_1_o1 = R_o1_oc * R_oc_c * r_c_c + T_o1_oc;
			double woxi_11 =  thet_1 - atan2(r_1_o1(1), r_1_o1(0));
			if(curvCoor(0)(0,0)(0) - 1.0E-12 <= woxi_11 
				&& woxi_11 <= curvCoor(0)(curvCoor(0).rows() - 1,0)(0) + 1.0E-12){
				WHEE_PHAS(temp_ij.ti, temp_ij.tj, 3, r_c_c);
			}
			//
			Matrix<long,8,2> inde;
			inde << temp_ij.ti - 1, temp_ij.tj - 1, 
				temp_ij.ti, temp_ij.tj - 1, 
				temp_ij.ti + 1, temp_ij.tj - 1, 
				temp_ij.ti - 1, temp_ij.tj, 
				temp_ij.ti + 1, temp_ij.tj, 
				temp_ij.ti - 1, temp_ij.tj + 1, 
				temp_ij.ti, temp_ij.tj + 1, 
				temp_ij.ti + 1, temp_ij.tj + 1;
			for(long ti = 0; ti < inde.rows(); ti ++){
				if(inde(ti,0) >= 0 && inde(ti,1) >= 0 
					&& inde(ti,0) < curvCoor(1).rows() && inde(ti,1) < curvCoor(1).cols() 
					&& F_sear(inde(ti,0), inde(ti,1)) == 0){
					breaList.push_back(INDE_INIT(inde(ti,0), inde(ti,1), 
						temp_ij.init_1, temp_ij.init_2));
				}
			}
		}
		coun_w ++;
		if(coun_w % 4000 == 0){
			cout<<coun_w<<endl;
		}
	}
	return 1;
}

long DEHW::TRANSITION_ZONE(long f_hr){
	double xi_11;
	if(f_hr == 1){
		xi_11 = curvCoor(0)(0,0)(0);//head
	}
	else{
		xi_11 = curvCoor(0)(curvCoor(0).rows() - 1,0)(0);//rear
	}
	//initial value
	double thet_c, thet_h;
	Vector3d r_2_2, r_1_1;
	double angl_fi, radi_fi, R_fmini, R_fmaxi;
	long numb_c = 1000;
	double thet_cL, thet_cH;
	WORM_CURV_2_CART(xi_11, a_h2 - d_f(0) / 2.0, r_1_1, thet_cL);
	WORM_CURV_2_CART(xi_11, d_f(1) / 2.0, r_1_1, thet_cH);
	thet_h = xi_11;
	long F_init = 0;
	for(long ti = 0; ti <= numb_c && F_init == 0; ti ++){
		thet_c = thet_cL + (thet_cH - thet_cL) / (double)numb_c * ti;
		Matrix<double,3,2> Dr_2_2;
		WHEE_TRAN(thet_c, thet_h, xi_11, r_2_2, Dr_2_2);
		WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
		if(-widtAngl <= angl_fi && angl_fi <= widtAngl 
			&& R_fmini <= radi_fi && radi_fi <= R_fmaxi){
			F_init = 1;
		}
	}
	if(F_init == 0){
		cout<<"ERROR in DEHW::TRANSITION_ZONE"<<f_hr<<"!"<<endl;
		return -1;
	}
	//closest point
	double miniDist = 1.0E20;
	list<INDE_INIT> breaList;
	for(long ti = 0; ti < curvCoor(1).rows() - 1; ti ++){
		for(long tj = 0; tj < curvCoor(1).cols() - 1; tj ++){
			double epsl_x = (curvCoor(1)(ti,tj)(0) - curvCoor(1)(ti+1,tj)(0)) / 4.0;
			double epsl_y = (curvCoor(1)(ti,tj+1)(1) - curvCoor(1)(ti,tj)(1)) / 4.0;
			if(curvCoor(1)(ti+1,tj)(0) - epsl_x <= angl_fi && 
				angl_fi <= curvCoor(1)(ti,tj)(0) + epsl_x && 
				curvCoor(1)(ti,tj)(1) - epsl_y <= radi_fi && 
				radi_fi <= curvCoor(1)(ti,tj+1)(1) + epsl_y){
				breaList.push_back(INDE_INIT(ti, tj, thet_c, thet_h));
				breaList.push_back(INDE_INIT(ti + 1, tj, thet_c, thet_h));
				breaList.push_back(INDE_INIT(ti, tj + 1, thet_c, thet_h));
				breaList.push_back(INDE_INIT(ti + 1, tj + 1, thet_c, thet_h));
			}
			double dist_ij = radi_fi * abs(angl_fi - curvCoor(1)(ti,tj)(0)) 
				+ abs(radi_fi - curvCoor(1)(ti,tj)(1));
			if(dist_ij < miniDist){
				miniDist = dist_ij;
			}
		}
	}
	cout<<"DEHW::TRANSITION_ZONE"<<f_hr<<", miniDist="<<miniDist
		<<", initial number="<<breaList.size()<<endl;
	//breadth-first search
	double epsl_d = 1.0E-9;
	long coun_w = 0;
	MatrixXi F_sear = MatrixXi::Zero(curvCoor(1).rows(), curvCoor(1).cols());
	while(!breaList.empty()){
		INDE_INIT temp_ij = breaList.front();
		breaList.pop_front();
		if(F_sear(temp_ij.ti, temp_ij.tj) == 1){
			continue;
		}
		//
		WHEE_CURV_2_CART_3(curvCoor(1)(temp_ij.ti, temp_ij.tj)(0), 
			curvCoor(1)(temp_ij.ti, temp_ij.tj)(1), 
			r_2_2, temp_ij.init_1, temp_ij.init_2, xi_11);
		//
		WHEE_G2L(r_2_2, angl_fi, radi_fi, R_fmini, R_fmaxi);
		double dist_ij = 
			radi_fi * abs(angl_fi - curvCoor(1)(temp_ij.ti, temp_ij.tj)(0)) 
			+ abs(radi_fi - curvCoor(1)(temp_ij.ti, temp_ij.tj)(1));
		if(dist_ij < epsl_d){
			F_sear(temp_ij.ti, temp_ij.tj) = 1;
			WHEE_PHAS(temp_ij.ti, temp_ij.tj, 3 + f_hr, r_2_2);
			Matrix<long,8,2> inde;
			inde << temp_ij.ti - 1, temp_ij.tj - 1, 
				temp_ij.ti, temp_ij.tj - 1, 
				temp_ij.ti + 1, temp_ij.tj - 1, 
				temp_ij.ti - 1, temp_ij.tj, 
				temp_ij.ti + 1, temp_ij.tj, 
				temp_ij.ti - 1, temp_ij.tj + 1, 
				temp_ij.ti, temp_ij.tj + 1, 
				temp_ij.ti + 1, temp_ij.tj + 1;
			for(long ti = 0; ti < inde.rows(); ti ++){
				if(inde(ti,0) >= 0 && inde(ti,1) >= 0 
					&& inde(ti,0) < curvCoor(1).rows() && inde(ti,1) < curvCoor(1).cols() 
					&& F_sear(inde(ti,0), inde(ti,1)) == 0){
					breaList.push_back(INDE_INIT(inde(ti,0), inde(ti,1), 
						temp_ij.init_1, temp_ij.init_2));
				}
			}
		}
		coun_w ++;
		if(coun_w % 4000 == 0){
			cout<<coun_w<<endl;
		}
	}
	return 1;
}

long DEHW::WORM_RELI(Vector3d &tempCoor, long ti, long tj){
	long reliLeng = 40;
	if(tj > curvCoor(0).cols() - reliLeng 
		|| ti < reliLeng || ti > curvCoor(0).rows() - reliLeng){
		Vector3d tempXYZ = tempCoor;
		Matrix<double,2,1> reliAmou;
		reliAmou << 14.0E-6, 18.0E-6;//0 - tip, 1 - end
		double reliExpo = 3.0;
		double tempReli = 0.0;
		if(tj > curvCoor(0).cols() - reliLeng 
			&& ti >= reliLeng && ti <= curvCoor(0).rows() - reliLeng){
			tempReli = 
				pow((tj - (curvCoor(0).cols() - reliLeng)) / (double)reliLeng, reliExpo) 
				* reliAmou(0);
		}
		else if(tj > curvCoor(0).cols() - reliLeng && ti < reliLeng){
			double tempAngl = 
				atan((tj - (curvCoor(0).cols() - reliLeng)) / (double)(reliLeng - ti));
			double tempRati = tempAngl / (PI / 2.0);
			double maxiAmou = reliAmou(1) 
				+ (- 1.0 + cos(tempRati * PI)) * (reliAmou(1) - reliAmou(0)) / 2.0;
			double tempRadi = sqrt(pow(tj - (curvCoor(0).cols() - reliLeng), 2.0) 
				+ pow(reliLeng - ti, 2.0));
			tempReli = pow(tempRadi / (double)reliLeng, reliExpo) * maxiAmou;
		}
		else if(tj > curvCoor(0).cols() - reliLeng && ti > curvCoor(0).rows() - reliLeng){
			double tempAngl = atan((tj - (curvCoor(0).cols() - reliLeng)) 
				/ (double)(ti - (curvCoor(0).rows() - reliLeng)));
			double tempRati = tempAngl / (PI / 2.0);
			double maxiAmou = reliAmou(1) 
				+ (- 1.0 + cos(tempRati * PI)) * (reliAmou(1) - reliAmou(0)) / 2.0;
			double tempRadi = sqrt(pow(tj - (curvCoor(0).cols() - reliLeng), 2.0) 
				+ pow(ti - (curvCoor(0).rows() - reliLeng), 2.0));
			tempReli = pow(tempRadi / (double)reliLeng, reliExpo) * maxiAmou;
		}
		else if(tj <= curvCoor(0).cols() - reliLeng && ti < reliLeng){
			tempReli = pow((reliLeng - ti) / (double)reliLeng, reliExpo) * reliAmou(1);
		}
		else if(tj <= curvCoor(0).cols() - reliLeng && ti > curvCoor(0).rows() - reliLeng){
			tempReli = 
				pow((ti - (curvCoor(0).rows() - reliLeng)) / (double)reliLeng, reliExpo) 
				* reliAmou(1);
		}
		if(abs(tempReli) > 1.0E-12){
			double tempRadi_0 = sqrt(pow(tempXYZ(0), 2.0) + pow(tempXYZ(1), 2.0));
			double tempRadi = a_h2 - tempRadi_0;
			tempRadi = sqrt(tempRadi * tempRadi + tempXYZ(2) * tempXYZ(2));
			double tempAngl = tempReli / tempRadi;
			double tempThet_0 = asin(tempXYZ(2) / tempRadi);
			double tempThet_1 = tempThet_0 + tempAngl;
			double tempRadi_1 = a_h2 - tempRadi * cos(tempThet_1);
			double tempFact = tempRadi_1 / tempRadi_0;
			tempXYZ(0) = tempFact * tempXYZ(0);
			tempXYZ(1) = tempFact * tempXYZ(1);
			tempXYZ(2) = tempXYZ(2) + tempRadi * (sin(tempThet_1) - sin(tempThet_0));
			tempCoor = tempXYZ;
		}
	}
	return 1;
}

long DEHW::WHEE_RELI(Vector3d &tempCoor, long ti, long tj){
	double reliLeng = 40;
	if(tj < reliLeng || ti < reliLeng || ti > curvCoor(1).rows() - reliLeng){
		Vector3d tempXYZ = tempCoor;
		double reliExpo = 3.0;
		Matrix<double,2,1> reliAmou;
		reliAmou << 12.0E-6, 16.0E-6;//0 - tip, 1 - end
		//
		double tempReli = 0.0;
		if(ti < reliLeng){
			if(tj >= reliLeng){
				tempReli = pow((reliLeng - ti) / (double)reliLeng, reliExpo) * reliAmou(1);
			}
			else{
				double tempAngl = atan((reliLeng - tj) / (double)(reliLeng - ti));
				double tempRati = tempAngl / (PI / 2.0);
				double maxiAmou = reliAmou(1) 
					+ (- 1.0 + cos(tempRati * PI)) * (reliAmou(1) - reliAmou(0)) / 2.0;
				double tempRadi = sqrt(pow((reliLeng - tj), 2.0) + pow((reliLeng - ti), 2.0));
				tempReli = pow((tempRadi / reliLeng), reliExpo) * maxiAmou;
			}
		}
		else if(ti > curvCoor(1).rows() - reliLeng){
			if(tj >= reliLeng){
				tempReli = 
					pow((ti - (curvCoor(1).rows() - reliLeng)) / (double)reliLeng, reliExpo) 
					* reliAmou(1);
			}
			else{
				double tempAngl = atan((reliLeng - tj) 
					/ (double)(ti - (curvCoor(1).rows() - reliLeng)));
				double tempRati = tempAngl / (PI / 2.0);
				double maxiAmou = reliAmou(1) 
					+ (- 1.0 + cos(tempRati * PI)) * (reliAmou(1) - reliAmou(0)) / 2.0;
				double tempRadi = sqrt(pow((reliLeng - tj), 2.0) 
					+ pow((ti - (curvCoor(1).rows() - reliLeng)), 2.0));
				tempReli = pow((tempRadi / reliLeng), reliExpo) * maxiAmou;
			}
		}
		else if(tj < reliLeng){
			tempReli = pow((reliLeng - tj) / (double)reliLeng, reliExpo) * reliAmou(0);
		}
		if(abs(tempReli) > 1.0E-12){
			double tempAngl = tempReli 
				/ sqrt(tempXYZ(0) * tempXYZ(0) + tempXYZ(1) * tempXYZ(1));
			Matrix<double,3,3> tempMatr;
			tempMatr << cos(tempAngl), -sin(tempAngl), 0.0,
				sin(tempAngl), cos(tempAngl), 0.0,
				0.0, 0.0, 1.0;
			tempXYZ = tempMatr * tempXYZ;
			tempCoor = tempXYZ;
		}
	}
	return 1;
}

long DEHW::WORM_GRID(){
	OUTPUT_TIME("WORM_GRID");
	//profile
	wormProf.resize(gridNumb(0,5) + 1, gridNumb(0,3) + 1);
	profCurv.resize(gridNumb(0,5) + 1, 1);
	for(long ti = 0; ti <= gridNumb(0,5); ti ++){
		for(long tj = 0; tj <= gridNumb(0,3); tj ++){
			long id_ti = ti * pow2[refiLeve];
			long id_tj = tj * pow2[refiLeve];
			wormProf(ti,tj).block(0,0,3,1) = cartCoor(0)(id_ti, id_tj);
			profCurv(ti) = - curvCoor(0)(id_ti, id_tj)(0);
			Matrix<double,3,1> tempXYZ = cartCoor(0)(id_ti, id_tj);
			double tempAngl = wormCurv(1);
			Matrix<double,3,3> tempMatr;
			tempMatr << cos(tempAngl), -sin(tempAngl), 0.0,
				sin(tempAngl), cos(tempAngl), 0.0,
				0.0, 0.0, 1.0;
			Matrix<double,3,1> tempXYZ_1 = tempMatr * tempXYZ;
			tempXYZ << tempXYZ_1(0), -tempXYZ_1(1), -tempXYZ(2);
			wormProf(gridNumb(0,5) - ti, tj).block(0,1,3,1) = tempMatr.transpose() * tempXYZ;
		}
	}
	//points
	Matrix3d tempRota_0;
	tempRota_0 << cos(analAngl(0)), -sin(analAngl(0)), 0.0, 
		sin(analAngl(0)), cos(analAngl(0)), 0.0, 
		0.0, 0.0, 1.0;
	Matrix3d tempRota_1;
	tempRota_1 << 1.0, 0.0, 0.0, 
		0.0, 0.0, 1.0, 
		0.0, -1.0, 0.0;
	tempRota_1 = tempRota_1 * tempRota_0;
	Vector3d tempTran;
	tempTran << - (a_h2 + centErro), 0.0, 0.0;
	Matrix<Matrix<Matrix<long,Dynamic,Dynamic>,Dynamic,Dynamic>,Dynamic,Dynamic> blocPoin;
	blocPoin.resize(z(0),1);
	blocPoin(0).resize(gridNumb(0,5) + 1, 4);
	for(long ti = 0; ti <= gridNumb(0,5); ti ++){
		blocPoin(0)(ti,0).resize(gridNumb(0,1) + 1, gridNumb(0,0) + 1);
		blocPoin(0)(ti,1).resize(gridNumb(0,2) + 1, gridNumb(0,0) / 2 + 1);
		blocPoin(0)(ti,2).resize(gridNumb(0,3) + 1, 2 * gridNumb(0,2) + 1);
		blocPoin(0)(ti,3).resize(gridNumb(0,2) + 1, gridNumb(0,0) / 2 + 1);
		//
		Matrix<double,Dynamic,Dynamic> rootProf_1, rootProf_2;
		WORM_ROOT(ti, 1, rootProf_1);
		WORM_ROOT(ti, -1, rootProf_2);
		Matrix<double,2,1> profRadi;
		profRadi << sqrt(rootProf_1(0,0) * rootProf_1(0,0) + rootProf_1(1,0) * rootProf_1(1,0)),
			sqrt(rootProf_2(0,0) * rootProf_2(0,0) + rootProf_2(1,0) * rootProf_2(1,0));
		Matrix<double,2,1> tranRadi;
		tranRadi(0) = profRadi(0) - PI / 4.0 * m_t;
		tranRadi(1) = profRadi(1) - PI / 4.0 * m_t;
		Matrix<double,3,4> blocCoor;
		blocCoor << 
			rootProf_1(0,0) * inneRadi(0) / profRadi(0), 
			rootProf_2(0,0) * inneRadi(0) / profRadi(1), 
			rootProf_2(0,0) * tranRadi(1) / profRadi(1), 
			rootProf_1(0,0) * tranRadi(0) / profRadi(0),
			rootProf_1(1,0) * inneRadi(0) / profRadi(0), 
			rootProf_2(1,0) * inneRadi(0) / profRadi(1), 
			rootProf_2(1,0) * tranRadi(1) / profRadi(1), 
			rootProf_1(1,0) * tranRadi(0) / profRadi(0),
			rootProf_1(2,0), rootProf_2(2,0), rootProf_2(2,0), rootProf_1(2,0);
		//0
		COORDINATE tempCoor;
		for(long tj = 0; tj <= gridNumb(0,1); tj ++){
			for(long tk = 0; tk <= gridNumb(0,0); tk ++){
				Matrix<double,3,2> downUp;
				downUp.block(0,0,3,1) = 
					(1.0 - (double)tk / gridNumb(0,0)) * blocCoor.block(0,0,3,1).eval() 
					+ (double)tk / gridNumb(0,0) * blocCoor.block(0,1,3,1).eval();
				downUp.block(0,1,3,1) = 
					(1.0 - (double)tk / gridNumb(0,0)) * blocCoor.block(0,3,3,1).eval() 
					+ (double)tk / gridNumb(0,0) * blocCoor.block(0,2,3,1).eval();
				Vector3d tempMatr = 
					(1.0 - (double)tj / gridNumb(0,1)) * downUp.block(0,0,3,1).eval() 
					+ (double)tj / gridNumb(0,1) * downUp.block(0,1,3,1).eval();
				tempCoor.COPY(tempRota_1 * tempMatr + tempTran);
				blocPoin(0)(ti,0)(tj,tk) = vem(0).TRY_ADD_POINT(tempCoor);
			}
		}
		//1
		Matrix<double,Dynamic,Dynamic> lineCoor;
		lineCoor.resize(3,2);
		lineCoor.block(0,0,3,1) = blocCoor.block(0,3,3,1).eval();
		lineCoor.block(0,1,3,1) = 0.5 * blocCoor.block(0,3,3,1).eval() 
			+ 0.5 * blocCoor.block(0,2,3,1).eval();
		for(long tj = 0; tj <= gridNumb(0,2); tj ++){
			for(long tk = 0; tk <= gridNumb(0,0) / 2; tk ++){
				Matrix<double,3,2> downUp;
				downUp.block(0,0,3,1) = 
					(1.0 - (double)tk / (gridNumb(0,0) / 2.0)) * lineCoor.block(0,0,3,1).eval() 
					+ (double)tk / (gridNumb(0,0)/2.0) * lineCoor.block(0,1,3,1).eval();
				downUp.block(0,1,3,1) = rootProf_1.block(0,tk,3,1);
				Vector3d tempMatr = 
					(1.0 - (double)tj / gridNumb(0,2)) * downUp.block(0,0,3,1).eval() 
					+ (double)tj / gridNumb(0,2) * downUp.block(0,1,3,1).eval();
				tempCoor.COPY(tempRota_1 * tempMatr + tempTran);
				blocPoin(0)(ti,1)(tj,tk) = vem(0).TRY_ADD_POINT(tempCoor);
			}
		}
		//2
		lineCoor.block(0,0,3,1) = lineCoor.block(0,1,3,1).eval();
		lineCoor.block(0,1,3,1) = 0.5 * wormProf(ti,gridNumb(0,3)).block(0,0,3,1) 
			+ 0.5 * wormProf(ti,gridNumb(0,3)).block(0,1,3,1);
		for(long tj = 0; tj <= gridNumb(0,3); tj ++){
			for(long tk = 0; tk <= gridNumb(0,2); tk ++){
				Matrix<double,3,2> downUp;
				downUp.block(0,0,3,1) = wormProf(ti,tj).block(0,0,3,1);
				downUp.block(0,1,3,1) = 
					(1.0 - (double)tj / gridNumb(0,3)) * lineCoor.block(0,0,3,1).eval() 
					+ (double)tj / gridNumb(0,3) * lineCoor.block(0,1,3,1).eval();
				Vector3d tempMatr = 
					(1.0 - (double)tk / gridNumb(0,2)) * downUp.block(0,0,3,1).eval() 
					+ (double)tk / gridNumb(0,2) * downUp.block(0,1,3,1).eval();
				tempCoor.COPY(tempRota_1 * tempMatr + tempTran);
				blocPoin(0)(ti,2)(tj,tk) = vem(0).TRY_ADD_POINT(tempCoor);
			}
			for(long tk = 1; tk <= gridNumb(0,2); tk ++){
				Matrix<double,3,2> downUp;
				downUp.block(0,0,3,1) = 
					(1.0 - (double)tj / gridNumb(0,3)) * lineCoor.block(0,0,3,1).eval() 
					+ (double)tj / gridNumb(0,3) * lineCoor.block(0,1,3,1).eval();
				downUp.block(0,1,3,1) = wormProf(ti,tj).block(0,1,3,1);
				Vector3d tempMatr = 
					(1.0 - (double)tk / gridNumb(0,2)) * downUp.block(0,0,3,1).eval() 
					+ (double)tk / gridNumb(0,2) * downUp.block(0,1,3,1).eval();
				tempCoor.COPY(tempRota_1 * tempMatr + tempTran);
				blocPoin(0)(ti,2)(tj, gridNumb(0,2) + tk) = vem(0).TRY_ADD_POINT(tempCoor);
			}
		}
		//3
		lineCoor.block(0,1,3,1) = blocCoor.block(0,2,3,1).eval();
		for(long tj = 0; tj <= gridNumb(0,2); tj ++){
			for(long tk = 0; tk <= gridNumb(0,0) / 2; tk ++){
				Matrix<double,3,2> downUp;
				downUp.block(0,0,3,1) = 
					(1.0 - (double)tk / (gridNumb(0,0) / 2.0)) * lineCoor.block(0,0,3,1).eval() 
					+ (double)tk / (gridNumb(0,0)/2.0) * lineCoor.block(0,1,3,1).eval();
				downUp.block(0,1,3,1) = rootProf_2.block(0,gridNumb(0,0) / 2 - tk,3,1);
				Vector3d tempMatr = 
					(1.0 - (double)tj / gridNumb(0,2)) * downUp.block(0,0,3,1).eval() 
					+ (double)tj / gridNumb(0,2) * downUp.block(0,1,3,1).eval();
				tempCoor.COPY(tempRota_1 * tempMatr + tempTran);
				blocPoin(0)(ti,3)(tj,tk) = vem(0).TRY_ADD_POINT(tempCoor);
			}
		}
	}
	for(long ti = 1; ti < z(0); ti ++){
		blocPoin(ti).resize(blocPoin(0).rows(), blocPoin(0).cols());
		double tempAngl = 2.0 * PI / z(0) * ti;
		Matrix<double,3,3> tempMatr;
		tempMatr << cos(tempAngl), 0.0, sin(tempAngl),
			0.0, 1.0, 0.0,
			-sin(tempAngl), 0.0, cos(tempAngl);
		for(long tj = 0; tj < blocPoin(0).rows(); tj ++){
			for(long tk = 0; tk < blocPoin(0).cols(); tk ++){
				blocPoin(ti)(tj,tk).resize(blocPoin(0)(tj,tk).rows(), blocPoin(0)(tj,tk).cols());
				for(long tm = 0; tm < blocPoin(0)(tj,tk).rows(); tm ++){
					for(long tn = 0; tn < blocPoin(0)(tj,tk).cols(); tn ++){
						Vector3d tempVect = vem(0).poinCoor[blocPoin(0)(tj,tk)(tm,tn)];
						COORDINATE tempCoor;
						tempCoor.COPY(tempMatr * (tempVect - tempTran) + tempTran);
						blocPoin(ti)(tj,tk)(tm,tn) = vem(0).TRY_ADD_POINT(tempCoor);
					}
				}
			}
		}
	}
	//volumes
	for(long ti = 0; ti < z(0); ti ++){
		for(long tj = 0; tj < blocPoin(ti).rows() - 1; tj ++){
			if(tj % 30 == 0){
				cout << "Volume layer " << tj << "/" << gridNumb(0,5) 
					<<"~"<< ti <<"/"<< z(0) << endl;
			}
			for(long tk = 0; tk < blocPoin(ti).cols(); tk ++){
				for(long tm = 0; tm < blocPoin(ti)(tj,tk).rows() - 1; tm ++){
					for(long tn = 0; tn < blocPoin(ti)(tj,tk).cols() - 1; tn ++){
						Matrix<long,8,1> tempCorn;
						tempCorn << blocPoin(ti)(tj,tk)(tm,tn + 1),
							blocPoin(ti)(tj,tk)(tm,tn),
							blocPoin(ti)(tj,tk)(tm + 1,tn),
							blocPoin(ti)(tj,tk)(tm + 1,tn + 1),
							blocPoin(ti)(tj + 1,tk)(tm,tn + 1),
							blocPoin(ti)(tj + 1,tk)(tm,tn),
							blocPoin(ti)(tj + 1,tk)(tm + 1,tn),
							blocPoin(ti)(tj + 1,tk)(tm + 1,tn + 1);
						long tempIden = vem(0).ADD_BLOCK_FROM_POINT(tempCorn);
						if(tk == 2 && tn == 0){
							MATRIX tempMatr_0;
							tempMatr_0.resize(4,1);
							tempMatr_0 << blocPoin(ti)(tj,tk)(tm,tn), 
								blocPoin(ti)(tj + 1,tk)(tm,tn), 
								blocPoin(ti)(tj + 1,tk)(tm + 1,tn),
								blocPoin(ti)(tj,tk)(tm + 1,tn);
							segmVolu(0).push_back(make_pair(tempMatr_0, -tempIden));
							Matrix<long,2,3> tempMatr_1;
							tempMatr_1 << tj, gridNumb(0,5), ti, 
								tm, gridNumb(0,3), ti;
							poinNumb(0).insert(map<COORDINATE,Matrix<long,2,3>>::value_type(
								vem(0).poinCoor[tempMatr_0(0)], tempMatr_1));
							tempMatr_1 << tj + 1, gridNumb(0,5), ti, 
								tm, gridNumb(0,3), ti;
							poinNumb(0).insert(map<COORDINATE,Matrix<long,2,3>>::value_type(
								vem(0).poinCoor[tempMatr_0(1)], tempMatr_1));
							tempMatr_1 << tj + 1, gridNumb(0,5), ti, 
								tm + 1, gridNumb(0,3), ti;
							poinNumb(0).insert(map<COORDINATE,Matrix<long,2,3>>::value_type(
								vem(0).poinCoor[tempMatr_0(2)], tempMatr_1));
							tempMatr_1 << tj, gridNumb(0,5), ti, 
								tm + 1, gridNumb(0,3), ti;
							poinNumb(0).insert(map<COORDINATE,Matrix<long,2,3>>::value_type(
								vem(0).poinCoor[tempMatr_0(3)], tempMatr_1));
						}
					}
				}
			}
		}
	}
	return 1;
}

long DEHW::WHEE_GRID(){
	OUTPUT_TIME("Wheel grid");
	//profile
	wheeProf.resize(gridNumb(1,4) + 1, gridNumb(1,3) + 1);
	Matrix<Matrix<double,2,2>,Dynamic,Dynamic> wheeCurv;
	wheeCurv.resize(gridNumb(1,4) + 1, gridNumb(1,3) + 1);
	Matrix<double,3,3> tempMatr_0;
	tempMatr_0 << cos(backAngl(1)), -sin(backAngl(1)), 0.0,
		-sin(backAngl(1)), -cos(backAngl(1)), 0.0,
		0.0, 0.0, -1.0;
	for(long ti = 0; ti <= gridNumb(1,4); ti ++){
		for(long tj = 0; tj <= gridNumb(1,3); tj ++){
			long id_ti = (gridNumb(1,4) - ti) * pow2[refiLeve];
			long id_tj = (gridNumb(1,3) - tj) * pow2[refiLeve];
			wheeProf(ti,tj).block(0,0,3,1) = cartCoor(1)(id_ti, id_tj);
			wheeProf(gridNumb(1,4) - ti,tj).block(0,1,3,1) 
				= tempMatr_0 * cartCoor(1)(id_ti, id_tj);
			wheeCurv(ti,tj).block(0,0,2,1) = curvCoor(1)(id_ti, id_tj);
			wheeCurv(ti,tj).block(0,1,2,1) = curvCoor(1)(id_ti, id_tj);
		}
	}
	//point
	double tempAngl_1 = analAngl(1) - 2.0 * PI / z(1) * 2.0;
	Matrix3d tempRota;
	tempRota << cos(tempAngl_1), -sin(tempAngl_1), 0.0,
		sin(tempAngl_1), cos(tempAngl_1), 0.0,
		0.0, 0.0, 1.0;
	Matrix<Matrix<Matrix<long,Dynamic,Dynamic>,Dynamic,Dynamic>,Dynamic,Dynamic> blocPoin;
	blocPoin.resize(gridNumb(1,5),1);
	blocPoin(0).resize(gridNumb(1,4) + 1, 4);
	for(long ti = 0; ti <= gridNumb(1,4); ti ++){
		blocPoin(0)(ti,0).resize(gridNumb(1,1) + 1, gridNumb(1,0) + 1);
		blocPoin(0)(ti,1).resize(gridNumb(1,2) + 1, gridNumb(1,0) / 2 + 1);
		blocPoin(0)(ti,2).resize(gridNumb(1,3) + 1, 2 * gridNumb(1,2) + 1);
		blocPoin(0)(ti,3).resize(gridNumb(1,2) + 1, gridNumb(1,0) / 2 + 1);
		//transition into "unfolded cone surface"
		double tempAlph_3 = wheeCurv(ti,0)(0,0);
		Matrix<Matrix<double,2,1>,Dynamic,Dynamic> wheeProf_;
		wheeProf_.resize(gridNumb(1,3) + 1, 2);
		for(long tk = 0; tk <= 1; tk ++){
			for(long tj = 0; tj <= gridNumb(1,3); tj ++){
				double r_2 = sqrt(pow(wheeProf(ti,tj)(0,tk), 2.0) 
					+ pow(wheeProf(ti,tj)(1,tk), 2.0));
				double alph_2 = atan2(wheeProf(ti,tj)(1,tk), wheeProf(ti,tj)(0,tk));
				double r_3 = wheeCurv(ti,tj)(1,tk);
				double r_1 = a_h2 / cos(tempAlph_3) - r_3;
				double alph_1 = r_2 * alph_2 / r_1;
				wheeProf_(tj,tk) << r_1 * cos(alph_1), r_1 * sin(alph_1);
			}
		}
		Matrix<double,Dynamic,Dynamic> rootProf_0, rootProf_1;
		double r_1f = a_h2 / cos(tempAlph_3) - (a_h2 - d_f(1) / 2.0);
		WHEE_ROOT(rootProf_0, wheeProf_, 0, r_1f, pitcAngl * cos(tempAlph_3));
		WHEE_ROOT(rootProf_1, wheeProf_, 1, r_1f, pitcAngl * cos(tempAlph_3));
		double tranRadi_0 = r_1f - PI / 4.0 * m_t;
		Matrix<double,2,1> tempAngl_0;
		tempAngl_0(0) = atan2(rootProf_0(1,0), rootProf_0(0,0));
		tempAngl_0(1) = atan2(rootProf_1(1,0), rootProf_1(0,0));
		Matrix<double,Dynamic,Dynamic> tranProf_0, tranProf_1, inneProf;
		tranProf_0.resize(2, gridNumb(1,0) + 1);
		tranProf_1.resize(3, gridNumb(1,0) + 1);
		inneProf.resize(3, gridNumb(1,0) + 1);
		for(long tj = 0; tj <= gridNumb(1,0); tj ++){
			double tempAngl_j = tempAngl_0(0) 
				+ (tempAngl_0(1) - tempAngl_0(0)) / gridNumb(1,0) * tj;
			tranProf_0.block(0,tj,2,1) << tranRadi_0 * cos(tempAngl_j), 
				tranRadi_0 * sin(tempAngl_j);
			tranProf_1.block(0,tj,3,1) = WHEE_CONE(tranProf_0.block(0,tj,2,1), tempAlph_3);
			double tempRadi = sqrt(pow(tranProf_1(0,tj), 2.0) + pow(tranProf_1(1,tj), 2.0));
			inneProf.block(0,tj,3,1) << inneRadi(1) / tempRadi * tranProf_1(0,tj), 
				inneRadi(1) / tempRadi * tranProf_1(1,tj), 
				tranProf_1(2,tj);
		}
		//0
		COORDINATE tempCoor;
		for(long tj = 0; tj <= gridNumb(1,1); tj ++){
			for(long tk = 0; tk <= gridNumb(1,0); tk ++){
				Vector3d tempVect = 
					(1.0 - (double)tj / gridNumb(1,1)) * inneProf.block(0,tk,3,1).eval() 
					+ (double)tj / gridNumb(1,1) * tranProf_1.block(0,tk,3,1).eval();
				tempCoor.COPY(tempRota * tempVect);
				blocPoin(0)(ti,0)(tj,tk) = vem(1).TRY_ADD_POINT(tempCoor);
			}
		}
		//1
		for(long tj = 0; tj <= gridNumb(1,2); tj ++){
			for(long tk = 0; tk <= gridNumb(1,0) / 2; tk ++){
				Matrix<double,2,1> tempPoin = 
					(1.0 - (double)tj / gridNumb(1,2)) * tranProf_0.block(0,tk,2,1).eval() 
					+ (double)tj / gridNumb(1,2) * rootProf_0.block(0,tk,2,1).eval();
				Vector3d tempVect = WHEE_CONE(tempPoin, tempAlph_3);
				tempCoor.COPY(tempRota * tempVect);
				blocPoin(0)(ti,1)(tj,tk) = vem(1).TRY_ADD_POINT(tempCoor);
			}
		}
		//2
		Matrix<double,Dynamic,Dynamic> lineCoor;
		lineCoor.resize(2, gridNumb(1,3) + 1);
		lineCoor.block(0,0,2,1) = tranProf_0.block(0,gridNumb(1,0) / 2,2,1);
		lineCoor.block(0,gridNumb(1,3),2,1) = 
			0.5 * (wheeProf_(gridNumb(1,3),0) + wheeProf_(gridNumb(1,3),1));
		for(long tj = 1; tj < gridNumb(1,3); tj ++){
			lineCoor.block(0,tj,2,1) = lineCoor.block(0,0,2,1) 
				+ (lineCoor.block(0,gridNumb(1,3),2,1) 
				- lineCoor.block(0,0,2,1)) / gridNumb(1,3) * tj;
		}
		for(long tj = 0; tj <= gridNumb(1,3); tj ++){
			for(long tk = 0; tk <= gridNumb(1,2); tk ++){
				Matrix<double,2,1> tempPoin = 
					(1.0 - (double)tk / gridNumb(1,2)) * wheeProf_(tj,0) 
					+ (double)tk / gridNumb(1,2) * lineCoor.block(0,tj,2,1).eval();
				Vector3d tempVect = WHEE_CONE(tempPoin, tempAlph_3);
				tempCoor.COPY(tempRota * tempVect);
				blocPoin(0)(ti,2)(tj,tk) = vem(1).TRY_ADD_POINT(tempCoor);
			}
			for(long tk = 1; tk <= gridNumb(1,2); tk ++){
				Matrix<double,2,1> tempPoin = 
					(1.0 - (double)tk / gridNumb(1,2)) * lineCoor.block(0,tj,2,1).eval() 
					+ (double)tk / gridNumb(1,2) * wheeProf_(tj,1);
				Vector3d tempVect = WHEE_CONE(tempPoin, tempAlph_3);
				tempCoor.COPY(tempRota * tempVect);
				blocPoin(0)(ti,2)(tj,gridNumb(1,2) + tk) = vem(1).TRY_ADD_POINT(tempCoor);
			}
		}
		//3
		for(long tj = 0; tj <= gridNumb(1,2); tj ++){
			for(long tk = 0; tk <= gridNumb(1,0) / 2; tk ++){
				Matrix<double,2,1> tempPoin = (1.0 - (double)tj / gridNumb(1,2)) 
					* tranProf_0.block(0,gridNumb(1,0) / 2 + tk,2,1).eval() 
					+ (double)tj / gridNumb(1,2) 
					* rootProf_1.block(0,gridNumb(1,0) / 2 - tk,2,1).eval();
				Vector3d tempVect = WHEE_CONE(tempPoin, tempAlph_3);
				tempCoor.COPY(tempRota * tempVect);
				blocPoin(0)(ti,3)(tj,tk) = vem(1).TRY_ADD_POINT(tempCoor);
			}
		}
	}
	for(long ti = 1; ti < gridNumb(1,5); ti ++){
		blocPoin(ti).resize(blocPoin(0).rows(), blocPoin(0).cols());
		double tempAngl = 2.0 * PI / z(1) * ti;
		Matrix<double,3,3> tempMatr;
		tempMatr << cos(tempAngl), -sin(tempAngl), 0.0,
			sin(tempAngl), cos(tempAngl), 0.0,
			0.0, 0.0, 1.0;
		for(long tj = 0; tj < blocPoin(0).rows(); tj ++){
			for(long tk = 0; tk < blocPoin(0).cols(); tk ++){
				blocPoin(ti)(tj,tk).resize(blocPoin(0)(tj,tk).rows(), blocPoin(0)(tj,tk).cols());
				for(long tm = 0; tm < blocPoin(0)(tj,tk).rows(); tm ++){
					for(long tn = 0; tn < blocPoin(0)(tj,tk).cols(); tn ++){
						Vector3d tempVect = vem(1).poinCoor[blocPoin(0)(tj,tk)(tm,tn)];
						COORDINATE tempCoor;
						tempCoor.COPY(tempMatr * tempVect);
						blocPoin(ti)(tj,tk)(tm,tn) = vem(1).TRY_ADD_POINT(tempCoor);
					}
				}
			}
		}
	}
	//volume
	for(long ti = 0; ti < gridNumb(1,5); ti ++){
		for(long tj = 0; tj < blocPoin(ti).rows() - 1; tj ++){
			if(tj % 10 == 0){
				cout<<"Volume layer "<< tj <<"/"<< blocPoin(ti).rows() - 1
					<<"~"<< ti <<"/"<< gridNumb(1,5) <<endl;
			}
			for(long tk = 0; tk < blocPoin(ti).cols(); tk ++){
				for(long tm = 0; tm < blocPoin(0)(tj,tk).rows() - 1; tm ++){
					for(long tn = 0; tn < blocPoin(0)(tj,tk).cols() - 1; tn ++){
						Matrix<long,8,1> tempCorn;
						tempCorn << blocPoin(ti)(tj,tk)(tm,tn + 1),
							blocPoin(ti)(tj,tk)(tm,tn),
							blocPoin(ti)(tj,tk)(tm + 1,tn),
							blocPoin(ti)(tj,tk)(tm + 1,tn + 1),
							blocPoin(ti)(tj + 1,tk)(tm,tn + 1),
							blocPoin(ti)(tj + 1,tk)(tm,tn),
							blocPoin(ti)(tj + 1,tk)(tm + 1,tn),
							blocPoin(ti)(tj + 1,tk)(tm + 1,tn + 1);
						long tempIden = vem(1).ADD_BLOCK_FROM_POINT(tempCorn);
						if(tk == 2 && tn == 0){
							MATRIX tempMatr_2;
							tempMatr_2.resize(4,1);
							tempMatr_2 << blocPoin(ti)(tj,tk)(tm,tn), 
								blocPoin(ti)(tj + 1,tk)(tm,tn), 
								blocPoin(ti)(tj + 1,tk)(tm + 1,tn),
								blocPoin(ti)(tj,tk)(tm + 1,tn);
							segmVolu(1).push_back(make_pair(tempMatr_2, -tempIden));
							Matrix<long,2,3> tempMatr_3;
							tempMatr_3 << tj, gridNumb(1,4), ti, 
								tm, gridNumb(1,3), ti;
							poinNumb(1).insert(map<COORDINATE,Matrix<long,2,3>>::value_type(
								vem(1).poinCoor[tempMatr_2(0)], tempMatr_3));
							tempMatr_3 << tj + 1, gridNumb(1,4), ti,  
								tm, gridNumb(1,3), ti;
							poinNumb(1).insert(map<COORDINATE,Matrix<long,2,3>>::value_type(
								vem(1).poinCoor[tempMatr_2(1)], tempMatr_3));
							tempMatr_3 << tj + 1, gridNumb(1,4), ti, 
								tm + 1, gridNumb(1,3), ti;
							poinNumb(1).insert(map<COORDINATE,Matrix<long,2,3>>::value_type(
								vem(1).poinCoor[tempMatr_2(2)], tempMatr_3));
							tempMatr_3 << tj, gridNumb(1,4), ti, 
								tm + 1, gridNumb(1,3), ti;
							poinNumb(1).insert(map<COORDINATE,Matrix<long,2,3>>::value_type(
								vem(1).poinCoor[tempMatr_2(3)], tempMatr_3));
						}
					}
				}
			}
		}
	}
	return 1;
}

long DEHW::WORM_ROOT_RADIUS(long flag, Matrix<double,2,3> tempPoin, 
	Matrix<double,2,1> &tempCent, double &tempRadi, Matrix<double,2,1> &tempAngl){
	Vector2d vect_1 = tempPoin.block(0,1,2,1) - tempPoin.block(0,0,2,1);
	vect_1 = vect_1.eval() / vect_1.norm();
	Vector2d vect_2 = tempPoin.block(0,2,2,1) - tempPoin.block(0,0,2,1);
	double tempLeng_1 = vect_2.dot(vect_1);
	double tempLeng_2 = sqrt(vect_2(0) * vect_2(0) + vect_2(1) * vect_2(1) 
		- tempLeng_1 * tempLeng_1);
	double targVari = tempLeng_1 / (R_f(0) - tempLeng_2);
	double middAngl = asin(targVari / sqrt(1.0 + targVari * targVari)) - atan(1.0 / targVari);
	tempRadi = R_f(0) - tempLeng_1 / cos(middAngl);
	Vector2d tempVect;
	tempVect << flag * vect_1(1), - flag * vect_1(0);
	tempCent = tempPoin.block(0,0,2,1) + tempRadi * tempVect;
	tempAngl(0) = atan2(- tempVect(1), - tempVect(0));
	tempAngl(1) = tempAngl(0) + flag * (PI / 2.0 - middAngl);
	return 1;
}

long DEHW::WORM_ROOT(long id, long flag, Matrix<double,Dynamic,Dynamic> &rootProf){
	rootProf.resize(3, gridNumb(0,0) / 2 + 1);
	//
	Matrix<double,3,3> tempPoin;
	if(flag == 1){
		tempPoin.block(0,0,3,1) = wormProf(id,0).block(0,0,3,1).eval();
		tempPoin.block(0,1,3,1) = wormProf(id,1).block(0,0,3,1).eval();
	}
	else{
		tempPoin.block(0,0,3,1) = wormProf(id,0).block(0,1,3,1).eval();
		tempPoin.block(0,1,3,1) = wormProf(id,1).block(0,1,3,1).eval();
	}
	tempPoin.block(0,2,3,1) << a_h2 * cos(profCurv(id)), a_h2 * sin(profCurv(id)), 0.0;
	Matrix<double,2,3> tempPoin_;
	Matrix<double,3,3> tempMatr;
	tempMatr << 0.0, 0.0, 1.0,
		-cos(profCurv(id)), -sin(profCurv(id)), 0.0,
		sin(profCurv(id)), -cos(profCurv(id)), 0.0;
	for(long ti = 0; ti < 3; ti ++){
		Matrix<double,3,1> tempPoin_i = tempMatr * tempPoin.block(0,ti,3,1).eval();
		tempPoin_(0,ti) = tempPoin_i(0);
		tempPoin_(1,ti) = tempPoin_i(1) + a_h2;
	}
	Matrix<double,2,1> tempCent;
	double tempRadi;
	Matrix<double,2,1> tempAngl;
	WORM_ROOT_RADIUS(flag, tempPoin_, tempCent, tempRadi, tempAngl);
	Matrix<double,2,1> tempPoinArce;
	tempPoinArce << tempCent(0) + tempRadi * cos(tempAngl(1)),
		tempCent(1) + tempRadi * sin(tempAngl(1));
	double tempAnglArce = atan2(tempPoinArce(1), tempPoinArce(0));
	//
	double tempAnglRoot;
	double tempAnglStar = acos(r_b2 / (d(1) / 2.0))
		- i_2h * profCurv(id) - tootThicAngl(0) / 2.0;
	if(flag == 1){
		tempAnglRoot = tempAnglStar + pitcAngl / 2.0;
	}
	else{
		tempAnglRoot = tempAnglStar - pitcAngl / 2.0;
	}
	//
	double sumLeng = flag * R_f(0) * (tempAnglRoot - tempAnglArce)
		+ flag * tempRadi * (tempAngl(1) - tempAngl(0));
	Matrix<double,3,1> tempTran;
	tempTran << a_h2 * cos(profCurv(id)), a_h2 * sin(profCurv(id)), 0.0;
	for(long ti = 0; ti <= gridNumb(0,0) / 2; ti ++){
		double leng_i = sumLeng / (gridNumb(0,0) / 2.0) * (double)ti;
		Matrix<double,3,1> poin_i;
		if(leng_i <= flag * R_f(0) * (tempAnglRoot - tempAnglArce)){
			double angl_i = tempAnglRoot - flag * leng_i / R_f(0);
			poin_i << R_f(0) * cos(angl_i), R_f(0) * sin(angl_i), 0.0;
		}
		else{
			leng_i = leng_i - flag * R_f(0) * (tempAnglRoot - tempAnglArce);
			double angl_i = tempAngl(1) - flag * leng_i / tempRadi;
			poin_i << tempCent(0) + tempRadi * cos(angl_i), 
				tempCent(1) + tempRadi * sin(angl_i), 0.0;
		}
		rootProf.block(0,ti,3,1) = tempMatr.transpose() * poin_i.eval() + tempTran;
	}
	return 1;
}

Matrix<double,3,1> DEHW::WHEE_CONE(Matrix<double,2,1> tempXY, double tempAlph_3){
	Matrix<double,3,1> tempXYZ;
	double r_1 = sqrt(pow(tempXY(0), 2.0) + pow(tempXY(1), 2.0));
	double alph_1 = atan2(tempXY(1), tempXY(0));
	double r_2 = r_1 * cos(tempAlph_3);
	double alph_2 = r_1 * alph_1 / r_2;
	double r_3 = a_h2 / cos(tempAlph_3) - r_1;
	tempXYZ << r_2 * cos(alph_2), r_2 * sin(alph_2), r_3 * sin(tempAlph_3);
	return tempXYZ;
}

long DEHW::WHEE_ROOT(Matrix<double,Dynamic,Dynamic> &rootProf, 
	Matrix<Matrix<double,2,1>,Dynamic,Dynamic> tempProf, long id, double r_f, double tempPitc){
	rootProf.resize(2, gridNumb(1,0) / 2 + 1);
	//
	Matrix<double,2,3> tempPoin;
	tempPoin.block(0,0,2,1) = tempProf(0,id);
	tempPoin.block(0,1,2,1) = tempProf(1,id);
	tempPoin.block(0,2,2,1) << 0.0, 0.0;
	Vector2d vect_1 = tempPoin.block(0,0,2,1) - tempPoin.block(0,1,2,1);
	vect_1 = vect_1.eval() / vect_1.norm();
	Vector2d vect_2 = tempPoin.block(0,2,2,1) - tempPoin.block(0,0,2,1);
	double tempLeng_1 = vect_2.dot(vect_1);
	double tempLeng_2 = sqrt(vect_2(0) * vect_2(0) + vect_2(1) * vect_2(1) 
		- tempLeng_1 * tempLeng_1);
	double targVari = tempLeng_1 / (r_f - tempLeng_2);
	double middAngl = asin(targVari / sqrt(1.0 + targVari * targVari)) - atan(1.0 / targVari);
	double tempRadi = tempLeng_1 / cos(middAngl) - r_f;
	double tempSign = (id == 0) ? 1.0 : -1.0;
	Vector2d tempVect;
	tempVect << -tempSign * vect_1(1), tempSign * vect_1(0);
	Matrix<double,2,1> tempCent = tempPoin.block(0,0,2,1) + tempRadi * tempVect;
	Matrix<double,2,1> tempAngl;
	tempAngl(0) = atan2(-tempVect(1), -tempVect(0));
	tempAngl(1) = tempAngl(0) + tempSign * (PI / 2.0 - middAngl);
	Matrix<double,2,1> tempPoinArce;
	tempPoinArce << tempCent(0) + tempRadi * cos(tempAngl(1)),
		tempCent(1) + tempRadi * sin(tempAngl(1));
	double tempAnglArce = atan2(tempPoinArce(1), tempPoinArce(0));
	//
	double tempAnglRoot;
	Matrix<double,2,1> tempPoinOppo = tempProf(0, 1 - id);
	tempAnglRoot = (atan2(tempPoinOppo(1), tempPoinOppo(0)) 
		+ atan2(tempPoin(1,0), tempPoin(0,0))) / 2.0;
	tempAnglRoot -= tempSign * tempPitc / 2.0;
	//
	double sumLeng = r_f * tempSign * (tempAnglArce - tempAnglRoot)
		+ tempRadi * tempSign * (tempAngl(1) - tempAngl(0));
	for(long ti = 0; ti <= gridNumb(1,0) / 2; ti ++){
		double leng_i = sumLeng / (gridNumb(1,0) / 2) * ti;
		if(leng_i <= r_f * tempSign * (tempAnglArce - tempAnglRoot)){
			double angl_i = tempAnglRoot + tempSign * leng_i / r_f;
			rootProf.block(0,ti,2,1) << r_f * cos(angl_i), r_f * sin(angl_i);
		}
		else{
			leng_i = leng_i - r_f * tempSign * (tempAnglArce - tempAnglRoot);
			double angl_i = tempAngl(1) - tempSign * leng_i / tempRadi;
			rootProf.block(0,ti,2,1) << tempCent(0) + tempRadi * cos(angl_i), 
				tempCent(1) + tempRadi * sin(angl_i);
		}
	}
	return 1;
}

bool DEHW::ADAPTIVE_REFINE(double tempCrit, long buckFact){
	OUTPUT_TIME("Contact search");
	BUCKET_SORT_R(buckFact);
	for(long tv = 0; tv < 2; tv ++){
		for(long ti = 0; ti < segmVolu(tv).size(); ti ++){
			MATRIX tempMatr = (segmVolu(tv)[ti]).first;
			for(long tj = 0; tj < tempMatr.rows(); tj ++){
				COORDINATE tempCoor = vem(tv).poinCoor[tempMatr(tj)];
				double tempDist = DISTANCE(tempCoor, 1 - tv);
				if(abs(tempDist) < (tv + 1) * tempCrit){
					(segmVolu(tv)[ti]).second = - (segmVolu(tv)[ti]).second;
					break;
				}
			}
		}
	}
	OUTPUT_TIME("Refine");
	Matrix<map<COORDINATE,Matrix<long,2,3>>,2,1> poinNumb_0;
	Matrix<list<long>,2,1> spliVolu;
	for(long tv = 0; tv < 2; tv ++){
		Matrix<long,Dynamic,Dynamic> spliFlag;
		map<MATRIX,COORDINATE> planSurf;
		for(long ti = 0; ti < segmVolu(tv).size(); ti ++){
			long tempIden = (segmVolu(tv)[ti]).second;
			if(tempIden > 0){
				spliVolu(tv).push_back(tempIden);
				MATRIX tempMatr_0 = (segmVolu(tv)[ti]).first;
				MATRIX tempMatr_1;
				tempMatr_1.resize(2,1);
				tempMatr_1 << tempMatr_0(0), tempMatr_0(1);
				COORDINATE tempCoor = PLNAE_2_SURFACE(tv, tempMatr_1, poinNumb_0);
				planSurf.insert(map<MATRIX,COORDINATE>::value_type(tempMatr_1, tempCoor));
				tempMatr_1.resize(2,1);
				tempMatr_1 << tempMatr_0(1), tempMatr_0(2);
				tempCoor = PLNAE_2_SURFACE(tv, tempMatr_1, poinNumb_0);
				planSurf.insert(map<MATRIX,COORDINATE>::value_type(tempMatr_1, tempCoor));
				tempMatr_1.resize(2,1);
				tempMatr_1 << tempMatr_0(2), tempMatr_0(3);
				tempCoor = PLNAE_2_SURFACE(tv, tempMatr_1, poinNumb_0);
				planSurf.insert(map<MATRIX,COORDINATE>::value_type(tempMatr_1, tempCoor));
				tempMatr_1.resize(2,1);
				tempMatr_1 << tempMatr_0(3), tempMatr_0(0);
				tempCoor = PLNAE_2_SURFACE(tv, tempMatr_1, poinNumb_0);
				planSurf.insert(map<MATRIX,COORDINATE>::value_type(tempMatr_1, tempCoor));
				tempMatr_1.resize(4,1);
				tempMatr_1 << tempMatr_0(0), tempMatr_0(1), tempMatr_0(2), tempMatr_0(3);
				tempCoor = PLNAE_2_SURFACE(tv, tempMatr_1, poinNumb_0);
				planSurf.insert(map<MATRIX,COORDINATE>::value_type(tempMatr_1, tempCoor));
			}
		}
		spliFlag.resize(spliVolu(tv).size(),1);
		for(long ti = 0; ti < spliVolu(tv).size(); ti ++){
			spliFlag(ti) = 4;
		}
		vem(tv).REFINE(spliVolu(tv), spliFlag, planSurf);
	}
	OUTPUT_TIME("Update data structure");
	poinNumb = poinNumb_0;
	Matrix<vector<pair<MATRIX,long>>,2,1> segmVolu_0;
	for(long tv = 0; tv < 2; tv ++){
		for(list<long>::const_iterator iter_0 = spliVolu(tv).begin(); 
			iter_0 != spliVolu(tv).end(); iter_0 ++){
			long tempVolu = *iter_0;
			map<long,MATRIX>::iterator iter_V = vem(tv).voluFace.find(tempVolu);
			for(long ti = 0; ti < (iter_V->second).rows(); ti ++){
				long tempFace = (iter_V->second)(ti);
				map<long,MATRIX>::iterator iter_F = vem(tv).faceLine.find(abs(tempFace));
				MATRIX refiLine = iter_F->second;
				refiLine.FLAG_SORT(tempFace);
				bool faceFlag = true;
				MATRIX faceMatr;
				faceMatr.resize(refiLine.rows(),1);
				for(long tj = 0; tj < refiLine.rows(); tj ++){
					map<long,MATRIX>::iterator iter_L = vem(tv).linePoin.find(abs(refiLine(tj)));
					if(refiLine(tj) > 0){
						faceMatr(tj) = (iter_L->second)(0);
					}
					else{
						faceMatr(tj) = (iter_L->second)(1);
					}
					map<long,COORDINATE>::iterator face_P = vem(tv).poinCoor.find(faceMatr(tj));
					map<COORDINATE,Matrix<long,2,3>>::iterator iter_1 = 
						poinNumb(tv).find(face_P->second);
					if(iter_1 == poinNumb(tv).end()){
						faceFlag = false;
						break;
					}
				}
				if(faceFlag){
					segmVolu_0(tv).push_back(make_pair(faceMatr, -tempVolu));
					break;
				}
			}
		}
	}
	segmVolu = segmVolu_0;
	return true;
}

COORDINATE DEHW::PLNAE_2_SURFACE(long tv, MATRIX tempMatr,
	Matrix<map<COORDINATE,Matrix<long,2,3>>,2,1> &poinNumb_0){
	Matrix<long,2,1> centNumb;
	centNumb << 0, 0;
	long tempToot;
	Matrix<long,2,3> tempIden;
	tempIden << 0, 0, 0,
		0, 0, 0;
	for(long ti = 0; ti < tempMatr.rows(); ti ++){
		COORDINATE tempCoor = vem(tv).poinCoor[tempMatr(ti)];
		map<COORDINATE,Matrix<long,2,3>>::iterator iter_0 = poinNumb(tv).find(tempCoor);
		Matrix<long,2,3> tempNumb = iter_0->second;
		centNumb(0) = centNumb(0) + tempNumb(0,0) * ((curvCoor(tv).rows() - 1) / tempNumb(0,1));
		centNumb(1) = centNumb(1) + tempNumb(1,0) * ((curvCoor(tv).cols() - 1) / tempNumb(1,1));
		tempIden(0,0) = tempIden(0,0) + 2 * tempNumb(0,0);
		tempIden(1,0) = tempIden(1,0) + 2 * tempNumb(1,0);
		if(ti == 0){
			tempToot = tempNumb(0,2);
			tempIden(0,1) = 2 * tempNumb(0,1);
			tempIden(1,1) = 2 * tempNumb(1,1);
			tempIden(0,2) = tempNumb(0,2);
			tempIden(1,2) = tempNumb(1,2);
		}
		tempNumb(0,0) = 2 * tempNumb(0,0);
		tempNumb(0,1) = 2 * tempNumb(0,1);
		tempNumb(1,0) = 2 * tempNumb(1,0);
		tempNumb(1,1) = 2 * tempNumb(1,1);
		poinNumb_0(tv).insert(map<COORDINATE,Matrix<long,2,3>>::value_type(tempCoor, tempNumb));
	}
	tempIden(0,0) = tempIden(0,0) /tempMatr.rows();
	tempIden(1,0) = tempIden(1,0) /tempMatr.rows();
	centNumb = centNumb / tempMatr.rows();
	Matrix<double,3,1> tempCoor_0;
	if(tv == 1){
		tempCoor_0 = cartCoor(tv)(curvCoor(tv).rows() - 1 - centNumb(0), 
			curvCoor(tv).cols() - 1 - centNumb(1));
		double tempAngl = analAngl(tv) - 2.0 * PI / z(1) * 2.0 + 2.0 * PI / z(1) * tempToot;
		Matrix3d tempRota;
		tempRota << cos(tempAngl), -sin(tempAngl), 0.0,
			sin(tempAngl), cos(tempAngl), 0.0,
			0.0, 0.0, 1.0;
		tempCoor_0 = tempRota * tempCoor_0;
	}
	else{
		tempCoor_0 = cartCoor(tv)(centNumb(0), centNumb(1));
		double tempAngl = analAngl(tv) + 2.0 * PI / z(0) * tempToot;
		Matrix3d tempRota;
		tempRota << cos(tempAngl), -sin(tempAngl), 0.0, 
			0.0, 0.0, 1.0, 
			-sin(tempAngl), -cos(tempAngl), 0.0;
		Vector3d tempTran;
		tempTran << -(a_h2 + centErro), 0.0, 0.0;
		tempCoor_0 = tempRota * tempCoor_0 + tempTran;
	}
	COORDINATE tempCoor_1;
	tempCoor_1.COPY(tempCoor_0);
	poinNumb_0(tv).insert(map<COORDINATE,Matrix<long,2,3>>::value_type(tempCoor_1, tempIden));
	return tempCoor_1;
}

bool DEHW::BUCKET_SORT_R(long buckFact){
	for(long tv = 0; tv < 2; tv ++){
		Matrix<double,Dynamic,Dynamic> locaCoor = MatrixXd::Zero(segmVolu(tv).size(),2);
		for(long ti = 0; ti < segmVolu(tv).size(); ti ++){
			MATRIX tempMatr = (segmVolu(tv)[ti]).first;
			COORDINATE tempCoor_0;
			tempCoor_0.COPY(MatrixXd::Zero(3,1));
			for(long tj = 0; tj < tempMatr.rows(); tj ++){
				COORDINATE tempCoor_1 = vem(tv).poinCoor[tempMatr(tj)];
				tempCoor_0 = tempCoor_0 + tempCoor_1;
			}
			tempCoor_0 = tempCoor_0 / tempMatr.rows();
			locaCoor(ti,0) = tempCoor_0(1);
			locaCoor(ti,1) = sqrt(pow(tempCoor_0(0) + (a_h2+centErro), 2.0)
				+ pow(tempCoor_0(2), 2.0));
		}
		buckNumb_R.block(tv,0,1,2) << gridNumb(0,4) /2 * buckFact, 
			gridNumb(0,3) / 2 * buckFact;
		double tempMini = locaCoor.block(0,0,locaCoor.rows(),1).minCoeff();
		double tempMaxi = locaCoor.block(0,0,locaCoor.rows(),1).maxCoeff();
		buckLoca_R(tv)(0,0) = tempMini - (tempMaxi-tempMini) / buckNumb_R(tv,0);
		buckLoca_R(tv)(0,1) = tempMaxi + (tempMaxi-tempMini) / buckNumb_R(tv,0);
		tempMini = locaCoor.block(0,1,locaCoor.rows(),1).minCoeff();
		tempMaxi = locaCoor.block(0,1,locaCoor.rows(),1).maxCoeff();
		buckLoca_R(tv)(1,0) = tempMini - (tempMaxi - tempMini) / buckNumb_R(tv,1);
		buckLoca_R(tv)(1,1) = tempMaxi + (tempMaxi - tempMini) / buckNumb_R(tv,1);
		buckLoca_R(tv)(0,2) = (buckLoca_R(tv)(0,1) - buckLoca_R(tv)(0,0)) / buckNumb_R(tv,0);
		buckLoca_R(tv)(1,2) = (buckLoca_R(tv)(1,1) - buckLoca_R(tv)(1,0)) / buckNumb_R(tv,1);
		buck_R(tv).resize(buckNumb_R(tv,0), buckNumb_R(tv,1));
		for(long ti = 0; ti < segmVolu(tv).size(); ti ++){
			long buck_r = (locaCoor(ti,0) - buckLoca_R(tv)(0,0)) / buckLoca_R(tv)(0,2);
			long buck_c = (locaCoor(ti,1) - buckLoca_R(tv)(1,0)) / buckLoca_R(tv)(1,2);
			buck_R(tv)(buck_r, buck_c).push_back(ti);
		}
	}
	return true;
}

double DEHW::DISTANCE(COORDINATE tempCoor, long voluIden){
	double infiDist = 1.0E12;
	Matrix<double,2,1> tempLoca;
	tempLoca << tempCoor(1), 
		sqrt(pow(tempCoor(0) + (a_h2 + centErro),2.0) + pow(tempCoor(2),2.0));
	long buck_r = (tempLoca(0) - buckLoca_R(voluIden)(0,0)) / buckLoca_R(voluIden)(0,2);
	long buck_c = (tempLoca(1) - buckLoca_R(voluIden)(1,0)) / buckLoca_R(voluIden)(1,2);
	if(buck_r < 0 || buck_r > buckNumb_R(voluIden,0) - 1 
		|| buck_c < 0 || buck_c > buckNumb_R(voluIden,1) - 1){
		return infiDist;
	}
	for(long ti = max(buck_r - 1, (long)0); 
		ti <= min(buck_r + 1, buckNumb_R(voluIden,0) - 1); ti ++){
		for(long tj = max(buck_c - 1, (long)0); 
			tj <= min(buck_c + 1, buckNumb_R(voluIden,1) - 1); tj ++){
			for(list<long>::iterator iter_0 = buck_R(voluIden)(ti, tj).begin();
				iter_0 != buck_R(voluIden)(ti, tj).end(); iter_0 ++){
				MATRIX tempMatr = (segmVolu(voluIden)[*iter_0]).first;
				COORDINATE tempVect_0, tempVect_1, tempVect_2, tempVect_3;
				tempVect_0 = vem(voluIden).poinCoor[tempMatr(0)];
				tempVect_1 = vem(voluIden).poinCoor[tempMatr(1)];
				tempVect_2 = vem(voluIden).poinCoor[tempMatr(2)];
				tempVect_3 = vem(voluIden).poinCoor[tempMatr(3)];
				double tempDist = 
					TRIANGLE_DISTANCE(tempCoor, tempVect_0, tempVect_1, tempVect_3);
				if(tempDist < infiDist){
					return tempDist;
				}
				tempDist = TRIANGLE_DISTANCE(tempCoor, tempVect_1, tempVect_2, tempVect_3);
				if(tempDist < infiDist){
					return tempDist;
				}
			}
		}
	}
	return infiDist;
}

double DEHW::TRIANGLE_DISTANCE(COORDINATE tempCoor, COORDINATE tempCoor_0, 
	COORDINATE tempCoor_1, COORDINATE tempCoor_2){
	double infiDist = 1.0E12;
	Vector3d mortPoin;
	tempCoor.COPY_TO(mortPoin);
	Matrix<Vector3d,3,1> nonmTria;
	tempCoor_0.COPY_TO(nonmTria(0));
	tempCoor_1.COPY_TO(nonmTria(1));
	tempCoor_2.COPY_TO(nonmTria(2));
	//projection
	Matrix<Vector3d,2,1> tangVect;
	tangVect(0) = nonmTria(1) - nonmTria(0);
	tangVect(1) = nonmTria(2) - nonmTria(0);
	Vector3d normVect = tangVect(0).cross(tangVect(1));
	normVect.normalize();
	Matrix<Vector3d,2,1> unitVect;
	unitVect(0) = tangVect(0).normalized();
	unitVect(1) = normVect.cross(unitVect(0));
	Vector2d mortProj;
	Matrix<Vector2d,3,1> nonmProj;
	nonmProj(0) << 0.0, 0.0;
	nonmProj(1) << tangVect(0).dot(unitVect(0)), tangVect(0).dot(unitVect(1));
	nonmProj(2) << tangVect(1).dot(unitVect(0)), tangVect(1).dot(unitVect(1));
	Vector3d tempVect = mortPoin - nonmTria(0);
	mortProj << tempVect.dot(unitVect(0)), tempVect.dot(unitVect(1));
	if(IN_TRIANGLE_2D(mortProj, nonmProj)){	
		double triaArea = TRIANGLE_AREA_2D(nonmProj(0), nonmProj(1), nonmProj(2));
		Matrix<double,3,1> tempShap, tempPoin;
		for(long tj = 0; tj < 3; tj ++){
			tempShap((tj + 2) % 3) = 
				TRIANGLE_AREA_2D(nonmProj(tj), nonmProj((tj + 1) % 3), mortProj) / triaArea;
		}
		tempPoin = MatrixXd::Zero(3,1);
		for(long tj = 0; tj < 3; tj ++){
			tempPoin = tempPoin	+ tempShap(tj) * nonmTria(tj);
		}
		Vector3d nonmMort = mortPoin - tempPoin;
		return nonmMort.dot(normVect);
	}
	return infiDist;
}

long DEHW::SOLVE(){
	//
	refiLeve = 4;//must be greater than 1
	TS_GRID();
	//
	analNumb = 40;
	for(long ta = 0; ta < analNumb; ta ++){
		//
		stringstream tempStre;
		tempStre << "Dehw/Resu_" << ta;
		mkdir(tempStre.str().c_str(), 0755);
		tempStre << "/";
		tempStre >> outpDire;
		tempStre.str("");
		tempStre.clear();
		//
		analAngl << 2.0 * PI / z(0) / analNumb * ta, 2.0 * PI / z(1) / analNumb * ta;
		WORM_GRID();
		WHEE_GRID();
		for(long ti = 0; ti < refiLeve; ti ++){
			ADAPTIVE_REFINE(distCrit(ti), pow2[ti]);
		}
		Matrix<set<long>,2,1> coupFree;
		for(long tv = 0; tv < 2; tv ++){
			vem(tv).VOLUME_2_ELEMENT();
			vem(tv).OUTPUT_ELEMENT(tv);
			//
			vem(tv).freeTran.resize(vem(tv).nodeCoor.rows(),1);
			vector<long> consFree;
			for(long ti = 0; ti < vem(tv).nodeCoor.rows(); ti ++){
				double tempX = vem(tv).nodeCoor(ti,0);
				double tempY = vem(tv).nodeCoor(ti,1);
				double tempZ = vem(tv).nodeCoor(ti,2);
				if(tv == 1){
					double tempR = sqrt(tempX * tempX + tempY * tempY);
					if(abs(tempR - inneRadi(1)) < 1.0E-10){
						double tempAngl = atan2(tempY, tempX);
						vem(tv).freeTran(ti) << cos(tempAngl), -sin(tempAngl), 0.0,
							sin(tempAngl), cos(tempAngl), 0.0,
							0.0, 0.0, 1.0;
						consFree.push_back(3 * ti + 0);
						consFree.push_back(3 * ti + 1);
						consFree.push_back(3 * ti + 2);
					}
					else{
						vem(tv).freeTran(ti) << 1.0, 0.0, 0.0,
							0.0, 1.0, 0.0,
							0.0, 0.0, 1.0;
					}
				}
				else{
					tempX += (a_h2 + centErro);
					double tempR = sqrt(tempX * tempX + tempZ * tempZ);
					if(abs(tempR - inneRadi(0)) < 1.0E-10){
						double tempAngl = atan2(tempX, tempZ);
						vem(tv).freeTran(ti) << sin(tempAngl), cos(tempAngl), 0.0,
							0.0, 0.0, 1.0,
							cos(tempAngl), -sin(tempAngl), 0.0;
						consFree.push_back(3 * ti + 0);
						consFree.push_back(3 * ti + 2);
						coupFree(tv).insert(3 * ti + 1);
					}
					else{
						vem(tv).freeTran(ti) << 1.0, 0.0, 0.0,
							0.0, 1.0, 0.0,
							0.0, 0.0, 1.0;
					}
				}
			}
			vem(tv).CONSTRAINT(consFree);
			vem(tv).COUPLE(coupFree(tv));
		}
		//
		OUTPUT_TIME("Load");
		Matrix<vector<pair<long,double>>,2,1> loadList;
		double tempLoad = inpuTorq / inneRadi(0) / coupFree(0).size();
		for(set<long>::iterator iter_0 = coupFree(0).begin(); 
			iter_0 != coupFree(0).end(); iter_0 ++){
			loadList(0).push_back(make_pair(*iter_0, tempLoad));
		}
		//
		OUTPUT_TIME("Non-mortar and mortar segment");
		list<pair<long,long>> tempList_0;
		for(long ti = 0; ti < vem(0).elemFace.rows(); ti ++){
			for(long tj = 0; tj < vem(0).elemFace(ti).rows(); tj ++){
				bool flag_ij = true;
				for(long tk = 0; tk < 3; tk ++){
					COORDINATE tempCoor;
					tempCoor(0) = vem(0).nodeCoor(vem(0).elemFace(ti)(tj,tk),0);
					tempCoor(1) = vem(0).nodeCoor(vem(0).elemFace(ti)(tj,tk),1);
					tempCoor(2) = vem(0).nodeCoor(vem(0).elemFace(ti)(tj,tk),2);
					map<COORDINATE,Matrix<long,2,3>>::iterator iter_0 
						= poinNumb(0).find(tempCoor);
					if(iter_0 == poinNumb(0).end()){
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
					COORDINATE tempCoor;
					tempCoor(0) = vem(1).nodeCoor(vem(1).elemFace(ti)(tj,tk),0);
					tempCoor(1) = vem(1).nodeCoor(vem(1).elemFace(ti)(tj,tk),1);
					tempCoor(2) = vem(1).nodeCoor(vem(1).elemFace(ti)(tj,tk),2);
					map<COORDINATE,Matrix<long,2,3>>::iterator iter_0 
						= poinNumb(1).find(tempCoor);
					if(iter_0 == poinNumb(1).end()){
						flag_ij = false;
						break;
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
		// OUTPUT_CONTACT();
		//
		OUTPUT_TIME("Spatial search and contact detection");
		Matrix<double,Dynamic,Dynamic> locaCoor;
		locaCoor.resize(nonmSegm.rows(),2);
		for(long ti = 0; ti < nonmSegm.rows(); ti ++){
			Matrix<double,1,3> tempMatr = (vem(0).nodeCoor.block(nonmSegm(ti,0),0,1,3)
				+ vem(0).nodeCoor.block(nonmSegm(ti,1),0,1,3)
				+ vem(0).nodeCoor.block(nonmSegm(ti,2),0,1,3)) / 3.0;
			locaCoor(ti,0) = tempMatr(0,1);
			locaCoor(ti,1) = sqrt(pow(tempMatr(0,0) + (a_h2 + centErro), 2.0)
				+ pow(tempMatr(0,2), 2.0));
		}
		MatrixXd tempFeat;
		tempFeat.resize(2,1);
		buckNumb << gridNumb(0,4) / 2 * pow2[refiLeve - 1], 
			gridNumb(0,3) / 2 * pow2[refiLeve - 1];
		Matrix<double,2,2> miniMaxi;
		miniMaxi(0,0) = locaCoor.block(0,0,locaCoor.rows(),1).minCoeff();
		miniMaxi(0,1) = locaCoor.block(0,0,locaCoor.rows(),1).maxCoeff();
		miniMaxi(1,0) = locaCoor.block(0,1,locaCoor.rows(),1).minCoeff();
		miniMaxi(1,1) = locaCoor.block(0,1,locaCoor.rows(),1).maxCoeff();
		tempFeat << (miniMaxi(0,1) - miniMaxi(0,0)) / buckNumb(0), 
			(miniMaxi(1,1) - miniMaxi(1,0)) / buckNumb(1);
		BUCKET_SORT(locaCoor, tempFeat, buckNumb(0), buckNumb(1));
		locaCoor.resize(mortSegm.rows(),2);
		for(long ti = 0; ti < mortSegm.rows(); ti ++){
			Matrix<double,1,3> tempMatr = (vem(1).nodeCoor.block(mortSegm(ti,0),0,1,3)
				+ vem(1).nodeCoor.block(mortSegm(ti,1),0,1,3)
				+ vem(1).nodeCoor.block(mortSegm(ti,2),0,1,3)) / 3.0;
			locaCoor(ti,0) = tempMatr(0,1);
			locaCoor(ti,1) = sqrt(pow(tempMatr(0,0) + (a_h2 + centErro), 2.0)
				+ pow(tempMatr(0,2), 2.0));
		}
		CONTACT_SEARCH(locaCoor);
		OUTPUT_CONTACT_ELEMENT();
		//sliding friction
		Vector3d omeg_h, omeg_2;
		omeg_h << 0.0, 1.0, 0.0;
		omeg_2 << 0.0, 0.0, i_2h;
		Matrix<Vector3d,Dynamic,Dynamic> fricDire;
		fricDire.resize(contElem.rows(), 1);
		for(long ti = 0; ti < contElem.rows(); ti ++){
			Vector3d posi_h, posi_2, velo_h, velo_2, velo_r;
			posi_h = contElem(ti).contPoin.block(0,0,3,1);
			posi_2 = contElem(ti).contPoin.block(0,1,3,1);
			posi_h(0) = posi_h(0) + (a_h2 + centErro);
			velo_h = omeg_h.cross(posi_h);
			velo_2 = omeg_2.cross(posi_2);
			velo_r = velo_h - velo_2;
			fricDire(ti) = fricCoef * velo_r.normalized();
		}
		// fricDire.resize(0,0);
		CONTACT_ANALYSIS(loadList, fricDire);
		//
		for(long ti = 0; ti < buckNumb(0); ti ++){
			for(long tj = 0; tj < buckNumb(1); tj ++){
				buck(ti,tj).clear();
			}
		}
		for(long tv = 0; tv < 2; tv ++){
			vem(tv).CLEAR();
			poinNumb(tv).clear();
			segmVolu(tv).clear();
			for(long ti = 0; ti < buckNumb_R(tv,0); ti ++){
				for(long tj = 0; tj < buckNumb_R(tv,1); tj ++){
					buck_R(tv)(ti,tj).clear();
				}
			}
		}
	}
	return 1;
}
