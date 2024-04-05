
#include "DEHW.h"
#include "BEAM.h"
#include "BLOCK.h"
#include "CYLINDER.h"

long SHOW_SLEC(DEHW &solu);
long GEOMETRY();
long GRID();
long LOAD();
long PATC_TEST();
long CYLI_HERT();
long PETSC_PREPARE(int argc, char **argv);

int main(int argc, char **argv){
	PETSC_PREPARE(argc, argv);
	cout<<"**********************************************************"<<endl;
	cout<<"*1 - DEHW drive: ICLOSE, MLLOFE, CILLOFE, MLLOSE, CILLOSE*"<<endl;
	cout<<"*2 - DEHW drive: grid discretization of tooth surface    *"<<endl;
	cout<<"*3 - DEHW drive: loaded contact analysis                 *"<<endl;
	cout<<"*4 - Static analysis of flexible beam                    *"<<endl;
	cout<<"*5 - Patch test of two blocks                            *"<<endl;
	cout<<"*6 - Loaded contact analysis of two cylindrical bodies   *"<<endl;
	cout<<"********************************************************"<<endl;
	cout<<"Please enter a number in 1~6: ";
	long meth;
	cin >> meth;
	switch(meth){
		case 1:{
			GEOMETRY();
			break;
		}
		case 2:{
			GRID();
			break;
		}
		case 3:{
			LOAD();
			break;
		}
		case 4:{
			BEAM solu;
			solu.SOLVE();
			break;
		}
		case 5:{
			PATC_TEST();
			break;
		}
		case 6:{
			CYLI_HERT();
			break;
		}
	}
#ifdef __linux__
    PetscCall(PetscFinalize());
#endif
	return 1;
}

long SHOW_SLEC(DEHW &solu){
	cout<<"*1 - Case 1, z/z=1/40, Δa=0, Δi=0                        *"<<endl;
	cout<<"*2 - Case 2, z/z=1/40, Δa=0.1, Δi=0                      *"<<endl;
	cout<<"*3 - Case 3, z/z=1/40, Δa=0.2, Δi=0.0014                 *"<<endl;
	cout<<"*4 - Case 4, z/z=1/45                                    *"<<endl;
	cout<<"*5 - Case 5, z/z=3/54                                    *"<<endl;
	cout<<"**********************************************************"<<endl;
	cout<<"Please enter a number in 1~5: ";
	long caid;
	cin >> caid;
	switch(caid){
		case 1:{
			solu.z << 1, 40;
			solu.a_h2 = 0.25;//m
			solu.modiTran = 0.0;
			solu.modiCent = 0.0;
			solu.r_b2 = 0.158/2.0;//m
			solu.beta_c = 11.0 * PI / 180.0;//rad
			//
			solu.z_k = 4.2;
			solu.d(0) = 0.082;//m
			solu.h_a_s << 0.6, 0.7;//smaller, unused
			solu.h_f_s << 0.95, 1.05;//larger, larger
			solu.R_a(1) = 0.0385;//m
			solu.offsR_a = 0.003;//m
			solu.wheeWidt = 0.06;//m
			solu.inneRadi << 0.018, 0.15;//m
			solu.inpuTorq = 180.0;//N*m
			solu.fricCoef = 0.08;
			break;
		}
		case 2:{
			solu.z << 1, 40;
			solu.a_h2 = 0.25;
			solu.modiTran = 0.1;
			solu.modiCent = 0.0;
			solu.r_b2 = 0.158/2.0;
			solu.beta_c = 11.0 * PI / 180.0;
			//
			solu.z_k = 4.2;
			solu.d(0) = 0.082;
			solu.h_a_s << 0.6, 0.7;//smaller, unused
			solu.h_f_s << 0.95, 1.05;//larger, larger
			solu.R_a(1) = 0.0385;
			solu.offsR_a = 0.003;
			solu.wheeWidt = 0.06;
			solu.inneRadi << 0.018, 0.15;
			solu.inpuTorq = 180.0;//N*m
			solu.fricCoef = 0.08;
			break;
		}
		case 3:{
			solu.z << 1, 40;
			solu.a_h2 = 0.25;
			solu.modiTran = 0.2;
			solu.modiCent = 0.0014;
			solu.r_b2 = 0.158/2.0;
			solu.beta_c = 11.0 * PI / 180.0;
			//
			solu.z_k = 4.2;
			solu.d(0) = 0.082;
			solu.h_a_s << 0.6, 0.7;//smaller, unused
			solu.h_f_s << 0.95, 1.05;//larger, larger
			solu.R_a(1) = 0.0385;
			solu.offsR_a = 0.003;
			solu.wheeWidt = 0.06;
			solu.inneRadi << 0.018, 0.15;
			solu.inpuTorq = 180.0;//N*m
			solu.fricCoef = 0.08;
			break;
		}
		case 4:{
			solu.z << 1, 45;
			solu.a_h2 = 0.16;
			solu.modiTran = 0.75;
			solu.modiCent = 0.0025;
			solu.r_b2 = 0.096/2.0;
			solu.beta_c = 11.0 * PI / 180.0;
			//
			solu.z_k = 4.7;
			solu.d(0) = 0.056;
			solu.h_a_s << 0.8, 0.7;//common, unused
			solu.h_f_s << 1.00, 1.00;//common, larger
			solu.R_a(1) = 0.024;
			solu.offsR_a = 0.0;
			solu.wheeWidt = 0.045;
			solu.inneRadi << 0.013, 0.105;
			solu.inpuTorq = 67.0;//N*m
			solu.fricCoef = 0.08;
			break;
		}
		case 5:{
			solu.z << 3, 54;
			solu.a_h2 = 0.28;
			solu.modiTran = 0.1;
			solu.modiCent = 0.001;
			solu.r_b2 = 0.18/2.0;
			solu.beta_c = 13.0 * PI / 180.0;
			//
			solu.z_k = 5.4;
			solu.d(0) = 0.13;
			solu.h_a_s << 0.8, 0.7;//common, unused
			solu.h_f_s << 1.00, 1.00;//common, larger
			solu.R_a(1) = 0.061;
			solu.offsR_a = 0.001;
			solu.wheeWidt = 0.08;
			solu.inneRadi << 0.035, 0.17;
			solu.inpuTorq = 700.0;//N*m
			solu.fricCoef = 0.08;
			break;
		}
	}
	return 1;
}

long GEOMETRY(){
	cout<<endl;
	cout<<"***DEHW drive: ICLOSE, MLLOFE, CILLOFE, MLLOSE, CILLOSE***"<<endl;
	DEHW solu;
	SHOW_SLEC(solu);
	solu.BASIC_PARAMETER();
	solu.GEOMETRY();
	return 1;
}

long GRID(){
	cout<<endl;
	cout<<"*****DEHW drive: grid discretization of tooth surface*****"<<endl;
	DEHW solu;
	SHOW_SLEC(solu);
	solu.BASIC_PARAMETER();
	solu.refiLeve = 4;
	solu.TS_GRID();
	// solu.ANALYTICAL_CONTACT();//accuracy not enough......
	// double sigm_H, b_H;
	// solu.AC_PRESSURE(6.292972372745343, 6.283185307179586, sigm_H, b_H);
	return 1;
}

long LOAD(){
	cout<<endl;
	cout<<"************DEHW drive: loaded contact analysis***********"<<endl;
	DEHW solu;
	SHOW_SLEC(solu);
	//
	cout<<endl;
	cout<<"**************DEHW drive: loaded contact analysis**************"<<endl;
	cout<<"*1 - Without tooth flank relief, without center distance error*"<<endl;
	cout<<"*2 - With tooth flank relief, without center distance error   *"<<endl;
	cout<<"*3 - Without tooth flank relief, with center distance error   *"<<endl;
	cout<<"*4 - With tooth flank relief, with center distance error      *"<<endl;
	cout<<"***************************************************************"<<endl;
	cout<<"Please enter a number in 1~4: ";
	long caid;
	cin >> caid;
	switch(caid){
		case 1:{
			solu.reliSwit = 0;
			solu.centErro = 0.0E-6;
			solu.distCrit << 38.0E-6, 34.0E-6, 30.0E-6, 26.0E-6, 22.0E-6, 18.0E-6;
			break;
		}
		case 2:{
			solu.reliSwit = 1;
			solu.centErro = 0.0E-6;
			solu.distCrit << 42.0E-6, 38.0E-6, 34.0E-6, 30.0E-6, 26.0E-6, 22.0E-6;
			break;
		}
		case 3:{
			solu.reliSwit = 0;
			solu.centErro = 10.0E-6;
			solu.distCrit << 42.0E-6, 38.0E-6, 34.0E-6, 30.0E-6, 26.0E-6, 22.0E-6;
			break;
		}
		case 4:{
			solu.reliSwit = 1;
			solu.centErro = 10.0E-6;
			solu.distCrit << 50.0E-6, 46.0E-6, 42.0E-6, 38.0E-6, 34.0E-6, 30.0E-6;
			break;
		}
	}
	//
	solu.BASIC_PARAMETER();
	solu.SOLVE();
	return 1;
}

long PATC_TEST(){
	cout<<endl;
	cout<<"*****************Patch test of two blocks*****************"<<endl;
	cout<<"*1 - Case I with sparse mesh                             *"<<endl;
	cout<<"*2 - Case I with refined mesh                            *"<<endl;
	cout<<"*3 - Case II with sparse mesh                            *"<<endl;
	cout<<"*4 - Case II with refined mesh                           *"<<endl;
	cout<<"*5 - Case III with sparse mesh                           *"<<endl;
	cout<<"*6 - Case III with refined mesh                          *"<<endl;
	cout<<"**********************************************************"<<endl;
	cout<<"Please enter a number in 1~6: ";
	long caid;
	cin >> caid;
	BLOCK solu(caid - 1);
	solu.SOLVE();
	return 1;
}

long CYLI_HERT(){
	cout<<endl;
	cout<<"****Loaded contact analysis of two cylindrical bodies*****"<<endl;
	cout<<"*1 - Case I, same length                                 *"<<endl;
	cout<<"*2 - Case II, different length, without corner rounding  *"<<endl;
	cout<<"*3 - Case III, different length, with corner rounding    *"<<endl;
	cout<<"**********************************************************"<<endl;
	cout<<"Please enter a number in 1~3: ";
	long caid;
	cin >> caid;
	CYLINDER solu(caid - 1);
	solu.SOLVE();
	return 1;
}

long PETSC_PREPARE(int argc, char **argv){
#ifdef __linux__
    char help[] = "Sparse Linear Solver of PETSc.\n\n";
    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, (char *)0, help));
    PetscMPIInt size;
    PetscCall(MPI_Comm_size(PETSC_COMM_WORLD,&size));
#endif
	return 1;
}
