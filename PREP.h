
#ifndef _PREP_H
#define _PREP_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <list>
#include <set>
#include <ctime>
#include <complex>
#include <random>
#include <sys/stat.h>
//for windows: <direct.h>

#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "Eigen/Sparse"
#include "interval/interval.hpp"

using namespace std;
using namespace Eigen;

const double PI = acos(-1.0);

//output string with current time
long OUTPUT_TIME(string outpStri);

//different numerical examples, different output directories
string outpDire = "";//extern
//concatenate outpDire with fileName
string DIRECTORY(string fileName);

//power of 2
const long pow2[] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};

//solve cubic equation with coefficients a, b, c, d by Cardano's method
//only real solutions, multiple solution is not specially recognized
long CUBIC(double a, double b, double c, double d, Matrix<double,Dynamic,Dynamic> &solu);

//sort every row or column of matrix
//dim: 0 - row, 1 - column; ascend/descend
MatrixXd SORT(const MatrixXd &x, int dim = 1, bool ascend = true);

long OUTPUT_TIME(string outpStri){
	time_t now = time(0);
	char *cnow = ctime(&now);
	cout<<outpStri<<": "<<cnow<<endl;
	return 1;
}

string DIRECTORY(string fileName){
	stringstream tempStre;
	tempStre << outpDire << fileName;
	return tempStre.str();
}

long CUBIC(double a, double b, double c, double d, Matrix<double,Dynamic,Dynamic> &solu){
	if(a == 0){
		cout<<"ERROR in CUBIC: a == 0"<<endl;
		return -1;
	}
	double p = (3.0 * a * c - b * b) / 3.0 / a / a;
	double q = (27.0 * a * a * d - 9.0 * a * b * c + 2.0 * b * b * b) / 27.0 / a / a / a;
	double Delt_r = (q / 2.0) * (q / 2.0) + (p / 3.0) * (p / 3.0) * (p / 3.0);
	complex<double> root_1, root_2;
	if(Delt_r >= 0.0){
		root_1 = complex(- q / 2.0 - sqrt(Delt_r), 0.0);
		root_2 = complex(- q / 2.0 + sqrt(Delt_r), 0.0);
	}
	else{
		root_1 = complex(- q / 2.0, - sqrt(- Delt_r));
		root_2 = complex(- q / 2.0, + sqrt(- Delt_r));
	}
	double ampl_1 = pow(abs(root_1), 1.0 / 3.0);
	double phas_1 = arg(root_1);
	if(phas_1 < 0.0){
		phas_1 += 2.0 * PI;
	}
	phas_1 /= 3.0;
	double ampl_2 = pow(abs(root_2), 1.0 / 3.0);
	double phas_2 = arg(root_2);
	if(phas_2 < 0.0){
		phas_2 += 2.0 * PI;
	}
	phas_2 /= 3.0;
	long soluNumb = 0;
	Matrix<double,9,1> solu_1;
	for(long ti = 0; ti < 3; ti ++){
		double phas_1ti = phas_1 + ti * 2.0 * PI / 3.0;
		root_1 = complex(ampl_1 * cos(phas_1ti), ampl_1 * sin(phas_1ti));
		for(long tj = 0; tj< 3; tj ++){
			double phas_2ti = phas_2 + tj * 2.0 * PI / 3.0;
			root_2 = complex(ampl_2 * cos(phas_2ti), ampl_2 * sin(phas_2ti));
			complex<double> root_ij = root_1 + root_2;
			complex<double> resu = root_ij * root_ij * root_ij + p * root_ij + q;
			if(abs(resu) < 1.0E-12 && abs(imag(root_ij)) < 1.0E-12){
				double miniDiff = 1.0E16;
				for(long tk = 0; tk < soluNumb; tk ++){
					if(abs(real(root_ij) - solu_1(tk)) < miniDiff){
						miniDiff = abs(real(root_ij) - solu_1(tk));
					}
				}
				if(miniDiff > 1.0E-14){
					solu_1(soluNumb) = real(root_ij);
					soluNumb ++;
				}
			}
		}
	}
	solu.resize(soluNumb, 1);
	for(long ti = 0; ti < soluNumb; ti ++){
		solu(ti) = solu_1(ti) - b / 3.0 / a;
	}
	return 1;
}

MatrixXd SORT(const MatrixXd &x, int dim, bool ascend)
{
    MatrixXd y(x);
    auto pred_ascend = [](double a0, double a1){return a0 < a1;};
    auto pred_descend = [](double a0, double a1){return a0 > a1;};

    if(dim == 0)
    {
        if(ascend)
        {
            for(auto row: y.rowwise())
            {
                sort(row.begin(), row.end(), pred_ascend);
            }
        }
        else
        {
            for(auto row: y.rowwise())
            {
                sort(row.begin(), row.end(), pred_descend);
            }
        }
    }
    else if(dim == 1 || dim == -1)
    {
        if(ascend)
        {
            for(auto col : y.colwise())
            {
                sort(col.begin(), col.end(), pred_ascend);
            }
        }
        else
        {
            for(auto col : y.colwise())
            {
                sort(col.begin(), col.end(), pred_descend);
            }
        }
    }
    return y;
}

#ifdef __linux__
#include <petscmat.h>
#include <petscksp.h>
bool ASSEMBLE(const SparseMatrix<double,RowMajor> &origStif, Mat &stif){
	for(long ti = 0; ti < origStif.rows(); ti ++){
		for(SparseMatrix<double,RowMajor>::InnerIterator iter_S(origStif, ti); 
			iter_S; ++ iter_S){
			PetscCall(MatSetValue(stif,
				iter_S.row(), iter_S.col(), iter_S.value(), INSERT_VALUES));
		}
	}
	
	return true;
}
#endif

#endif
