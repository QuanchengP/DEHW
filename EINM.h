
#ifndef _EINM_H
#define _EINM_H

#include "PREP.h"

//Simplify instantiation of intervals
typedef filib::interval<double, filib::native_switched, filib::i_mode_extended> FINTERVAL;
//Simplify access to methods of class fp_traits<>
typedef filib::fp_traits<FINTERVAL::value_type> TRAITS;

//extended interval Newton method
template<typename T>
class EINM{
public:
	FINTERVAL (*F) (const FINTERVAL&, const T&);
	FINTERVAL (*DF) (const FINTERVAL&, const T&); //Derivative df of function f
	T fuda;
	FINTERVAL searRang; //search range
	double epsl_rX;
	double epsl_X;
	double epsl_f; //Stopping criteria
	//
	list<pair<long, FINTERVAL>> done; //list of solutions
	//
	EINM();
	long FM0DFN0(FINTERVAL X, long veri); //0 is in F(mid(X)), but not in DF(X)
	long FM0DF0(FINTERVAL X, long &f_poe); //0 is in both F(mid(X)) and DF(X)
	//one iteration step when 0 is not in DF(X)
	long STEP_DFN0(double poe, FINTERVAL &X, long &veri);
	//one iteration step when 0 is not in F(poe), but in DF(X)
	long STEP_FN0DF0(double poe, FINTERVAL X, list<FINTERVAL> &X_NEW);
	//if two solution intervals intersect, then merge them
	long ABSORB(long veri, const FINTERVAL &X);
	long SOLVE();//extended interval Newton method
	long OUTPUT();
};

template<typename T>
EINM<T>::EINM(){
	filib::fp_traits<double, filib::native_switched>::setup(); //Do initializations
	
	searRang = FINTERVAL::ENTIRE();
}

template<typename T>
long EINM<T>::FM0DFN0(FINTERVAL X, long veri){
	double xtil = mid(X);
	long f_a = 0;
	long f_b = 0;
	long n = 0;
	while(true){
		if(f_a != 1){
			FINTERVAL FA = (*F)(X.inf(), fuda);
			if(filib::in(0.0, FA)){
				f_a = 1;
			}
			else{
				STEP_DFN0(X.inf(), X, veri);
				if(filib::isEmpty(X)){
					break;
				}
			}
		}
		if(f_b != 1){
			FINTERVAL FB = (*F)(X.sup(), fuda);
			if(filib::in(0.0, FB)){
				f_b = 1;
			}
			else{
				if(filib::in(xtil, X)){
					STEP_DFN0(X.sup(), X, veri);
				}
				else{
					STEP_DFN0(mid(X), X, veri);
				}
				if(filib::isEmpty(X)){
					break;
				}
			}
		}
		if(f_a == 1 && f_b == 1){
			if(veri == 1){
				ABSORB(1, X);
			}
			else{
				ABSORB(2, X);
			}
			break;
		}
		STEP_DFN0(mid(X), X, veri);
		if(filib::isEmpty(X)){
			break;
		}
		n ++;
		if(n >= 4 && filib::in(0.0, (*F)(X, fuda))){
			if(veri == 1){
				ABSORB(1, X);
			}
			else{
				ABSORB(3, X);
			}
			break;
		}
	}
	return 1;
}

template<typename T>
long EINM<T>::FM0DF0(FINTERVAL X, long &f_poe){
	if(!filib::in(0.0, (*F)(FINTERVAL(X.inf()), fuda))){
		f_poe = 1;
	}
	else if(!filib::in(0.0, (*F)(FINTERVAL(X.sup()), fuda))){
		f_poe = 2;
	}
	else if(!filib::in(0.0, (*F)(FINTERVAL((3.0*X.inf()+X.sup())/4.0), fuda))
		|| !filib::in(0.0, (*F)(FINTERVAL((X.inf()+3.0*X.sup())/4.0), fuda))){
		f_poe = 3;
	}
	else{
		ABSORB(5, X);
		f_poe = 4;
	}
	return 1;
}

template<typename T>
long EINM<T>::STEP_DFN0(double poe, FINTERVAL &X, long &veri){
	FINTERVAL X_ = (poe - (*F)(FINTERVAL(poe), fuda) / (*DF)(X, fuda)).intersect(X);
	if(!filib::isEmpty(X_) && X_.interior(X)){
		veri = 1;
	}
	X = X_;
	return 1;
}

template<typename T>
long EINM<T>::STEP_FN0DF0(double poe, FINTERVAL X, list<FINTERVAL> &X_NEW){
	FINTERVAL FPOE = (*F)(FINTERVAL(poe), fuda);
	FINTERVAL DFX = (*DF)(X, fuda);
	if(FPOE.inf() > 0.0){
		double p_n = poe - FPOE.inf() / DFX.inf();
		double q_n = poe - FPOE.inf() / DFX.sup();
		if(DFX.inf() == 0.0){
			FINTERVAL X_NEW_S = (FINTERVAL(TRAITS::ninfinity(), q_n)).intersect(X);
			if(!filib::isEmpty(X_NEW_S)){
				X_NEW.push_back(X_NEW_S);
			}
		}
		else if(DFX.sup() == 0.0){
			FINTERVAL X_NEW_S = (FINTERVAL(p_n, TRAITS::infinity())).intersect(X);
			if(!filib::isEmpty(X_NEW_S)){
				X_NEW.push_back(X_NEW_S);
			}
		}
		else{
			FINTERVAL X_NEW_S = (FINTERVAL(TRAITS::ninfinity(), q_n)).intersect(X);
			//when DFX.inf() and DFX.sup() are too large, it is possible that p_n == q_n
			if(!filib::isEmpty(X_NEW_S) && !(p_n == q_n && filib::diam(X_NEW_S) == 0)){
				X_NEW.push_back(X_NEW_S);
			}
			X_NEW_S = (FINTERVAL(p_n, TRAITS::infinity())).intersect(X);
			if(!filib::isEmpty(X_NEW_S) && !(p_n == q_n && filib::diam(X_NEW_S) == 0)){
				X_NEW.push_back(X_NEW_S);
			}
		}
	}
	else if(FPOE.sup() < 0.0){
		double p_n = poe - FPOE.sup() / DFX.inf();
		double q_n = poe - FPOE.sup() / DFX.sup();
		if(DFX.inf() == 0.0){
			FINTERVAL X_NEW_S = (FINTERVAL(q_n, TRAITS::infinity())).intersect(X);
			if(!filib::isEmpty(X_NEW_S)){
				X_NEW.push_back(X_NEW_S);
			}
		}
		else if(DFX.sup() == 0.0){
			FINTERVAL X_NEW_S = (FINTERVAL(TRAITS::ninfinity(), p_n)).intersect(X);
			if(!filib::isEmpty(X_NEW_S)){
				X_NEW.push_back(X_NEW_S);
			}
		}
		else{
			FINTERVAL X_NEW_S = (FINTERVAL(TRAITS::ninfinity(), p_n)).intersect(X);
			if(!filib::isEmpty(X_NEW_S) && !(p_n == q_n && filib::diam(X_NEW_S) == 0)){
				X_NEW.push_back(X_NEW_S);
			}
			X_NEW_S = (FINTERVAL(q_n, TRAITS::infinity())).intersect(X);
			if(!filib::isEmpty(X_NEW_S) && !(p_n == q_n && filib::diam(X_NEW_S) == 0)){
				X_NEW.push_back(X_NEW_S);
			}
		}
	}
	else{
		cout<<"ERROR in EINM::STEP!"<<endl;
	}
	return 1;
}

template<typename T>
long EINM<T>::ABSORB(long veri, const FINTERVAL &X){
	//"done" is a sorted list without overlapping intervals
	list<pair<long, FINTERVAL>>::iterator iter_d = done.begin();
	bool abso = false;
	while(iter_d != done.end()){
		if(!filib::disjoint((*iter_d).second, X)){
			//Overlapping intervals
			(*iter_d).second = filib::hull((*iter_d).second, X);
			if(veri == 1 || (*iter_d).first == 1){
				(*iter_d).first = 1;
			}
			else if(2 <= veri && veri <=3 && 
				2 <= (*iter_d).first && (*iter_d).first<=3){
				(*iter_d).first = 4;
			}
			else if(veri == 5 && (*iter_d).first == 5){
				(*iter_d).first = 5;
			}
			else if(veri == 6 && (*iter_d).first == 6){
				(*iter_d).first = 6;
			}
			else{
				(*iter_d).first = 7;
			}
			abso = true;
			break;
		}
		if(X.inf() < (*iter_d).second.inf()){
			//Sort from left to right
			done.insert(iter_d, pair<long, FINTERVAL>(veri, X));
			abso = true;
			break;
		}
		iter_d ++;
	}
	if(!abso){
		done.push_back(pair<long, FINTERVAL>(veri, X));
	}
	
	return 1;
}

template<typename T>
long EINM<T>::SOLVE(){
	list<FINTERVAL> todo;
	todo.push_back(searRang);
	while(!todo.empty()){
		FINTERVAL X = todo.back();
		todo.pop_back();
		if(!filib::in(0.0, (*F)(X, fuda))){
			continue;
		}
		if(!filib::in(0.0, (*DF)(X, fuda))){
			long veri = 0;
			while(true){
				STEP_DFN0(mid(X), X, veri);
				if(filib::isEmpty(X)){
					break;
				}
				else if(filib::in(0.0, (*F)(FINTERVAL(mid(X)), fuda))){
					FM0DFN0(X, veri);
					break;
				}
			}
		}
		else{
			long f_poe = 0;
			if(!filib::in(0.0, (*F)(FINTERVAL(mid(X)), fuda))){
				if(filib::mag((*F)(X, fuda)) < epsl_f &&
					((!filib::in(0.0, X) && filib::diam(X) / filib::mag(X) < epsl_rX) 
					|| (filib::in(0.0, X) && filib::diam(X) < epsl_X))){
					ABSORB(6, X);
					continue;
				}
			}
			else{
				FM0DF0(X, f_poe);
				if(f_poe == 3){
					todo.push_back(FINTERVAL(X.inf(), mid(X)));
					todo.push_back(FINTERVAL(mid(X), X.sup()));
				}
				if(f_poe == 3 || f_poe == 4){
					continue;
				}
			}
			double poe;
			if(f_poe == 0){
				poe = mid(X);
			}
			else if(f_poe == 1){
				poe = X.inf();
			}
			else if(f_poe == 2){
				poe = X.sup();
			}
			else{
				cout<<"ERROR in EINM::SOLVE!"<<endl;
			}
			list<FINTERVAL> X_NEW;
			STEP_FN0DF0(poe, X, X_NEW);
			if(X_NEW.size() == 0){
				continue;
			}
			else if(X_NEW.size() == 2){
				todo.push_back(X_NEW.back());
				X_NEW.pop_back();
				todo.push_back(X_NEW.back());
				X_NEW.pop_back();
				continue;
			}
			FINTERVAL X_NEW_S = X_NEW.back();
			if(filib::diam(X_NEW_S) < filib::diam(X) / 2.0){
				todo.push_back(X_NEW_S);
			}
			else{
				todo.push_back(FINTERVAL(X.inf(), mid(X)));
				todo.push_back(FINTERVAL(mid(X), X.sup()));
			}
		}
	}
	
	return 1;
}

template<typename T>
long EINM<T>::OUTPUT(){
	cout<<"There may be zeros in the interval(s):"<<endl;
	for(list<pair<long, FINTERVAL>>::const_iterator iter_d = done.begin(); 
		iter_d != done.end(); iter_d ++){
		cout<<(*iter_d).first<<" - "<<(*iter_d).second<<endl;
	}
	
	return 1;
}

#endif
