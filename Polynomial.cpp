//============================================================================
// Name        : Polynomial.cpp
// Author      : sc
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Polynomial.h"

/*struct poly_t
{
	double *a_coeffs;       // array of coefficients
	unsigned int degree;    // the degree of the polynomial
};

void init_poly(poly_t &p, double const init_coeffs[],
	unsigned int const init_degree);
void destroy_poly(poly_t &p);
unsigned int poly_degree(poly_t const &p);
double poly_coeff(poly_t const &p, unsigned int n);
double poly_val(poly_t const &p, double x);
void poly_add( poly_t &p, poly_t const &q );
void poly_multiply( poly_t &p, poly_t const &q );
double poly_divide( poly_t &p, double r );
void poly_diff( poly_t &p );
double poly_approx_int( poly_t const &p, double a, double b, unsigned int n );*/

#ifndef MARMOSET_TESTING
int main();
#endif

#ifndef MARMOSET_TESTING
int main()
{
	/*poly_t data{nullptr,0};
	double d[5]{-5,1,3,-1,2};
	init_poly(data,d,4);
	poly_t data2{nullptr,0};
	double d2[8]{2,12,4,-12,7, 1,3,2};
    init_poly(data2,d2,7);
    poly_t data3{nullptr,0};
    double d3[5]{-5,1,3,-1,-2};
    	init_poly(data3,d3,4);
	std::cout << poly_val(data,3)<<std::endl;
	for (std::size_t k{ 0 }; k < data.degree + 1; ++k) {
				std::cout<<data.a_coeffs[k]<<" ";
			}
	std::cout<<std::endl;
	for (std::size_t k{ 0 }; k < data3.degree + 1; ++k) {
					std::cout<<data3.a_coeffs[k]<<" ";
				}
			std::cout<<std::endl;
	poly_add(data,data3);
	for (std::size_t k{ 0 }; k < data.degree + 1; ++k) {
			std::cout<<data.a_coeffs[k]<<" ";
		}
	std::cout<<std::endl;
	poly_add(data,data3);
		for (std::size_t k{ 0 }; k < data.degree + 1; ++k) {
				std::cout<<data.a_coeffs[k]<<" ";
			}
		std::cout<<std::endl;
	for (std::size_t k{ 0 }; k < data2.degree + 1; ++k) {
				std::cout<<data2.a_coeffs[k]<<" ";
			}
	std::cout<<std::endl;
	poly_add(data,data2);
	for (std::size_t k{ 0 }; k < data.degree + 1; ++k) {
		std::cout<<data.a_coeffs[k]<<" ";
		}
	destroy_poly(data);
	    poly_t data{nullptr,0};
		double d[4]{6,-3,0,7};
		init_poly(data,d,3);

		for (std::size_t k{ 0 }; k < data.degree + 1; ++k) {
			std::cout<<data.a_coeffs[k]<<" ";
		   }
		std::cout<<std::endl;

		std::cout<<std::endl;
		 poly_diff( data );
		for (std::size_t k{ 0 }; k < data.degree + 1; ++k) {
			std::cout<<data.a_coeffs[k]<<" ";
		    }
		std::cout<<std::endl;
		poly_t data2{nullptr,0};
				double d2[1]{9};
				init_poly(data2,d2,0);
				for (std::size_t k{ 0 }; k < data2.degree + 1; ++k) {
							std::cout<<data2.a_coeffs[k]<<" ";
						   }
						std::cout<<std::endl;
						poly_diff( data2 );
		for (std::size_t k{ 0 }; k < data2.degree + 1; ++k) {
					std::cout<<data2.a_coeffs[k]<<" ";
				    }
		destroy_poly(data);
	poly_t data{nullptr,0};
			double d[2]{-4,4};
	init_poly(data,d,1);
	poly_t data2{nullptr,0};
				double d2[2]{1,1};
		init_poly(data2,d2,1);
		poly_multiply( data, data2 );
		for (std::size_t k{ 0 }; k < data.degree + 1; ++k) {
							std::cout<<data.a_coeffs[k]<<" ";
			}*/
	poly_t data{nullptr,0};
				double d[2]{17.5011,0.170403};
		init_poly(data,d,1);
		poly_t data2{nullptr,0};
						double d2[5]{0.0132435,152.889,0.0585781,2.14551,955.224};
				init_poly(data2,d2,4);
				poly_multiply( data, data2 );
				for (std::size_t k{ 0 }; k < data.degree + 1; ++k) {
											std::cout<<data.a_coeffs[k]<<" ";
							}

	return 0;
}
#endif

void init_poly(poly_t &p, double const init_coeffs[],
	unsigned int const init_degree) {
	if (p.a_coeffs != nullptr) {
		destroy_poly(p);
	}
	p.degree = init_degree;
	p.a_coeffs = new double[p.degree + 1];
	for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
		p.a_coeffs[k] = init_coeffs[k];
	}
}

void destroy_poly(poly_t &p) {
	delete[] p.a_coeffs;
	p.a_coeffs = nullptr;
	p.degree = 0;
}

unsigned int poly_degree(poly_t const &p) {
	if (p.a_coeffs == nullptr) {
		throw 0;
	}
	return p.degree;
}

double poly_coeff(poly_t const &p, unsigned int n) {
	if (p.a_coeffs == nullptr) {
		throw 0;
	}
	if (n > p.degree) {
		return 0;
	}
	else {
		return p.a_coeffs[n];
	}
}

double poly_val(poly_t const &p, double x) {
	if (p.a_coeffs == nullptr) {
	    throw 0;
	}
	double ans = 0;
	for (int i = p.degree; i >= 0; i--) {
        ans = p.a_coeffs[i] + x * ans;
	}
	return ans;
}

void poly_add( poly_t &p, poly_t const &q ){
	if (p.a_coeffs == nullptr || q.a_coeffs == nullptr) {
		throw 0;
	}
	poly_t tem{nullptr,0};
	if(q.degree > p.degree){
			tem.degree = q.degree;
			tem.a_coeffs=new double[tem.degree + 1];
				for (std::size_t k{ 0 }; k < tem.degree + 1; ++k) {
							tem.a_coeffs[k] = 0;
						}
				for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
						tem.a_coeffs[k] = p.a_coeffs[k];
					}
				destroy_poly(p);
				p.degree=tem.degree;
				p.a_coeffs=new double[p.degree + 1];
				for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
					 p.a_coeffs[k] = tem.a_coeffs[k] + q.a_coeffs[k];
				 }
		}else if(q.degree < p.degree){
			tem.degree =p.degree;
			tem.a_coeffs=new double[tem.degree + 1];
			for (std::size_t k{ 0 }; k < tem.degree + 1; ++k) {
						   tem.a_coeffs[k] = 0;
					}
			for (std::size_t k{ 0 }; k < q.degree + 1; ++k) {
					  tem.a_coeffs[k] = q.a_coeffs[k];
				}
            for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
	             p.a_coeffs[k] = tem.a_coeffs[k] + p.a_coeffs[k];
    			 }
		}else{
			unsigned int real_degree = p.degree;
			while(p.a_coeffs[real_degree] + q.a_coeffs[real_degree] == 0){
				real_degree--;
			}
			tem.degree = real_degree;
			tem.a_coeffs=new double[tem.degree + 1];
			for (std::size_t k{ 0 }; k < tem.degree + 1; ++k) {
					   tem.a_coeffs[k] = 0;
				}
			for (std::size_t k{ 0 }; k < tem.degree + 1; ++k) {
					tem.a_coeffs[k] = p.a_coeffs[k];
			    }
				destroy_poly(p);
				p.degree=tem.degree;
			    p.a_coeffs=new double[p.degree + 1];
			 for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
				   p.a_coeffs[k] = tem.a_coeffs[k] + q.a_coeffs[k];
			    }
		}
	   destroy_poly(tem);
}

void poly_subtract( poly_t &p, poly_t const &q ){
	if (p.a_coeffs == nullptr || q.a_coeffs == nullptr) {
			throw 0;
		}
		poly_t tem{nullptr,0};
		if(q.degree > p.degree){
				tem.degree = q.degree;
				tem.a_coeffs=new double[tem.degree + 1];
					for (std::size_t k{ 0 }; k < tem.degree + 1; ++k) {
								tem.a_coeffs[k] = 0;
							}
					for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
							tem.a_coeffs[k] = p.a_coeffs[k];
						}
					destroy_poly(p);
					p.degree=tem.degree;
					p.a_coeffs=new double[p.degree + 1];
					for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
						 p.a_coeffs[k] = tem.a_coeffs[k] - q.a_coeffs[k];
					 }
			}else if(q.degree < p.degree){
				tem.degree =p.degree;
				tem.a_coeffs=new double[tem.degree + 1];
				for (std::size_t k{ 0 }; k < tem.degree + 1; ++k) {
							   tem.a_coeffs[k] = 0;
						}
				for (std::size_t k{ 0 }; k < q.degree + 1; ++k) {
						  tem.a_coeffs[k] = q.a_coeffs[k];
					}
	            for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
		             p.a_coeffs[k] = p.a_coeffs[k] - tem.a_coeffs[k];
	    			 }
			}else{
				unsigned int real_degree = p.degree;
				while(p.a_coeffs[real_degree] - q.a_coeffs[real_degree] == 0){
					real_degree--;
				}
				tem.degree = real_degree;
				tem.a_coeffs=new double[tem.degree + 1];
				for (std::size_t k{ 0 }; k < tem.degree + 1; ++k) {
						   tem.a_coeffs[k] = 0;
					}
				for (std::size_t k{ 0 }; k < tem.degree + 1; ++k) {
						tem.a_coeffs[k] = p.a_coeffs[k];
				    }
					destroy_poly(p);
					p.degree=tem.degree;
				    p.a_coeffs=new double[p.degree + 1];
				 for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
					   p.a_coeffs[k] = tem.a_coeffs[k] - q.a_coeffs[k];
				    }
			}
		   destroy_poly(tem);
}

void poly_multiply( poly_t &p, poly_t const &q ){
	if (p.a_coeffs == nullptr || q.a_coeffs == nullptr) {
				throw 0;
			}
			poly_t tem{nullptr,0};
			if((p.a_coeffs[0]==0 && p.degree==0) || (q.a_coeffs[0]==0 && q.degree==0)){
				destroy_poly(p);
			    p.degree=0;
			    p.a_coeffs=new double[p.degree + 1];
			    p.a_coeffs[0]=0;
			}else{
				tem.degree = p.degree + q.degree;
				tem.a_coeffs=new double[tem.degree + 1];
				for (std::size_t k{ 0 }; k < tem.degree + 1; ++k) {
						tem.a_coeffs[k] = 0.0;
					}
			    for (std::size_t i = 0; i < p.degree + 1; i++){
					for (std::size_t j = 0; j < q.degree + 1; j++){
					    tem.a_coeffs[i+j] += p.a_coeffs[i] * q.a_coeffs[j];
					    }
					}
					destroy_poly(p);
				    p.degree=tem.degree;
					p.a_coeffs=new double[p.degree + 1];
				for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
				    p.a_coeffs[k] = tem.a_coeffs[k];
				    }
			}
			destroy_poly(tem);
}

double poly_divide( poly_t &p, double r ){
	if (p.a_coeffs == nullptr) {
		    throw 0;
		}
	if(p.degree > 0){
	poly_t tem{nullptr,0};
	tem.degree = p.degree;
    tem.a_coeffs=new double[tem.degree + 1];
	 int flag = p.degree;
	 for (int k = flag; k >= 0; --k){
		 if(k == flag){
			 tem.a_coeffs[k] = p.a_coeffs[k];
		   }else{
			 tem.a_coeffs[k] = p.a_coeffs[k] + tem.a_coeffs[k+1] * r;
		   }
		 }
	  double d=0.0;
	  d = tem.a_coeffs[0];
	  destroy_poly(p);
	  p.degree=tem.degree - 1;
	  p.a_coeffs=new double[p.degree + 1];
	  for (int k = p.degree; k >= 0; --k) {
	       p.a_coeffs[k] = tem.a_coeffs[k+1];
	       }
	  destroy_poly(tem);
	  return d;
	}else{
		double d=p.a_coeffs[0];
		p.a_coeffs[0]=0;
		return d;
	}
}

void poly_diff( poly_t &p ){
	if (p.a_coeffs == nullptr) {
			    throw 0;
	}
	poly_t tem{nullptr,0};
	if(p.degree > 0){
	tem.degree = p.degree-1;
	tem.a_coeffs=new double[tem.degree + 1];
	 for (int k = tem.degree; k >= 0; --k) {
		tem.a_coeffs[k] = p.a_coeffs[k+1]*(k+1);
	    }
	 destroy_poly(p);
	 p.degree=tem.degree;
	 p.a_coeffs=new double[p.degree + 1];
	 for (std::size_t k{ 0 }; k < p.degree + 1; ++k) {
	        p.a_coeffs[k] = tem.a_coeffs[k];
	      }
	}else{
		p.a_coeffs[0] = 0;
	}
	 destroy_poly(tem);
}

double poly_approx_int( poly_t const &p, double a, double b, unsigned int n ){
	if (p.a_coeffs == nullptr) {
	     throw 0;
		}
	double sum=0.0;
	double h=(b-a)/n;
	for (std::size_t k{ 0 }; k <= n; ++k) {
        if(k==0 || k==n){
        	 sum = sum + poly_val(p, a + k*h);
        }else{
        sum = sum + 2*poly_val(p, a + k*h);
        }
		  }
	return h / 2 * sum;
}






