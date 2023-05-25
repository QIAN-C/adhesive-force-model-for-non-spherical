#ifndef MECHSYS_LBM_USERDATA_H
#define MECHSYS_LBM_USERDATA_H


struct InputData
{
     Array<Cell *> Left;
    Array<Cell *> Right;
    Array<Cell *> LeftR;
    Array<Cell *> RightL;
	
	size_t        nx;
	size_t	      ny;
	size_t        nz;

	double		  dx;
	double	      dt;

	double        Length;
    double        Width;
    double        Depth;
	
	double        rp;
	double        rp2;
	double        rp3;
	
	double        dp;
	double      rho_p;
	double      rho_p2;
	double      rho_p3;
	double		  rho_g;

	double        E;
	double        Sig;
	double        G;
	double        Gama;

	double        E2;
	double        Sig2;
	double        G2;
	double        Gama2;

	double        E3;
	double        Sig3;
	double        G3;
	double        Gama3;

	Vec3_t		  up_in;
	Vec3_t		  uf_in;
        Vec3_t    uf_in_particle;

        size_t 		  Np_tol;
	size_t 		  Np_tol_in;
        size_t 		  Np_tol2;
	size_t 		  Np_tol_in2;
	size_t 		  Np_tol_in3;
	
	size_t 		  Np_cap;
	size_t 		  Np_ecp;
	
	double		  nu;
	size_t   inter_time;
        size_t   inter_time2;
	size_t   inter_number;
	
};

#endif
