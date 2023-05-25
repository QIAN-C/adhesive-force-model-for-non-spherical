/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_DEM_INTERACTON_H
#define MECHSYS_DEM_INTERACTON_H

// Std Lib
#include <math.h>
//#include <map>
#include <unordered_map>
#include <vector>
#include <utility>

#include <iostream>
#include <complex>
using namespace std;


// MechSys
#include <mechsys/dem/particle.h>

#include<fstream>
#include<iostream>



namespace DEM
{

// typedefs
typedef std::unordered_map<size_t,Vec3_t> FrictionMap_t;
//typedef std::map<std::pair<int,int>,Vec3_t> FrictionMap_t;
typedef Array<std::pair<int,int> > ListContacts_t;

class Interacton   //General class for interactons
{
public:

    // Constructor and destructor
    Interacton () {}           ///< Default constructor
    virtual ~Interacton () {}  ///< Destructor

    // Methods
    virtual bool UpdateContacts   (double alpha) =0;    ///< Update contacts by verlet algorithm
    virtual bool CalcForce        (double dt = 0.0) =0; ///< Calculates the contact force between particles
    virtual void UpdateParameters () =0;                ///< Update the parameters

    // Data
    Particle     * __restrict__ P1;        ///< First particle
    Particle     * __restrict__ P2;        ///< Second particle
    size_t         I1;        ///< Index of the first particle
    size_t         I2;        ///< Index of the second particle
    Vec3_t         F1;        ///< Provisional force  for particle 1
    Vec3_t         F2;        ///< Provisional force  for particle 2
    Vec3_t         T1;        ///< Provisional torque for particle 1
    Vec3_t         T2;        ///< Provisional torque for particle 2
#ifdef USE_THREAD
    pthread_mutex_t lck;              ///< to protect variables in multithreading
#endif
};

class CInteracton: public Interacton // Interacton for collision
{
public:
    // Constructor
    CInteracton (Particle * Pt1, Particle * Pt2); ///< Constructor requires pointers to both particles
    CInteracton () {};

    // Methods
    virtual bool UpdateContacts (double alpha);    ///< Update contacts by verlet algorithm
    virtual bool CalcForce      (double dt = 0.0); ///< Calculates the contact force between particles
    virtual void UpdateParameters ();              ///< Update the parameters

    // Data
     double     R;	 ///< diameter ***
    double     E;	 ///< elastic module ***
    double     Sig;      ///< Possion's ratio ***
    double     G;	 ///< shear module ***
    double     Gama;     ///< Surface energy ***
    Vec3_t     x_t;          ///< slide distance ***
    Vec3_t     x_t2;          ///< slide distance ***
    Vec3_t     thei_t ;   ///< rolling distance ***
     Vec3_t     thei_tw;     //twisiting distance
    Vec3_t	   SFr; 	 ///< Vector of static friction*************
    
    bool           First;     ///< Is it the first collision?
    double         Kn;        ///< Normal stiffness
    double         Kt;        ///< Tengential stiffness
    double         Gn;        ///< Normal viscous coefficient
    double         Gt;        ///< Tangential viscous coefficient
    double         Mu;        ///< Microscpic coefficient of friction
    double         Epot;      ///< Potential elastic energy
    double         dEvis;     ///< Energy dissipated in viscosity at time step
    double         dEfric;    ///< Energy dissipated by friction at time step
    size_t         Nc;        ///< Number of contacts
    size_t         Nsc;       ///< Number of sliding contacts
    size_t         Nr;        ///< Number of rolling contacts (only for spheres)
    Vec3_t         Fn;        ///< Normal force between elements
    Vec3_t         Fnet;      ///< Net normal force
    Vec3_t         Ftnet;     ///< Net tangential force
    Vec3_t         Xc;        ///< Net Position of the contact
    Vec3_t         Branch;    ///< Branch vector
    
	  Vec3_t		  Fdvv; 						  ///< Static Friction displacement for the vertex vertex pair***********
	  Vec3_t		 Fdr;							 ///< Rolling displacement *******
	  double		 beta;							 ///< Rolling stiffness coefficient**********
	  double		 eta;							 ///< Plastic moment coefficient**********

	
    ListContacts_t Lee;       ///< List of edge-edge contacts 
    ListContacts_t Lvf;       ///< List of vertex-face contacts 
    ListContacts_t Lfv;       ///< List of face-vertex contacts
    ListContacts_t Lvt;       ///< List of vertex-torus contacts
    ListContacts_t Ltv;       ///< List of torus-vertex contacts
    ListContacts_t Lvc;       ///< List of vertex-cylinder contacts
    ListContacts_t Lcv;       ///< List of cylinder-vertex contacts
      ListContacts_t Lvv;       ///< List of cylinder-vertex contacts**************************
      ListContacts_t Lve;       ///< List of cylinder-vertex contacts**************************
      ListContacts_t Lev;       ///< List of cylinder-vertex contacts**************************
    FrictionMap_t  Fdee;      ///< Static friction displacement for pair of edges
    FrictionMap_t  Fdvf;      ///< Static friction displacement for pair of vertex-face
    FrictionMap_t  Fdfv;      ///< Static friction displacement for pair of face-vertex
    FrictionMap_t  Fdvt;      ///< Static friction displacement for pair of vertex-torus
    FrictionMap_t  Fdtv;      ///< Static friction displacement for pair of torus-vertex
    FrictionMap_t  Fdvc;      ///< Static friction displacement for pair of vertex-cylinder
    FrictionMap_t  Fdcv;      ///< Static friction displacement for pair of cylinder-vertex
     FrictionMap_t  Fdv_v;      ///< Static friction displacement for pair of cylinder-vertex*****************************
     FrictionMap_t  Fdve;      ///< Static friction displacement for pair of cylinder-vertex*****************************
     FrictionMap_t  Fdev;      ///< Static friction displacement for pair of cylinder-vertex*****************************
protected:
    template<typename FeatureA_T, typename FeatureB_T>
    bool _update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, ListContacts_t & L, double dt, double ee);
    template<typename FeatureA_T, typename FeatureB_T>
    void _update_contacts        (FeatureA_T & A, FeatureB_T & B, ListContacts_t & L, double alpha);
};

class CInteractonSphere: public CInteracton // Collision interacton for spheres
{
public:
    // Methods 
    CInteractonSphere (Particle * Pt1, Particle * Pt2); ///< Constructor requires pointers to both particles
    bool UpdateContacts (double alpha);                 ///< Update contacts by verlet algorithm
    bool CalcForce (double dt = 0.0);                   ///< Calculates the contact force between particles
    void UpdateParameters ();                           ///< Update the parameters

    // Data
    //ListContacts_t Lvv;                            ///< List of edge-edge contacts 
    //FrictionMap_t  Fdvv;      ///< Static friction displacement for pair of edges
    Vec3_t         Fdvv;                           ///< Static Friction displacement for the vertex vertex pair
    Vec3_t         Fdr;                            ///< Rolling displacement 
    double         beta;                           ///< Rolling stiffness coefficient
    double         eta;                            ///< Plastic moment coefficient
    double         Bn;                             ///< Elastic normal constant for the cohesion
    double         Bt;                             ///< Elastic tangential constant for the cohesion
    double         eps;                            ///< Maximun strain before fracture
    bool           cohesion;                       ///< A flag to determine if cohesion must be used for spheres or not

protected:
    void _update_rolling_resistance(double dt);               ///< Calculates the rolling resistance torque
};

class BInteracton: public Interacton // Interacton for cohesion
{
public:
    BInteracton (Particle * Pt1, Particle * Pt2, size_t Fi1, size_t Fi2); ///< Constructor requires pointers to both particles and alos the indixes of the size they share

    // Methods
    bool UpdateContacts (double alpha);    ///< Update contacts by verlet algorithm
    bool CalcForce      (double dt = 0.0); ///< Calculates the contact force between particles
    void UpdateParameters ();              ///< Update the parameters in case they change

    // Data
    size_t IF1;                            ///< Index of the shared face for particle 1
    size_t IF2;                            ///< Index of the shared face for particle 2
    double Area;                           ///< Area of the shared side
    double Bn;                             ///< Elastic normal constant for the cohesion
    double Bt;                             ///< Elastic tangential constant for the cohesion
    double Bm;                             ///< Elastic for the cohesion torque
    double Gn;                             ///< Dissipative normal constant
    double Gt;                             ///< Dissipative tangential constant
    double L0;                             ///< Equilibrium distance
    double An;                             ///< Angular displacement
    double eps;                            ///< Maximun strain before fracture
    double eta;                            ///< ratio between normal and tangential breaking threshold
    bool   valid;                          ///< Check if the bound has not been broken
    double s1,t1;                          ///< Planar coordinates for face F1
    double s2,t2;                          ///< Planar coordinates for face F2
    Vec3_t Fnet;                           ///< Net Normal force excerted by the interacton
    Vec3_t Ftnet;                          ///< Net Tangential force excerted by the interacton
    Vec3_t xnet;                           ///< Position where the force is applied

};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

// Collision interacton

inline CInteracton::CInteracton (Particle * Pt1, Particle * Pt2)
{
	//std::cout << "**************************1*******************************" << std::endl;
    First           = true;
    P1              = Pt1;
    P2              = Pt2;
    I1              = P1->Index;
    I2              = P2->Index;
    //double r1       = pow(P1->Props.V,1.0/3.0);
    //double r2       = pow(P2->Props.V,1.0/3.0);
    //Kn              = (r1+r2)*ReducedValue(Pt1->Props.Kn,Pt2->Props.Kn);
    //Kt              = (r1+r2)*ReducedValue(Pt1->Props.Kt,Pt2->Props.Kt);
    Kn              = ReducedValue(Pt1->Props.Kn,Pt2->Props.Kn);
    Kt              = ReducedValue(Pt1->Props.Kt,Pt2->Props.Kt);
    double me       = ReducedValue(Pt1->Props.m ,Pt2->Props.m );
    Gn              = 2*ReducedValue(Pt1->Props.Gn,Pt2->Props.Gn);
    Gt              = 2*ReducedValue(Pt1->Props.Gt,Pt2->Props.Gt);


/*
        std::cout << "Pt1->Props.x=                 " <<Pt1->x  <<  std::endl; //###################################2022.01.
        std::cout << "Pt2->Props.x=                 " <<Pt2->x  <<  std::endl; //###################################2022.01.
        std::cout << "Pt2->Props.Dmax=                 " <<Pt2->Dmax  <<  std::endl; //###################################2022.01.
        std::cout << "Pt2->Props.R=                 " <<Pt2->Props.R  <<  std::endl; //###################################2022.01.
        std::cout << "Pt2->Props.particle_radius=                 " <<Pt2->Props.particle_radius  <<  std::endl; //###################################2022.01.
        std::cout << "Pt1->Props.particle_radius=                 " <<Pt1->Props.particle_radius  <<  std::endl; //###################################2022.01.
	*/
           


	x_t        =0.0;///***
	thei_t     =0.0;///***
	thei_tw    =0.0;///***

	R	= ReducedValue(P1->Props.R,P2->Props.R);  ///***

//*****************************************irregular by qian, 2022.10.13******************************************************//
     double  P1_R,P2_R;
	 P1_R=pow(3*P1->Props.V/(4*M_PI),(1.0/3.0));
	 P2_R=pow(3*P2->Props.V/(4*M_PI),(1.0/3.0));
	 R	= ReducedValue(P1_R,P2_R);  ///***

//*****************************************irregular by qian, 2022.10.13******************************************************//

	
	E	= 1.0/((1-Pt1->Props.Sig*Pt1->Props.Sig)/Pt1->Props.E + (1-Pt2->Props.Sig*Pt2->Props.Sig)/Pt2->Props.E);///***
	G	= 1.0/((2-Pt1->Props.Sig)/Pt1->Props.G + (2-Pt2->Props.Sig)/Pt2->Props.G);///***
	Gama= ReducedValue(Pt1->Props.Gama,Pt2->Props.Gama); //??????????????????

    if (Gn < -0.001)
    {
        if (fabs(Gn)>1.0) throw new Fatal("CInteracton the restitution coefficient is greater than 1");
        Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
        Gt = 0.0;
    }
    Gn *= me;
    Gt *= me;
    //Mu              = 2*ReducedValue(Pt1->Props.Mu,Pt2->Props.Mu);
    if (Pt1->Props.Mu>1.0e-12&&Pt2->Props.Mu>1.0e-12)
    {
        if (!Pt1->IsFree())      Mu = Pt1->Props.Mu;
        else if (!Pt2->IsFree()) Mu = Pt2->Props.Mu;
        else                     Mu = std::max(Pt1->Props.Mu,Pt2->Props.Mu);
    }
    else 
    {
        Mu          = 0.0;
    }
    //if (Pt2->Tag<-1) std::cout << Mu << " " << Pt1->Tag << " " << Pt2->Tag << std::endl;

	  Nc = 0;   //############
	  Nsc = 0;//############
	
	  Epot = 0.0;//############
	  Fdr  = OrthoSys::O;//############
	  Fdvv = OrthoSys::O;//############

	
    CalcForce(0.0);

#ifdef USE_THREAD
    pthread_mutex_init(&lck,NULL);
#endif
}

inline bool CInteracton::UpdateContacts (double alpha)
{
	//std::cout << "**************************2*******************************" << std::endl;

    Lee.Resize(0);
    Lvf.Resize(0);
    Lfv.Resize(0);
    Lvt.Resize(0);
    Ltv.Resize(0);
    Lvc.Resize(0);
    Lcv.Resize(0);
	  Lvv.Resize(0);//***************
	 // Lve.Resize(0);//***************
	 // Lev.Resize(0);//***************
    if (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*alpha)
    {
        _update_contacts (P1->Edges     ,P2->Edges     ,Lee,alpha);
        _update_contacts (P1->Verts     ,P2->Faces     ,Lvf,alpha);
        _update_contacts (P1->Faces     ,P2->Verts     ,Lfv,alpha);
        _update_contacts (P1->Verts     ,P2->Tori      ,Lvt,alpha);
        _update_contacts (P1->Tori      ,P2->Verts     ,Ltv,alpha);
        _update_contacts (P1->Verts     ,P2->Cylinders ,Lvc,alpha);
        _update_contacts (P1->Cylinders ,P2->Verts     ,Lcv,alpha);
		_update_contacts (P1->Verts ,P2->Verts	   ,Lvv,alpha);//***************
	//	_update_contacts (P1->Verts ,P2->Edges 	   ,Lve,alpha);//***************2022.11.23delete
		//_update_contacts (P1->Edges ,P2->Verts 	   ,Lev,alpha);//***************2022.11.23delete//     ||Lve.Size()>0||Lev.Size()>0
        if (Lee.Size()>0||Lvf.Size()>0||Lfv.Size()>0||Ltv.Size()>0||Lvt.Size()>0||Lcv.Size()>0||Lvc.Size()>0||Lvv.Size()>0) return true;
        else return false;
    }
    else return false;
}

inline bool CInteracton::CalcForce (double dt)
{
	//std::cout << "**************************3*******************************" << std::endl;

    bool   overlap = false;
    Epot   = 0.0;
    dEvis  = 0.0;
    dEfric = 0.0;
    Nc     = 0;
    Nsc    = 0;
    Nr     = 0;
    Fnet   = OrthoSys::O;
    Ftnet  = OrthoSys::O;
    Xc     = OrthoSys::O;
    F1     = OrthoSys::O;
    F2     = OrthoSys::O;
    T1     = OrthoSys::O;
    T2     = OrthoSys::O;
    Branch = P1->x-P2->x;
    if (norm(Branch) > P1->Dmax + P2->Dmax) return false;
    if (_update_disp_calc_force (P1->Edges     ,P2->Edges     ,Fdee,Lee,dt,1)) overlap = true;
    if (_update_disp_calc_force (P1->Verts     ,P2->Faces     ,Fdvf,Lvf,dt,0)) overlap = true;
    if (_update_disp_calc_force (P1->Faces     ,P2->Verts     ,Fdfv,Lfv,dt,2)) overlap = true;
    if (_update_disp_calc_force (P1->Verts     ,P2->Tori      ,Fdvt,Lvt,dt,0)) overlap = true;
    if (_update_disp_calc_force (P1->Tori      ,P2->Verts     ,Fdtv,Ltv,dt,0)) overlap = true;
    if (_update_disp_calc_force (P1->Verts     ,P2->Cylinders ,Fdvc,Lvc,dt,0)) overlap = true;
    if (_update_disp_calc_force (P1->Cylinders ,P2->Verts     ,Fdcv,Lcv,dt,0)) overlap = true;
		if (_update_disp_calc_force (P1->Verts ,P2->Verts	  ,Fdv_v,Lvv,dt,4)) overlap = true;//**************
	 //   if (_update_disp_calc_force (P1->Verts ,P2->Edges	  ,Fdve,Lve,dt,0)) overlap = true;//**************2022.11.23delete
	  //  if (_update_disp_calc_force (P1->Edges ,P2->Verts	  ,Fdve,Lev,dt,0)) overlap = true;//**************2022.11.23delete

    return overlap;
}

inline void CInteracton::UpdateParameters ()
{
	//std::cout << "**************************4*******************************" << std::endl;

    Kn              = 2*ReducedValue(P1->Props.Kn,P2->Props.Kn);
    Kt              = 2*ReducedValue(P1->Props.Kt,P2->Props.Kt);
    double me       = ReducedValue(P1->Props.m ,P2->Props.m );
    Gn              = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn);
    Gt              = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt);
    if (Gn < -0.001)
    {
        if (fabs(Gn)>1.0) throw new Fatal("CInteracton the restitution coefficient is greater than 1");
        Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
        Gt = 0.0;
    }
    Gn *= me;
    Gt *= me;
    if (P1->Props.Mu>1.0e-12&&P2->Props.Mu>1.0e-12)
    {
        if (!P1->IsFree())      Mu = P1->Props.Mu;
        else if (!P2->IsFree()) Mu = P2->Props.Mu;
        else                    Mu = std::max(P1->Props.Mu,P2->Props.Mu);
    }
    else 
    {
        Mu          = 0.0;
    }
}

template<typename FeatureA_T, typename FeatureB_T>
inline bool CInteracton::_update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, ListContacts_t & L, double dt, double ee)
{
 //std::cout << "**************************5*******************************" << std::endl;
    // update
    std::cout << "L.Size()=**"<<L.Size() << std::endl;
    for (size_t k=0; k<L.Size(); ++k)
    {
        size_t i = L[k].first;
        size_t j = L[k].second;

			 double realpha=0.25*(P1->Props.R+P2->Props.R);
			 if(norm(P2->x-P1->x)>P2->Props.particle_radius+P1->Props.particle_radius+ realpha)      
						{
                          continue;   //#################################################2022.01.13
			           }
           // if (!Overlap ((*A[i]), (*B[j]),P1->Props.R + realpha,P2->Props.R + realpha)) continue;//#################################################2022.01.13
	    	if (Distance ((*A[i]), (*B[j]))>P1->Props.R+P2->Props.R+2*realpha) continue;//#################################################2022.01.13


        Vec3_t xi, xf;
        Distance ((*A[i]), (*B[j]), xi, xf);
		std::cout << "xi="<<xi<< " xf="<< xf<< " P1->x="<< P1->x<< " P2->x="<< P2->x << std::endl;	  
        double dist  = norm(xf-xi);
/*********************************************************************************************************************
in order to prohbit the fiber threading phenomenon by qian   2022.01.08
*********************************************************************************************************************/
//The wrong algorithm, qian 2022.10.13
/*                
      if(norm(xf-xi)<=P1->Props.R)      
			{
                           dist  = -norm(xf-xi);
			}

	     if(norm(xf-xi)<=P2->Props.R)       
				{
		                   dist  = -norm(xf-xi);
				}     */
/*********************************************************************************************************************

*********************************************************************************************************************/
                 double delta = P1->Props.R + P2->Props.R - dist;

                   // dist  = norm(xf-xi);  //******************************************************************************************************



		 double a0	= pow((9*M_PI*Gama*R*R/E),(1.0/3.0)); //equilibrium contact radius**************
		 double delta_c= a0*a0/(2.0*R*pow(6.0,1.0/3.0)); // critical overlap********************
						
         if (delta>-delta_c)
        {
        /*
#ifdef USE_CHECK_OVERLAP
            //if (delta > 0.4*std::min(P1->Props.R,P2->Props.R))
            if (delta > 0.4*(P1->Props.R+P2->Props.R))
            {
                std::cout << std::endl; 
                std::cout << "Maximun overlap between " << P1->Index         << " and " << P2->Index <<  std::endl; 
                std::cout << "Is the first collision? " << First             << std::endl; 
                std::cout << "Overlap                 " << delta             <<  std::endl; 
                std::cout << "Particle's tags         " << P1->Tag           << " and " << P2->Tag   <<  std::endl; 
                std::cout << "Memory address          " << P1                << " and " << P2        <<  std::endl; 
                std::cout << "Position particle 1     " << P1->x             << std::endl;
                std::cout << "Position particle 2     " << P2->x             << std::endl;
                std::cout << "Velocity particle 1     " << P1->v             << std::endl;
                std::cout << "Velocity particle 2     " << P2->v             << std::endl;
                std::cout << "Ang Velocity particle 1 " << P1->w             << std::endl;
                std::cout << "Ang Velocity particle 2 " << P2->w             << std::endl;
                std::cout << "Mass particle 1         " << P1->Props.m       << std::endl;
                std::cout << "Mass particle 2         " << P2->Props.m       << std::endl;
                std::cout << "Diameter particle 1     " << P1->Dmax          << std::endl;
                std::cout << "Diameter particle 2     " << P2->Dmax          << std::endl;
                std::cout << "Sradius particle 1      " << P1->Props.R       << std::endl;
                std::cout << "Sradius particle 2      " << P2->Props.R       << std::endl;
                std::cout << "Number of faces  1      " << P1->Faces.Size()  << std::endl;
                std::cout << "Number of faces  2      " << P2->Faces.Size()  << std::endl;
                //P1->Tag = 10000;
                //P2->Tag = 10000;
                return true;
                //throw new Fatal("Interacton::_update_disp_calc_force: Maximun overlap detected between particles %d(%d) and %d(%d)",P1->Index,P1->Tag,P2->Index,P2->Tag);
            }
#endif  */
            // Count a contact
            Nc++;
            //First = false;


			if(P1->Props.fixedparticle==true||P2->Props.fixedparticle==true)
								{
								P1->Props.ready_fixedparticle=true;
									P2->Props.ready_fixedparticle=true;
								}

			



            // update force
            Vec3_t n = (xf-xi)/dist;
std::ofstream outfile123;
		  outfile123.open("JKR_123.txt", std::ios::app);
		 outfile123<< delta<<"					   "<<dist<<"					   "<< xf<<" 					"<<xi<<"					"<<n<< std::endl;
			 //outfile << delta<<" 				"<< ratio<<"				 "<<Fratio<<"				  "<<Fne<< std::endl;
			 outfile123.close();

			
            Vec3_t x = xi+n*((P1->Props.R*P1->Props.R-P2->Props.R*P2->Props.R+dist*dist)/(2*dist));
            Xc += x;
            Vec3_t t1,t2,x1,x2;
            Rotation(P1->w,P1->Q,t1);
            Rotation(P2->w,P2->Q,t2);
            x1 = x - P1->x;
            x2 = x - P2->x;
          
			Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));    //
            Vec3_t vt = vrel - dot(n,vrel)*n;
			
//#############################################################################################################################
	
			if(!First)
		   {

			if(delta<=0.0){

			            //double a0	= pow((9*M_PI*Gama*R*R/E),(1.0/3.0)); //equilibrium contact radius
						//double delta_c= a0*a0/(2.0*R*pow(6.0,1.0/3.0)); // critical overlap
						double a	= pow(fabs(delta)*R,0.5); //contact radius	
			
				
						double e_alpha =0.5;
						double alpha=1.2728-4.2783*e_alpha+11.087*pow(e_alpha,2)-22.348*pow(e_alpha,3)+27.467*pow(e_alpha,4)-18.022*pow(e_alpha,5)+4.8218*pow(e_alpha,6);//#########################################################
						alpha=0.2;
						 Kt 	  = 8.0*G*a; //tangential stiffness coefficient
						double kn	  = 4.0*E*a/3.0;
						Gn	   = alpha*pow(P1->Props.m*kn,0.5); //normal friction/damping coefficient
						Gt	   = Gn;  //tangential friction/damping coefficient
			
						double Fc	= 3*M_PI*Gama*R; // critical force
						 double ratio = delta/delta_c;
						 double Fratio=0.0;
										   
										  
									
										   double Ftable[90][2]={ 
										   -0.99905224	,   -0.584512255 ,
										   -0.99004279	,  -0.648545190 ,
										   -0.98535915	, -0.667864607 ,
										   -0.98241445	,   -0.678344691, 
										   -0.97888911	,   -0.689710370 ,
										   -0.97051696	,   -0.713063228 ,
										   -0.95881144	,   -0.740078782 ,
										   -0.94238978	,   -0.771142447 ,
										   -0.93931370	,   -0.776323257 ,
										   -0.92550233	,   -0.797696855 ,
										   -0.89533079	,   -0.836322973 ,
										   -0.88886621	,   -0.843500608 ,
										   -0.86272505	,   -0.869503278 ,
										   -0.82199688	,   -0.902409047 ,
										   -0.80553361	,   -0.913632383 ,
										   -0.79985731	,   -0.917264033 ,
										   -0.78679635	,   -0.925190035 ,
										   -0.73329718	,  -0.952169085 ,
										   -0.71108845	,   -0.961131882 ,
										   -0.69893312	,   -0.965550173 ,
										   -0.68488742	,   -0.970252836 ,
										   -0.64840157	,   -0.980597232 ,
										   -0.62541686	,   -0.985835430 ,
										   -0.59928035	,   -0.990692257 ,
										   -0.58231967	,   -0.993255659 ,
										   -0.56957678	,   -0.994891655 ,
										   -0.54055165	,   -0.997733469 ,
										   -0.52054777	,   -0.999010223 ,
										   -0.49197014	,   -0.999922834 ,
										   -0.48644843	,   -0.999980168 ,
										   -0.47261207	,   -0.999959918 ,
										   -0.46528018	,   -0.999855832 ,
										   -0.45123374	,   -0.999479774 ,
										   -0.43947734	,   -0.998990174 ,
										   -0.42648406	,   -0.998267961 ,
										   -0.41054886	,   -0.997128678 ,
										   -0.39204589	,   -0.995464557 ,
										   -0.38560246	,   -0.994800868 ,
										   -0.36979780	,   -0.992992744 ,
										   -0.35617870	,   -0.991233050 ,
										   -0.34043559	,   -0.988971576 ,
										   -0.32476582	,   -0.986483816 ,
										   -0.29791107	,   -0.981685750 ,
										   -0.28051517	,   -0.978227038 ,
										   -0.26587058	,   -0.975107097 ,
										   -0.23648696	,   -0.968287637 ,
										   -0.20653917	,   -0.960591612 ,
										   -0.19559450	,   -0.957596350 ,
										   -0.17905038	,   -0.952887041 ,
										   -0.15811715	,   -0.946620575 ,
										   -0.12448193	,   -0.935849603 ,
										   -0.09940752	,   -0.927272127 ,
										   -0.09067325	,   -0.924176837 ,
										   -0.08126348	,   -0.920780891 ,
										   -0.07993753	,   -0.920297287 ,
										   -0.07873777	,   -0.919858628 ,
										   -0.07664266	,   -0.919090162 ,
										   -0.07213905	,   -0.917427760 ,
										   -0.06881591	,   -0.916191927 ,
										   -0.06636782	,   -0.915276546 ,
										   -0.06233643	,   -0.913759981 ,
										   -0.05963870	,   -0.912738775 ,
										   -0.05602120	,   -0.911361439 ,
										   -0.05329057	,   -0.910315749 ,
										   -0.05199593	,   -0.909818161 ,
										   -0.04969871	,   -0.908932381 ,
										   -0.04797466	,   -0.908265212 ,
										   -0.04635744	,   -0.907637519 ,
										   -0.04225428	,   -0.906036885 ,
										   -0.04008153	,   -0.905184614 ,
										   -0.03883383	,   -0.904693739 ,
										   -0.03556953	,   -0.903404445 ,
										   -0.03248330	,   -0.902178797 ,
										   -0.03019098	,   -0.901264244 ,
										   -0.02840943	,   -0.900551007 ,
										   -0.02531810	,   -0.899308294 ,
										   -0.01945889	,   -0.896935189 ,
										   -0.01116457	,   -0.893536360 ,
										   -0.00936052	,   -0.892791007 ,
										   -0.00887078	,   -0.892588293 ,
										   -0.00806020	,   -0.892252427 ,
										   -0.00728260	,   -0.891929813 ,
										   -0.00685426	,   -0.891751931 ,
										   -0.00606659	,   -0.891424503 ,
										   -0.00593068	,   -0.891367968 ,
										   -0.00506905	,   -0.891009245 ,
										   -0.00457149	,   -0.890801872 ,
										   -0.00326904	 ,-0.890258262 ,
										   -0.00249866	 , -0.889936195 ,
										   -0.00119656	 ,  -0.889390941 
										     };
												   for (size_t i=0; i<90; i++)
												   {
													   if (ratio > Ftable[i][0])
													   {
														   Fratio=Ftable[i][1]+(Ftable[i+1][1]-Ftable[i][1])/(Ftable[i+1][0]-Ftable[i][0])*(ratio-Ftable[i][0]);
														   //std::cout << " i=	" <<i<< " Ftable=  " << Ftable[i][1]<< " ratio=    " << ratio<< std::endl;
														   continue;
													   }
											   }
							  Vec3_t Fne = Fratio*Fc*n; //Modification 1   attention: the force direction
							   Vec3_t Fnd = Gn*dot(n,vrel)*n;
							  

		
							   Vec3_t F = Fne + Fnd ; 
							 F1	+= (-F);
				             F2	+=  (F);

		std::ofstream outfile0;
									   outfile0.open("JKR_n.txt", std::ios::app);
									   if(n(0)>0)
										   {
										   if(F(0)>0)
											   {
									   outfile0 << delta<<",	"<<dot(n,vrel)<<",	 "<<norm(F)<<std::endl;
									   //outfile << delta<<"				 "<< ratio<<"				  "<<Fratio<<"				   "<<Fne<< std::endl;
												}
										   else{
									   outfile0 << delta<<",	"<<dot(n,vrel)<<",	 "<<-norm(F)<<std::endl;
										   }
										   }
									   else{
												  if(F(0)>0)
																				   {
																		   outfile0 << delta<<",	"<<-dot(n,vrel)<<",   "<<norm(F)<<std::endl;
																		   //outfile << delta<<"				 "<< ratio<<"				  "<<Fratio<<"				   "<<Fne<< std::endl;
																					}
																			   else{
																		   outfile0 << delta<<",	"<<-dot(n,vrel)<<",   "<<-norm(F)<<std::endl;
																			   }
		
										   }
									   outfile0.close(); 
								
							 
        Vec3_t T, Tt;
		 
        Tt = cross (x1,F) ;  //
        Quaternion_t q;
        Conjugate (P1->Q,q);
        Rotation  (Tt,q,T);
         T1 += (-T);
		
        Tt = cross (x2,F) ;  // 
        Conjugate (P2->Q,q);
         Rotation  (Tt,q,T);
         T2 += T;
				}
				}
		if(delta>0.0){

		    First=false;
		
			double a	= pow(delta*R,0.5); //contact radius	
			
			 //******************************************************************************contact radius of edge-edge 2022.10.06 by qian, begin******************************//

				        std::ofstream outfilez;
						outfilez.open("JKR_a.txt", std::ios::app);
						 outfilez << delta<<"				   "<<a<< std::endl;
						//outfile << delta<<"				  "<< ratio<<"				   "<<Fratio<<" 				"<<Fne<< std::endl;
						outfilez.close();



			 
				  //******************************************************************************contact radius of edge-edge 2022.10.06 by qian, begin******************************//
				  //ee=0;
					  if( ee==1)
					  {
					  std::cout << " ee=  " <<ee<<	  std::endl;
				  
					  double AREA=0.0;
					  double AREA_distance_ratio=0.0;
				  
					  double distance_ratio=0.0;
					  double AREA_max=0.0;
				  
					  AREA_max=(2*P1->Props.R) * (2*P2->Props.R);
				  
					  distance_ratio=delta/(P1->Props.R + P2->Props.R);  //delta
					  std::cout << " delta=   " <<delta<< " distance_ratio=   " <<distance_ratio<<	std::endl;
				  
					  double d_table[31][2]={					  
									0 ,  0,
									0.033333333, 0.102090762,
									0.066666667, 0.198901358,
									0.1, 0.290383371,
									0.133333333, 0.376485668,
									0.166666667, 0.457148603,
									0.2, 0.532313289,
									0.233333333, 0.601760814,
									0.266666667, 0.665727973,
									0.3, 0.723887387,
									0.333333333, 0.776125638,
									0.366666667, 0.822490667,
									0.4, 0.862760641,
									0.433333333, 0.896853571,
									0.466666667, 0.924425674,
									0.5, 0.945490119,
									0.533333333, 0.959760927,
									0.566666667, 0.971723823,
									0.6, 0.980088574,
									0.633333333, 0.986379018,
									0.666666667, 0.990893703,
									0.7, 0.994147524,
									0.733333333, 0.996414766,
									0.766666667, 0.997930109,
									0.8, 0.99884411,
									0.833333333, 0.999457702,
									0.866666667, 0.9997651,
									0.9, 0.999927923,
									0.933333333, 0.999982074,
									0.966666667, 0.99999752,
									1,	 1,
									};
				  
									   for (size_t i=0; i<31; i++)
										  {
												if (distance_ratio > d_table[i][0])
												{
													AREA_distance_ratio=d_table[i][1]+(d_table[i+1][1]-d_table[i][1])/(d_table[i+1][0]-d_table[i][0])*(distance_ratio-d_table[i][0]);
													std::cout << " i=  " <<i<< " d_table[i][1]=  " << d_table[i][1]<< " d_table[i+1][1]=	" << d_table[i+1][1]<< " AREA_distance_ratio=	" << AREA_distance_ratio<< std::endl;
													continue;
												}
										  }
						   AREA=AREA_distance_ratio*AREA_max;
				  
				  
						std::cout << " AREA_distance_ratio=   " <<AREA_distance_ratio<< " AREA_max=   " <<AREA_max<< " AREA=  " <<AREA<<  std::endl;
				  
				  
						double angle_ratio=0.0;
						double AREA_angle_ratio=1.0;
						double cylinder_max=0.0;
				  
					   angle_ratio= Distance_angle ((*A[i]), (*B[j]), angle_ratio, cylinder_max); 
					   cylinder_max=Distance_cylinder_max ((*A[i]), (*B[j]), angle_ratio, cylinder_max); 
					   std::cout << " angle_ratio= " <<angle_ratio<<  std::endl;
					  std::cout << " cylinder_max= " <<cylinder_max<< std::endl;
				  
					  double A_table[22][2]={
						  1,   1,
						  0.984807753, 1.015426612,
						  0.939692621, 1.064177772,
						  0.866025404, 1.154700539,
						  0.766044443, 1.305407289,
						  0.64278761,  1.555723827,
						  0.5, 2,
						  0.422618262, 2.366201583,
						  0.342020143, 2.9238044,
						  0.309016994, 3.236067978,
						  0.275637356, 3.627955279,
						  0.258819045, 3.863703306,
						  0.224951054, 4.445411482,
						  0.207911691, 4.809734344,
						  0.190808995, 5.240843064,
						  0.173648178, 5.758770483,
						  0.156434465, 6.392453221,
						  0.139173101, 7.185296534,
						  0.121869343, 8.205509048,
						  0.104528463, 9.566772233,
						  0.087155743, 11.47371325,
						  0,   0,
						  };
				  
					  for (size_t i=0; i<22; i++)
											   {
												   if (angle_ratio < A_table[i][0])
												   {
													   AREA_angle_ratio=A_table[i][1]+(A_table[i+1][1]-A_table[i][1])/(A_table[i+1][0]-A_table[i][0])*(angle_ratio-A_table[i][0]);
													   //std::cout << " i=	" <<i<< " Ftable=  " << Ftable[i][1]<< " ratio=    " << ratio<< std::endl;
													   continue;
												   }
										   }
							std::cout << " angle_ratio=   " <<angle_ratio<< " AREA=   " <<AREA<< " AREA_angle_ratio=  " << AREA_angle_ratio<< std::endl;
				  
							AREA=AREA*AREA_angle_ratio;
							double cylinder_overlap;
							cylinder_overlap=cylinder_max*2*pow(fabs(delta)*R,0.5);
							
					  std::cout<< " delta=	" <<delta << " AREA=  " <<AREA<< " cylinder_max=  " << cylinder_max<< " cylinder_overlap=  " << cylinder_overlap<< std::endl;
					  
						if(AREA>cylinder_overlap)
						  {
							AREA=cylinder_overlap;
							}
				  
				  
						   a=pow(AREA/M_PI,0.5);
						   delta=a*a/R;

						   
				                      std::ofstream outfileq;
										 outfileq.open("JKR_b.txt", std::ios::app);
										  outfileq << delta<<"				   "<<a<< std::endl;
										 //outfile << delta<<"				   "<< ratio<<" 				"<<Fratio<<"				 "<<Fne<< std::endl;
										 outfileq.close();

						   
                 //*********************************************************************************jkr original formula*************************************************************************************//
				// double a0	= pow((9*M_PI*Gama*R*R/E),(1.0/3.0)); //equilibrium contact radius**************
						 
				// Fne=4*E*pow(a,3)/(3*R)-2*M_PI*a*a*pow(4*Gama*E/(M_PI*a),0.5)*n;
	   
				  
				  }
				  
				  
				  
				  //******************************************************************************contact radius of edge-edge 2022.10.06 by qian, end******************************//	  
				  //******************************************************************************contact radius of edge-edge 2022.10.06 by qian, end******************************//	
				  
			double e_alpha =0.8;
			double alpha=1.2728-4.2783*e_alpha+11.087*pow(e_alpha,2)-22.348*pow(e_alpha,3)+27.467*pow(e_alpha,4)-18.022*pow(e_alpha,5)+4.8218*pow(e_alpha,6);//#########################################################
			//alpha=0.5;
			 Kt 	  = 8.0*G*a; //tangential stiffness coefficient
			double kn	  = 4.0*E*a/3.0;
			Gn	   = alpha*pow(P1->Props.m*kn,0.5); //normal friction/damping coefficient
			Gt	   = Gn;  //tangential friction/damping coefficient

			 
			//Vec3_t Fc	= 3*M_PI*Gama*R*n; // critical force
			double Fc	= 3*M_PI*Gama*R; // critical force
			 double ratio = delta/delta_c;
		     double Fratio=0.0;
							   
							 
						
							   double Ftable[461][2]={ 
							   -1.0,      -0.95 ,
							   -0.7233934943 , -0.956312873 ,
							   -0.7120345372 , -0.960773949,
							   -0.6111813145 ,	-0.988620724,
							   -0.5622965558 ,	-0.995718097,
							   -0.5124456352 ,	-0.999375691,
							   -0.4785866453 ,-0.999997157,
							   -0.4208657037 ,	-0.997897821,
							   -0.3718303641 ,	-0.993239474,
							   -0.3209987346 ,-0.985851105,
							   -0.2862049717 ,-0.979388125,
							   -0.2362498370 ,-0.968229621,
							   -0.1863787347 ,-0.95499984,
							   -0.1244077761 ,-0.935824917,
							   -0.0854958065 ,-0.922316149,
							   -0.0440067886 ,-0.906721951,
							   -0.0410077963 ,-0.905548342,
							   -0.0351212136 ,-0.903226807,
							   -0.0311653133 ,-0.901653405,
							   -0.0255545350 ,-0.899403571,
							   -0.0222252979 ,-0.898058529,
							   -0.0194314907 ,-0.896924039,
							   -0.0168916310 ,-0.895888113,
							   -0.0125633717 ,-0.894112789,
							   -0.0094744550 ,-0.892838144,
							   -0.0064534441 ,-0.891585367,
							   -0.0033922364 ,-0.89030973,
							   -0.0030641204 ,-0.890172631,
							   -0.0026476899 ,-0.889998529,
							   0.000653961 ,   -0.888614101    ,
							   0.001070147 ,   -0.888439076    ,
							   0.001486415 ,   -0.888263903    ,
							   0.001902765 ,   -0.888088581    ,
							   0.002319197 ,   -0.88791311 ,
							   0.00273571  ,   -0.88773749 ,
							   0.003152305 ,   -0.887561721    ,
							   0.003568982 ,   -0.887385802    ,
							   0.00398574  ,   -0.887209735    ,
							   0.00440258  ,   -0.887033519    ,
							   0.004819502 ,   -0.886857154    ,
							   0.005236506 ,   -0.886680639    ,
							   0.005653591 ,   -0.886503976    ,
							   0.006070759 ,   -0.886327163    ,
							   0.006488008 ,   -0.886150201    ,
							   0.006905338 ,   -0.88597309 ,
							   0.007322751 ,   -0.885795829    ,
							   0.007740245 ,   -0.88561842 ,
							   0.008157821 ,   -0.885440861    ,
							   0.008575479 ,   -0.885263152    ,
							   0.008993218 ,   -0.885085295    ,
							   0.009411039 ,   -0.884907288    ,
							   0.009828942 ,   -0.884729131    ,
							   0.010246927 ,   -0.884550825    ,
							   0.010664993 ,   -0.88437237 ,
							   0.011083141 ,   -0.884193765    ,
							   0.011501371 ,   -0.884015011    ,
							   0.011919682 ,   -0.883836107    ,
							   0.012338076 ,   -0.883657054    ,
							   0.012756551 ,   -0.883477851    ,
							   0.013175107 ,   -0.883298498    ,
							   0.013593746 ,   -0.883118996    ,
							   0.014012466 ,   -0.882939344    ,
							   0.014431267 ,   -0.882759542    ,
							   0.014850151 ,   -0.882579591    ,
							   0.015269116 ,   -0.88239949 ,
							   0.015688163 ,   -0.882219239    ,
							   0.016107292 ,   -0.882038838    ,
							   0.016526502 ,   -0.881858288    ,
							   0.016945794 ,   -0.881677587    ,
							   0.017365168 ,   -0.881496737    ,
							   0.017784623 ,   -0.881315737    ,
							   0.01820416  ,   -0.881134587    ,
							   0.018623779 ,   -0.880953287    ,
							   0.01904348  ,   -0.880771837    ,
							   0.019463262 ,   -0.880590237    ,
							   0.019883126 ,   -0.880408487    ,
							   0.020303071 ,   -0.880226587    ,
							   0.020723098 ,   -0.880044536    ,
							   0.021143207 ,   -0.879862336    ,
							   0.021563398 ,   -0.879679986    ,
							   0.02198367  ,   -0.879497485    ,
							   0.022404024 ,   -0.879314834    ,
							   0.02282446  ,   -0.879132033    ,
							   0.023244977 ,   -0.878949082    ,
							   0.023665576 ,   -0.87876598 ,
							   0.024086257 ,   -0.878582729    ,
							   0.024507019 ,   -0.878399326    ,
							   0.024927863 ,   -0.878215774    ,
							   0.025348789 ,   -0.878032071    ,
							   0.025769796 ,   -0.877848218    ,
							   0.026190885 ,   -0.877664214    ,
							   0.026612056 ,   -0.87748006 ,
							   0.027033308 ,   -0.877295755    ,
							   0.027454642 ,   -0.8771113  ,
							   0.027876058 ,   -0.876926694    ,
							   0.028297555 ,   -0.876741938    ,
							   0.028719134 ,   -0.876557031    ,
							   0.029140794 ,   -0.876371974    ,
							   0.029562536 ,   -0.876186766    ,
							   0.02998436  ,   -0.876001407    ,
							   0.030828253 ,   -0.875630237    ,
							   0.031672472 ,   -0.875258465    ,
							   0.032517018 ,   -0.874886089    ,
							   0.03336189  ,   -0.87451311 ,
							   0.034207089 ,   -0.874139527    ,
							   0.035052614 ,   -0.873765341    ,
							   0.035898466 ,   -0.87339055 ,
							   0.036744644 ,   -0.873015155    ,
							   0.037591149 ,   -0.872639156    ,
							   0.03843798  ,   -0.872262552    ,
							   0.039285138 ,   -0.871885343    ,
							   0.040132622 ,   -0.87150753 ,
							   0.040980432 ,   -0.871129111    ,
							   0.041828569 ,   -0.870750086    ,
							   0.042677033 ,   -0.870370456    ,
							   0.043525822 ,   -0.869990219    ,
							   0.044374938 ,   -0.869609377    ,
							   0.04522438  ,   -0.869227929    ,
							   0.046074149 ,   -0.868845873    ,
							   0.046924244 ,   -0.868463211    ,
							   0.047774665 ,   -0.868079943    ,
							   0.048625413 ,   -0.867696067    ,
							   0.049476487 ,   -0.867311583    ,
							   0.050327887 ,   -0.866926492    ,
							   0.051179613 ,   -0.866540793    ,
							   0.052031666 ,   -0.866154487    ,
							   0.052884044 ,   -0.865767572    ,
							   0.053736749 ,   -0.865380048    ,
							   0.05458978  ,   -0.864991916    ,
							   0.055443138 ,   -0.864603175    ,
							   0.056296821 ,   -0.864213825    ,
							   0.057150831 ,   -0.863823866    ,
							   0.058005167 ,   -0.863433297    ,
							   0.058859828 ,   -0.863042119    ,
							   0.059714816 ,   -0.862650331    ,
							   0.061425771 ,   -0.861864924    ,
							   0.063138029 ,   -0.861077075    ,
							   0.064851592 ,   -0.860286782    ,
							   0.066566459 ,   -0.859494043    ,
							   0.068282629 ,   -0.858698858    ,
							   0.070000104 ,   -0.857901224    ,
							   0.071718882 ,   -0.857101139    ,
							   0.073438964 ,   -0.856298603    ,
							   0.075160349 ,   -0.855493612    ,
							   0.076883038 ,   -0.854686166    ,
							   0.07860703  ,   -0.853876263    ,
							   0.080332325 ,   -0.853063901    ,
							   0.082058923 ,   -0.852249079    ,
							   0.083786825 ,   -0.851431794    ,
							   0.085516029 ,   -0.850612045    ,
							   0.087246535 ,   -0.849789831    ,
							   0.088978345 ,   -0.84896515	   ,
							   0.090711457 ,   -0.848138	   ,
							   0.092445871 ,   -0.847308379    ,
							   0.094181588 ,   -0.846476287    ,
							   0.095918607 ,   -0.84564172	   ,
							   0.097656928 ,   -0.844804677    ,
							   0.099396551 ,   -0.843965158    ,
							   0.103751305 ,   -0.84185551	   ,
							   0.108114193 ,   -0.839730343    ,
							   0.112485215 ,   -0.83758963	   ,
							   0.11686437  ,   -0.835433346    ,
							   0.121251654 ,   -0.833261463    ,
							   0.125647067 ,   -0.831073957    ,
							   0.130050607 ,   -0.828870801    ,
							   0.134462272 ,   -0.826651968    ,
							   0.138882062 ,   -0.824417434    ,
							   0.143309973 ,   -0.822167171    ,
							   0.147746005 ,   -0.819901153    ,
							   0.152190156 ,   -0.817619355    ,
							   0.156642424 ,   -0.815321749    ,
							   0.179025471 ,   -0.803595709    ,
							   0.201611215 ,   -0.791470587    ,
							   0.224399463 ,   -0.778943124    ,
							   0.247390023 ,   -0.766010064    ,
							   0.270582706 ,   -0.752668154    ,
							   0.293977329 ,   -0.73891414	   ,
							   0.317573707 ,   -0.724744775    ,
							   0.341371662 ,   -0.710156811    ,
							   0.365371015 ,   -0.695147004    ,
							   0.389571593 ,   -0.679712109    ,
							   0.413973223 ,   -0.663848888    ,
							   0.438575736 ,   -0.647554102    ,
							   0.463378964 ,   -0.630824514    ,
							   0.488382744 ,   -0.613656889    ,
							   0.513586913 ,   -0.596047997    ,
							   0.53899131  ,   -0.577994606    ,
							   0.564595779 ,   -0.559493488    ,
							   0.590400164 ,   -0.540541417    ,
							   0.616404311 ,   -0.521135168    ,
							   0.642608069 ,   -0.501271519    ,
							   0.66901129  ,   -0.480947249    ,
							   0.695613826 ,   -0.460159139    ,
							   0.722415532 ,   -0.438903973    ,
							   0.749416266 ,   -0.417178534    ,
							   0.776615886 ,   -0.394979609    ,
							   0.804014254 ,   -0.372303987    ,
							   0.831611231 ,   -0.349148458    ,
							   0.859406683 ,   -0.325509812    ,
							   0.887400476 ,   -0.301384844    ,
							   0.915592478 ,   -0.276770348    ,
							   0.943982558 ,   -0.251663121    ,
							   0.972570589 ,   -0.22605996	   ,
							   1.001356444 ,   -0.199957667    ,
							   1.030339997 ,   -0.173353041    ,
							   1.059521126 ,   -0.146242886    ,
							   1.088899707 ,   -0.118624006    ,
							   1.118475621 ,   -0.090493207    ,
							   1.148248749 ,   -0.061847296    ,
							   1.178218973 ,   -0.032683082    ,
							   1.208386179 ,   -0.002997376    ,
							   1.211413729 ,   0			   ,
							   1.24179748  ,   0.030263031 ,
							   1.272377974 ,   0.061054249 ,
							   1.303155098 ,   0.092376839 ,
							   1.334128743 ,   0.124233985 ,
							   1.365298799 ,   0.15662887  ,
							   1.39666516  ,   0.189564675 ,
							   1.428227718 ,   0.22304458  ,
							   1.45998637  ,   0.257071765 ,
							   1.491941012 ,   0.291649405 ,
							   1.52409154  ,   0.326780678 ,
							   1.556437856 ,   0.362468759 ,
							   1.588979857 ,   0.39871682  ,
							   1.621717446 ,   0.435528035 ,
							   1.654650524 ,   0.472905575 ,
							   1.687778997 ,   0.510852609 ,
							   1.721102767 ,   0.549372307 ,
							   1.754621741 ,   0.588467836 ,
							   1.788335826 ,   0.628142362 ,
							   1.822244929 ,   0.668399052 ,
							   1.85634896  ,   0.709241068 ,
							   1.890647829 ,   0.750671575 ,
							   1.925141446 ,   0.792693734 ,
							   1.959829724 ,   0.835310706 ,
							   1.994712575 ,   0.878525651 ,
							   2.029789914 ,   0.922341727 ,
							   2.065061655 ,   0.966762093 ,
							   2.100527715 ,   1.011789904 ,
							   2.13618801  ,   1.057428317 ,
							   2.172042458 ,   1.103680486 ,
							   2.208090977 ,   1.150549564 ,
							   2.244333487 ,   1.198038705 ,
							   2.280769908 ,   1.246151059 ,
							   2.317400162 ,   1.294889777 ,
							   2.354224171 ,   1.344258009 ,
							   2.391241856 ,   1.394258904 ,
							   2.428453143 ,   1.444895608 ,
							   2.465857956 ,   1.496171269 ,
							   2.50345622  ,   1.548089033 ,
							   2.54124786  ,   1.600652045 ,
							   2.579232805 ,   1.653863448 ,
							   2.617410981 ,   1.707726386 ,
							   2.655782317 ,   1.762244    ,
							   2.694346742 ,   1.817419433 ,
							   2.733104186 ,   1.873255824 ,
							   2.77205458  ,   1.929756313 ,
							   2.811197854 ,   1.986924039 ,
							   2.85053394  ,   2.04476214  ,
							   2.890062772 ,   2.103273752 ,
							   2.929784283 ,   2.162462012 ,
							   2.969698407 ,   2.222330056 ,
							   3.009805077 ,   2.282881018 ,
							   3.050104231 ,   2.344118031 ,
							   3.090595803 ,   2.406044229 ,
							   3.13127973  ,   2.468662744 ,
							   3.17215595  ,   2.531976707 ,
							   3.213224399 ,   2.595989249 ,
							   3.254485018 ,   2.660703499 ,
							   3.295937744 ,   2.726122587 ,
							   3.337582517 ,   2.792249641 ,
							   3.379419278 ,   2.859087789 ,
							   3.421447967 ,   2.926640158 ,
							   3.463668526 ,   2.994909873 ,
							   3.506080896 ,   3.063900061 ,
							   3.54868502  ,   3.133613845 ,
							   3.591480841 ,   3.20405435  ,
							   3.634468303 ,   3.2752247   ,
							   3.677647349 ,   3.347128016 ,
							   3.721017924 ,   3.41976742  ,
							   3.764579974 ,   3.493146034 ,
							   3.808333443 ,   3.567266979 ,
							   3.852278279 ,   3.642133374 ,
							   3.896414427 ,   3.717748338 ,
							   3.940741834 ,   3.794114991 ,
							   3.985260449 ,   3.871236449 ,
							   4.02997022  ,   3.94911583  ,
							   4.074871094 ,   4.027756251 ,
							   4.119963021 ,   4.107160828 ,
							   4.16524595  ,   4.187332676 ,
							   4.210719832 ,   4.268274909 ,
							   4.256384617 ,   4.349990643 ,
							   4.302240254 ,   4.43248299  ,
							   4.348286697 ,   4.515755063 ,
							   4.394523896 ,   4.599809975 ,
							   4.440951803 ,   4.684650837 ,
							   4.487570371 ,   4.77028076  ,
							   4.534379554 ,   4.856702855 ,
							   4.581379303 ,   4.943920232 ,
							   4.628569574 ,   5.031936    ,
							   4.67595032  ,   5.120753268 ,
							   4.723521496 ,   5.210375144 ,
							   4.771283056 ,   5.300804736 ,
							   4.819234956 ,   5.392045151 ,
							   4.867377152 ,   5.484099496 ,
							   4.915709601 ,   5.576970876 ,
							   4.964232257 ,   5.670662397 ,
							   5.012945078 ,   5.765177164 ,
							   5.061848022 ,   5.860518281 ,
							   5.110941046 ,   5.956688853 ,
							   5.160224108 ,   6.053691982 ,
							   5.209697165 ,   6.151530772 ,
							   5.259360178 ,   6.250208324 ,
							   5.309213104 ,   6.349727741 ,
							   5.359255903 ,   6.450092123 ,
							   5.409488535 ,   6.551304572 ,
							   5.459910959 ,   6.653368189 ,
							   5.510523136 ,   6.756286071 ,
							   5.561325027 ,   6.86006132  ,
							   5.612316592 ,   6.964697034 ,
							   5.663497792 ,   7.070196311 ,
							   5.71486859  ,   7.176562249 ,
							   5.766428947 ,   7.283797946 ,
							   5.818178825 ,   7.391906498 ,
							   5.870118187 ,   7.500891002 ,
							   5.922246996 ,   7.610754554 ,
							   5.974565214 ,   7.721500249 ,
							   6.027072806 ,   7.833131183 ,
							   6.079769734 ,   7.94565045  ,
							   6.132655962 ,   8.059061145 ,
							   6.185731456 ,   8.17336636  ,
							   6.238996178 ,   8.28856919  ,
							   6.292450095 ,   8.404672727 ,
							   6.346093172 ,   8.521680064 ,
							   6.399925372 ,   8.639594292 ,
							   6.453946663 ,   8.758418504 ,
							   6.50815701  ,   8.87815579  ,
							   6.56255638  ,   8.998809242 ,
							   6.617144737 ,   9.120381948 ,
							   6.67192205  ,   9.242877    ,
							   6.726888285 ,   9.366297487 ,
							   6.78204341  ,   9.490646498 ,
							   6.83738739  ,   9.615927121 ,
							   6.892920196 ,   9.742142445 ,
							   6.948641793 ,   9.869295558 ,
							   7.004552151 ,   9.997389547 ,
							   7.060651237 ,   10.1264275  ,
							   7.116939021 ,   10.2564125  ,
							   7.173415471 ,   10.38734764 ,
							   7.230080556 ,   10.519236   ,
							   7.286934246 ,   10.65208067 ,
							   7.343976509 ,   10.78588473 ,
							   7.401207317 ,   10.92065127 ,
							   7.458626638 ,   11.05638337 ,
							   7.516234443 ,   11.19308411 ,
							   7.574030703 ,   11.33075659 ,
							   7.632015387 ,   11.46940388 ,
							   7.690188467 ,   11.60902906 ,
							   7.748549915 ,   11.74963523 ,
							   7.8070997   ,   11.89122545 ,
							   7.865837795 ,   12.03380282 ,
							   7.924764171 ,   12.17737041 ,
							   7.9838788   ,   12.32193131 ,
							   8.043181654 ,   12.46748859 ,
							   8.102672706 ,   12.61404534 ,
							   8.162351927 ,   12.76160464 ,
							   8.222219291 ,   12.91016956 ,
							   8.282274771 ,   13.05974319 ,
							   8.342518339 ,   13.2103286  ,
							   8.402949968 ,   13.36192888 ,
							   8.463569633 ,   13.5145471  ,
							   8.524377306 ,   13.66818634 ,
							   8.585372962 ,   13.82284968 ,
							   8.646556574 ,   13.97854019 ,
							   8.707928117 ,   14.13526096 ,
							   8.769487565 ,   14.29301505 ,
							   8.831234892 ,   14.45180556 ,
							   8.893170073 ,   14.61163554 ,
							   8.955293083 ,   14.77250809 ,
							   9.017603897 ,   14.93442626 ,
							   9.08010249  ,   15.09739315 ,
							   9.142788837 ,   15.26141182 ,
							   9.205662915 ,   15.42648535 ,
							   9.268724698 ,   15.59261682 ,
							   9.331974162 ,   15.75980929 ,
							   9.395411285 ,   15.92806584 ,
							   9.45903604  ,   16.09738955 ,
							   9.522848406 ,   16.26778348 ,
							   9.586848358 ,   16.43925072 ,
							   9.651035874 ,   16.61179432 ,
							   9.715410929 ,   16.78541738 ,
							   9.779973502 ,   16.96012295 ,
							   9.844723568 ,   17.13591411 ,
							   9.909661105 ,   17.31279393 ,
							   9.974786092 ,   17.49076548 ,
							   10.0400985  ,   17.66983184 ,
							   10.70352623 ,   19.52138495 ,
							   11.38567292 ,   21.48579431 ,
							   12.08651837 ,   23.5661265  ,
							   12.80604356 ,   25.76544563 ,
							   13.54423061 ,   28.08681356 ,
							   14.30106265 ,   30.53328999 ,
							   15.07652377 ,   33.10793259 ,
							   15.87059892 ,   35.81379712 ,
							   16.68327386 ,   38.65393751 ,
							   17.51453509 ,   41.63140599 ,
							   18.36436981 ,   44.74925316 ,
							   19.23276587 ,   48.01052807 ,
							   20.11971172 ,   51.41827831 ,
							   21.02519635 ,   54.97555007 ,
							   21.9492093  ,   58.68538821 ,
							   22.89174058 ,   62.55083632 ,
							   23.85278066 ,   66.57493679 ,
							   24.83232044 ,   70.76073086 ,
							   25.83035121 ,   75.11125865 ,
							   26.84686467 ,   79.62955925 ,
							   27.88185284 ,   84.31867073 ,
							   28.93530809 ,   89.18163019 ,
							   30.0072231  ,   94.22147381 ,
							   31.09759085 ,   99.44123688 ,
							   32.2064046  ,   104.8439538 ,
							   33.33365788 ,   110.4326583 ,
							   34.47934444 ,   116.2103832 ,
							   35.64345832 ,   122.1801604 ,
							   36.82599374 ,   128.3450215 ,
							   38.02694515 ,   134.7079969 ,
							   39.24630721 ,   141.2721168 ,
							   40.48407475 ,   148.0404104 ,
							   41.74024281 ,   155.0159065 ,
							   43.01480657 ,   162.2016331 ,
							   44.30776141 ,   169.6006177 ,
							   45.61910285 ,   177.2158873 ,
							   46.94882656 ,   185.0504683 ,
							   48.29692835 ,   193.1073866 ,
							   49.66340418 ,   201.3896675 ,
							   52.45146239 ,   218.6424162 ,
							   55.31297131 ,   236.8329076 ,
							   58.24790288 ,   255.9853281 ,
							   61.25623071 ,   276.123857  ,
							   64.33792995 ,   297.2726674 ,
							   67.49297716 ,   319.4559264 ,
							   70.72135021 ,   342.6977954 ,
							   74.02302817 ,   367.0224306 ,
							   77.39799123 ,   392.453983  ,
							   80.84622061 ,   419.016599  ,
							   84.36769848 ,   446.7344203 ,
							   87.96240792 ,   475.6315845 ,
							   91.63033283 ,   505.7322249 ,
							   95.37145788 ,   537.060471  ,
							   99.18576848 ,   569.6404485 ,
							   103.0732507 ,   603.4962799 , 

							   125.0000639 ,   806.2004    ,							   
							   150.0000044 ,   1059.989    ,
							   174.9998873 ,   1335.913    ,
                               199.9998    ,   1632.321    ,
                               229.9999886 ,   2013.199    ,
                               255.9999928 ,   2364.157    ,
                               279.9999972 ,   2704.387    ,
                               300.9999992 ,   3014.343    ,
                               399.9999606 ,   4618.133    ,
                               499.9991627 ,   6454.288    ,


							   };
									   for (size_t i=0; i<461; i++)
									   {
										   if (ratio > Ftable[i][0])
										   {
											   Fratio=Ftable[i][1]+(Ftable[i+1][1]-Ftable[i][1])/(Ftable[i+1][0]-Ftable[i][0])*(ratio-Ftable[i][0]);
											   //std::cout << " i=	" <<i<< " Ftable=  " << Ftable[i][1]<< " ratio=    " << ratio<< std::endl;
											   continue;
										   }
								   }
				  Vec3_t Fne = Fratio*Fc*n; //Modification 1	  attention: the force direction
				  Vec3_t Fnd = Gn*dot(n,vrel)*n;
				  Vec3_t Fn = Fne+Fnd;
				
		  std::ofstream outfile;                    
					  outfile.open("JKR_0.txt", std::ios::app);
					   outfile<<ee<<"					   "<< L.Size()<<"					   "<< k<<"					   "<< delta<<"					   "<<a<<"					   "<< ratio<<" 					"<<Fratio<<" 					"<<Fc<<" 					"<<n<<"					 "<<norm(Fne)<<"					   "<< Gn<<" 					"<< dot(n,vrel)<<" 					"<< Fnd<<" 					"<< norm(P2->v)<<"					   "<< norm(P1->v)<<"					   "<< norm(P2->w)<<"					   "<< norm(P1->w)<<"					   "<< std::endl;
					  //outfile << delta<<" 				"<< ratio<<"				 "<<Fratio<<"				  "<<Fne<< std::endl;
					  outfile.close();
			
			
			
			
			std::ofstream outfile0;
										   outfile0.open("JKR_n.txt", std::ios::app);
										   if(n(0)>0)
											   {
											   if(Fn(0)>0)
												   {
										   outfile0 << delta<<",	"<<dot(n,vrel)<<",	 "<<norm(Fn)<<std::endl;
										   //outfile << delta<<"				 "<< ratio<<"				  "<<Fratio<<"				   "<<Fne<< std::endl;
													}
											   else{
										   outfile0 << delta<<",	"<<dot(n,vrel)<<",	 "<<-norm(Fn)<<std::endl;
											   }
											   }
										   else{
													  if(Fn(0)>0)
																					   {
																			   outfile0 << delta<<",	"<<-dot(n,vrel)<<",   "<<norm(Fn)<<std::endl;
																			   //outfile << delta<<"				 "<< ratio<<"				  "<<Fratio<<"				   "<<Fne<< std::endl;
																						}
																				   else{
																			   outfile0 << delta<<",	"<<-dot(n,vrel)<<",   "<<-norm(Fn)<<std::endl;
																				   }
			
											   }
										   outfile0.close(); 


			
	        Vec3_t ts = vt/norm(vt);
	        x_t     +=(vt*dt);
			x_t     -= dot(x_t,n)*n;
			Vec3_t Ft =Kt* x_t+Gt*vt;
	        if (norm(Ft)>Mu*(norm(Fn)+2*Fc))  //Mu=0.3
	        {
	            
	            Ft = Mu*(norm(Fn)+2*Fc)*ts;
	        }       

		         Vec3_t F = Fn + Ft;
				 F1	+= (-F);
				 F2	+=  (F);

                   std::ofstream outfileF;
					 outfileF.open("JKR_F.txt", std::ios::app);
					  outfileF<<ee<<"					   "<< L.Size()<<"					   "<< k<<"					   "<< norm(Fn)<<" 					  "<<norm(Ft)<<" 				   "<<norm(F)<<"					"<<norm(F1)<<"					"<< std::endl;
					 //outfile << delta<<"				   "<< ratio<<" 				"<<Fratio<<"				 "<<Fne<< std::endl;
					 outfileF.close();


				 
/*
		   SFr 	 += dt*vt;
		   SFr		-= dot(SFr,n)*n;
		   Vec3_t tan= SFr;
		   if(norm(tan)>0.0) tan/=norm(tan);
		   if(norm(SFr)>Mu*(norm(Fn)+2*Fc)/Kt)
		   {
			   SFr = Mu*(norm(Fn)+2*Fc)/Kt*tan;
		   }

         Vec3_t F = Fne + Fnd + Kt*SFr + Gt*vt;
        

 */


if(Util::IsNan(norm(F))) 
								   {
							
									std::cout << "Fne = " << Fne		<< std::endl;
							
									std::cout << "Fnd = " << Fnd		<< std::endl;
									
									std::cout << "Ft = " << Ft		<< std::endl;
													  
								   }


		
		 
		 //Rolling resistance
		 double miu_R=2/3*Mu;
		 double yita_R=miu_R*norm(Fn);
		 double KR = 4*Fc*pow((a/a0),1.5);
		 
		 Vec3_t Vr = P1->Props.R*P2->Props.R*cross(Vec3_t(t1 - t2),n)/(P1->Props.R+P2->Props.R);
		 Fdr += Vr*dt;
		 Fdr -= dot(Fdr,n)*n;
		 Vec3_t tanr = Fdr;
		 if (norm(tanr)>0.0) tanr/=norm(tanr);
		 Vec3_t Mr = KR*Fdr+yita_R*Vr;
	
		  if (norm(Mr)>KR*0.1*2*M_PI*R)  
		 {
			 Mr = KR*0.1*2*M_PI*R*tanr;
			
		 }
		 
   

		 

		 /*******************************      twisting torque           *******************************/
			//Vec3_t omiga_n = dot(n,(P1->w-P2->w))*n;
			
        Vec3_t omiga_n = dot(n,(t1-t2))*n;
		 thei_tw+=(omiga_n*dt);
        Vec3_t Mt = ((Kt*a*a)/2)*(thei_tw)  +  ((Gn*a*a)/2)*omiga_n;
		 if (norm(Mt)>2/3*a*Mu*(norm(Fne)+2*Fc))	//Mu=0.3
				{
					Mt = 2/3*a*Mu*(norm(Fne)+2*Fc)*omiga_n/norm(omiga_n);
				}
				/*****************************************************/

		                               





		  /*******************************      rolling torque           *******************************/
/*
		 double miu_R=2/3*Mu;
		 double yita_R=miu_R*norm(Fne);
		 
         double KR = 4*Fc*pow((a/a0),1.5);

		Vec3_t omiga = t1-t2;
		Vec3_t Vl = -R*cross(omiga,n);	//??????????//negative ???-0.5*((P2->Props.R-P1->Props.R)/(P2->Props.R+P1->Props.R))*vt
		Vec3_t tR =  Vl/norm(Vl);
		 
		  thei_t+=(Vl*dt);
		 Vec3_t Mr = KR*thei_t;       //+  yita_R*Vl
		 
		 //19–64 milliradians=(19–64)*0.0573
		 //double angle = 60*0.0573*M_PI/180.0;
        if (norm(thei_t)>0.01*2*M_PI*P1->Props.R)	//Mu=0.3      20*10-3rad
				{
		    	 Mr = (KR*0.01*2*M_PI*P1->Props.R)*tR;		
				}
				
				/*****************************************************/


        Vec3_t T, Tt;
		 
        Tt = cross (x1,F) + Mt + P1->Props.R*cross(n,Mr); //+ cross(Mr,n)
        Quaternion_t q;
        Conjugate (P1->Q,q);
        Rotation  (Tt,q,T);
        T1 += (-T);
		
        Tt = cross (x2,F) + Mt+ P2->Props.R*cross(n,Mr); // +cross(Mr,n)
        Conjugate (P2->Q,q);
        Rotation  (Tt,q,T);
        T2 += T;

			}//##########################
        }
    }
    return false;
}

template<typename FeatureA_T, typename FeatureB_T>
inline void CInteracton::_update_contacts (FeatureA_T & A, FeatureB_T & B, ListContacts_t & L, double alpha)
{

	
		



	if(P1->Props.fixedparticle==true)
	if(P2->Props.fixedparticle==true)  return;

    for (size_t i=0; i<A.Size(); ++i)
    for (size_t j=0; j<B.Size(); ++j)
    {



		double realpha=0.25*(P1->Props.R+P2->Props.R);
        if (!Overlap ((*A[i]), (*B[j]),P1->Props.R + realpha,P2->Props.R + realpha)) continue;//#################################################2022.01.13
        if (Distance ((*A[i]), (*B[j]))<=P1->Props.R+P2->Props.R+2*realpha)
        {
            std::pair<int,int> p;
            p = std::make_pair(i,j);
            L.Push(p);
			
			P1->Props.interact=true;//2022.03.26  by qian
		    P2->Props.interact=true;//2022.03.26  by qian

			//2022.04.01 in order to make particle fixed with fixed particle
				if(P1->Props.fixedparticle==true||P2->Props.fixedparticle==true)
					{
					P1->Props.ready_fixedparticle=true;
						P2->Props.ready_fixedparticle=true;
					}

			
		
        }
		else{
			                      P1->Props.interact=false;//2022.03.26	by qian
								  P2->Props.interact=false;//2022.03.26	by qian

								  P1->Props.ready_fixedparticle=false; //2022.04.01
			                      P2->Props.ready_fixedparticle=false; //2022.04.01
			}
		
    }




	
}

//Collision interacton for spheres

inline CInteractonSphere::CInteractonSphere (Particle * Pt1, Particle * Pt2)
{
	//std::cout << "**************************6*******************************" << std::endl;

    First           = true;
    P1   = Pt1;
    P2   = Pt2;
    I1 = P1->Index;
    I2 = P2->Index;
    Kn   = 2*ReducedValue(Pt1->Props.Kn,Pt2->Props.Kn);
    Kt   = 2*ReducedValue(Pt1->Props.Kt,Pt2->Props.Kt);
    double me       = ReducedValue(Pt1->Props.m ,Pt2->Props.m );
    Gn              = 2*ReducedValue(Pt1->Props.Gn,Pt2->Props.Gn);
    Gt              = 2*ReducedValue(Pt1->Props.Gt,Pt2->Props.Gt);


	x_t 	   =0.0;///***
	thei_t	   =0.0;///***
	thei_tw    =0.0;///***
		

	R	= ReducedValue(P1->Props.R,P2->Props.R);  ///***
	E	= 1.0/((1-Pt1->Props.Sig*Pt1->Props.Sig)/Pt1->Props.E + (1-Pt2->Props.Sig*Pt2->Props.Sig)/Pt2->Props.E);///***
	G	= 1.0/((2-Pt1->Props.Sig)/Pt1->Props.G + (2-Pt2->Props.Sig)/Pt2->Props.G);///***
	Gama= ReducedValue(Pt1->Props.Gama,Pt1->Props.Gama); //??????????????????

	
    if (Gn < 0.0)
    {
        if (fabs(Gn)>1.0) throw new Fatal("CInteractonSphere the restitution coefficient is greater than 1");
        Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
        Gt = 0.0;
    }
    Gn *= me;
    Gt *= me;
    //Mu   = 2*ReducedValue(Pt1->Props.Mu,Pt2->Props.Mu);
    if (Pt1->Props.Mu>1.0e-12&&Pt2->Props.Mu>1.0e-12)
    {
        if (!Pt1->IsFree())      Mu = Pt1->Props.Mu;
        else if (!Pt2->IsFree()) Mu = Pt2->Props.Mu;
        else                     Mu = std::max(Pt1->Props.Mu,Pt2->Props.Mu);
    }
    else 
    {
        Mu          = 0.0;
    }
    beta = 2*ReducedValue(Pt1->Props.Beta,Pt2->Props.Beta);
    eta  = 2*ReducedValue(Pt1->Props.Eta,Pt2->Props.Eta);
    Bn   = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn);
    Bt   = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt);
    eps  = 2*ReducedValue(P1->Props.eps,P2->Props.eps);

    if (eps<0.0)
    {
        cohesion = true;
        if (Bn<Kn) throw new Fatal("CInteractonSphere::CInteractonSphere if cohesion is applied, Bn must be larger than Kn");
    }
    else
    {
        cohesion = false;
    }
    
    Nc = 0;
    Nsc = 0;

    Epot = 0.0;
    Fdr  = OrthoSys::O;
    Fdvv = OrthoSys::O;

    //std::pair<int,int> p;
    //p = std::make_pair(0,0);
    //Lvv.Push(p);

    CalcForce(0.0);
#ifdef USE_THREAD
    pthread_mutex_init(&lck,NULL);
#endif
}

inline void CInteractonSphere::_update_rolling_resistance(double dt)
{
	//std::cout << "**************************7*******************************" << std::endl;

    Vec3_t t1,t2;
    Rotation(P1->w,P1->Q,t1);
    Rotation(P2->w,P2->Q,t2);
    Vec3_t Normal = Fn/norm(Fn);
    Vec3_t Vr = P1->Props.R*P2->Props.R*cross(Vec3_t(t1 - t2),Normal)/(P1->Props.R+P2->Props.R);
    Fdr += Vr*dt;
    Fdr -= dot(Fdr,Normal)*Normal;
    Vec3_t tan = Fdr;
    if (norm(tan)>0.0) tan/=norm(tan);
    double Kr = beta*Kt;
    if (norm(Fdr)>eta*Mu*norm(Fn)/Kr)
    {
        Fdr = eta*Mu*norm(Fn)/Kr*tan;
        Nr++;
    }
    Vec3_t Ft = -Kr*Fdr;

    Vec3_t Tt = P1->Props.R*cross(Normal,Ft);
    Vec3_t T;
#ifdef USE_THREAD
    //pthread_mutex_lock(&lck);
    //pthread_mutex_lock(&P1->lck);
    //pthread_mutex_lock(&P2->lck);
    //std::lock_guard<std::mutex> lk1(P1->mtex);
    //std::lock_guard<std::mutex> lk2(P2->mtex);
#endif
    Quaternion_t q;
    Conjugate (P1->Q,q);
    Rotation  (Tt,q,T);
    //P1->T += T;
    T1 += T;

    Tt = P2->Props.R*cross(Normal,Ft);
    Conjugate (P2->Q,q);
    Rotation  (Tt,q,T);
    //P2->T -= T;
    T2 -= T;
#ifdef USE_THREAD
    //pthread_mutex_unlock(&lck);
    //pthread_mutex_unlock(&P1->lck);
    //pthread_mutex_unlock(&P2->lck);
#endif
}

inline bool CInteractonSphere::CalcForce(double dt)
{
	//std::cout << "**************************8*******************************" << std::endl;

    Epot   = 0.0;
    dEvis  = 0.0;
    dEfric = 0.0;
    Nc     = 0;
    Nsc    = 0;
    Nr     = 0;
    Fnet   = OrthoSys::O;
    Ftnet  = OrthoSys::O;
    Xc     = OrthoSys::O;
    F1     = OrthoSys::O;
    F2     = OrthoSys::O;
    T1     = OrthoSys::O;
    T2     = OrthoSys::O;
    Branch = P1->x-P2->x;

	if(norm(P2->x-P1->x)>P2->Props.particle_radius+P1->Props.particle_radius)   return false;
    Vec3_t xi = P1->x;
    Vec3_t xf = P2->x;
    double dist = norm(Branch);
/*********************************************************************************************************************
in order to prohbit the fiber threading phenomenon by qian   2022.01.08
*********************************************************************************************************************/
        if(norm(P2->x-P1->x)<=P2->Props.particle_radius)      
			{
                           dist  = -norm(xf-xi);
			}

if(norm(P2->x-P1->x)<=P1->Props.particle_radius)	  
			   {
						  dist	= -norm(xf-xi);
			   }

/*********************************************************************************************************************

*********************************************************************************************************************/
    double delta = P1->Props.R + P2->Props.R - dist;


     dist  = norm(xf-xi);


     double a0	= pow((9*M_PI*Gama*R*R/E),(1.0/3.0)); //equilibrium contact radius**************
		 double delta_c= a0*a0/(2.0*R*pow(6.0,1.0/3.0)); // critical overlap********************
						
        if (delta>-delta_c)
        {
        /*
#ifdef USE_CHECK_OVERLAP
            //if (delta > 0.4*std::min(P1->Props.R,P2->Props.R))
            if (delta > 0.4*(P1->Props.R+P2->Props.R))
            {
                std::cout << std::endl; 
                std::cout << "Maximun overlap between " << P1->Index         << " and " << P2->Index <<  std::endl; 
                std::cout << "Is the first collision? " << First             << std::endl; 
                std::cout << "Overlap                 " << delta             <<  std::endl; 
                std::cout << "Particle's tags         " << P1->Tag           << " and " << P2->Tag   <<  std::endl; 
                std::cout << "Memory address          " << P1                << " and " << P2        <<  std::endl; 
                std::cout << "Position particle 1     " << P1->x             << std::endl;
                std::cout << "Position particle 2     " << P2->x             << std::endl;
                std::cout << "Velocity particle 1     " << P1->v             << std::endl;
                std::cout << "Velocity particle 2     " << P2->v             << std::endl;
                std::cout << "Ang Velocity particle 1 " << P1->w             << std::endl;
                std::cout << "Ang Velocity particle 2 " << P2->w             << std::endl;
                std::cout << "Mass particle 1         " << P1->Props.m       << std::endl;
                std::cout << "Mass particle 2         " << P2->Props.m       << std::endl;
                std::cout << "Diameter particle 1     " << P1->Dmax          << std::endl;
                std::cout << "Diameter particle 2     " << P2->Dmax          << std::endl;
                std::cout << "Sradius particle 1      " << P1->Props.R       << std::endl;
                std::cout << "Sradius particle 2      " << P2->Props.R       << std::endl;
                std::cout << "Number of faces  1      " << P1->Faces.Size()  << std::endl;
                std::cout << "Number of faces  2      " << P2->Faces.Size()  << std::endl;
                //P1->Tag = 10000;
                //P2->Tag = 10000;
                return true;
                //throw new Fatal("Interacton::_update_disp_calc_force: Maximun overlap detected between particles %d(%d) and %d(%d)",P1->Index,P1->Tag,P2->Index,P2->Tag);
            }
#endif   */
            // Count a contact
            Nc++;
            //First = false;

            // update force
            Vec3_t n = (xf-xi)/dist;
			
			//std::cout << " n=	  " << n<< std::endl;
            Vec3_t x = xi+n*((P1->Props.R*P1->Props.R-P2->Props.R*P2->Props.R+dist*dist)/(2*dist));
            Xc += x;
            Vec3_t t1,t2,x1,x2;
            Rotation(P1->w,P1->Q,t1);
            Rotation(P2->w,P2->Q,t2);
            x1 = x - P1->x;
            x2 = x - P2->x;
            Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
            Vec3_t vt = vrel - dot(n,vrel)*n;

//#############################################################################################################################
	
			if(!First)
		   {

			if(delta<=0.0){

						double a	= pow(fabs(delta)*R,0.5); //contact radius	
						double e_alpha =0.8;
						double alpha=1.2728-4.2783*e_alpha+11.087*pow(e_alpha,2)-22.348*pow(e_alpha,3)+27.467*pow(e_alpha,4)-18.022*pow(e_alpha,5)+4.8218*pow(e_alpha,6);//#########################################################
						
						 Kt 	  = 8.0*G*a; //tangential stiffness coefficient
						double kn	  = 4.0*E*a/3.0;
						Gn	   = alpha*pow(P1->Props.m*kn,0.5); //normal friction/damping coefficient
						Gt	   = Gn;  //tangential friction/damping coefficient
			
						 
						//Vec3_t Fc	= 3*M_PI*Gama*R*n; // critical force
						double Fc	= 3*M_PI*Gama*R; // critical force
						 double ratio = delta/delta_c;
						 double Fratio=0.0;
										   
										  
									
										   double Ftable[90][2]={ 
										   -0.99905224	,   -0.584512255 ,
										   -0.99004279	,  -0.648545190 ,
										   -0.98535915	, -0.667864607 ,
										   -0.98241445	,   -0.678344691, 
										   -0.97888911	,   -0.689710370 ,
										   -0.97051696	,   -0.713063228 ,
										   -0.95881144	,   -0.740078782 ,
										   -0.94238978	,   -0.771142447 ,
										   -0.93931370	,   -0.776323257 ,
										   -0.92550233	,   -0.797696855 ,
										   -0.89533079	,   -0.836322973 ,
										   -0.88886621	,   -0.843500608 ,
										   -0.86272505	,   -0.869503278 ,
										   -0.82199688	,   -0.902409047 ,
										   -0.80553361	,   -0.913632383 ,
										   -0.79985731	,   -0.917264033 ,
										   -0.78679635	,   -0.925190035 ,
										   -0.73329718	,  -0.952169085 ,
										   -0.71108845	,   -0.961131882 ,
										   -0.69893312	,   -0.965550173 ,
										   -0.68488742	,   -0.970252836 ,
										   -0.64840157	,   -0.980597232 ,
										   -0.62541686	,   -0.985835430 ,
										   -0.59928035	,   -0.990692257 ,
										   -0.58231967	,   -0.993255659 ,
										   -0.56957678	,   -0.994891655 ,
										   -0.54055165	,   -0.997733469 ,
										   -0.52054777	,   -0.999010223 ,
										   -0.49197014	,   -0.999922834 ,
										   -0.48644843	,   -0.999980168 ,
										   -0.47261207	,   -0.999959918 ,
										   -0.46528018	,   -0.999855832 ,
										   -0.45123374	,   -0.999479774 ,
										   -0.43947734	,   -0.998990174 ,
										   -0.42648406	,   -0.998267961 ,
										   -0.41054886	,   -0.997128678 ,
										   -0.39204589	,   -0.995464557 ,
										   -0.38560246	,   -0.994800868 ,
										   -0.36979780	,   -0.992992744 ,
										   -0.35617870	,   -0.991233050 ,
										   -0.34043559	,   -0.988971576 ,
										   -0.32476582	,   -0.986483816 ,
										   -0.29791107	,   -0.981685750 ,
										   -0.28051517	,   -0.978227038 ,
										   -0.26587058	,   -0.975107097 ,
										   -0.23648696	,   -0.968287637 ,
										   -0.20653917	,   -0.960591612 ,
										   -0.19559450	,   -0.957596350 ,
										   -0.17905038	,   -0.952887041 ,
										   -0.15811715	,   -0.946620575 ,
										   -0.12448193	,   -0.935849603 ,
										   -0.09940752	,   -0.927272127 ,
										   -0.09067325	,   -0.924176837 ,
										   -0.08126348	,   -0.920780891 ,
										   -0.07993753	,   -0.920297287 ,
										   -0.07873777	,   -0.919858628 ,
										   -0.07664266	,   -0.919090162 ,
										   -0.07213905	,   -0.917427760 ,
										   -0.06881591	,   -0.916191927 ,
										   -0.06636782	,   -0.915276546 ,
										   -0.06233643	,   -0.913759981 ,
										   -0.05963870	,   -0.912738775 ,
										   -0.05602120	,   -0.911361439 ,
										   -0.05329057	,   -0.910315749 ,
										   -0.05199593	,   -0.909818161 ,
										   -0.04969871	,   -0.908932381 ,
										   -0.04797466	,   -0.908265212 ,
										   -0.04635744	,   -0.907637519 ,
										   -0.04225428	,   -0.906036885 ,
										   -0.04008153	,   -0.905184614 ,
										   -0.03883383	,   -0.904693739 ,
										   -0.03556953	,   -0.903404445 ,
										   -0.03248330	,   -0.902178797 ,
										   -0.03019098	,   -0.901264244 ,
										   -0.02840943	,   -0.900551007 ,
										   -0.02531810	,   -0.899308294 ,
										   -0.01945889	,   -0.896935189 ,
										   -0.01116457	,   -0.893536360 ,
										   -0.00936052	,   -0.892791007 ,
										   -0.00887078	,   -0.892588293 ,
										   -0.00806020	,   -0.892252427 ,
										   -0.00728260	,   -0.891929813 ,
										   -0.00685426	,   -0.891751931 ,
										   -0.00606659	,   -0.891424503 ,
										   -0.00593068	,   -0.891367968 ,
										   -0.00506905	,   -0.891009245 ,
										   -0.00457149	,   -0.890801872 ,
										   -0.00326904	 ,-0.890258262 ,
										   -0.00249866	 , -0.889936195 ,
										   -0.00119656	 ,  -0.889390941 
										     };
												   for (size_t i=0; i<90; i++)
												   {
													   if (ratio > Ftable[i][0])
													   {
														   Fratio=Ftable[i][1]+(Ftable[i+1][1]-Ftable[i][1])/(Ftable[i+1][0]-Ftable[i][0])*(ratio-Ftable[i][0]);
														   //std::cout << " i=	" <<i<< " Ftable=  " << Ftable[i][1]<< " ratio=    " << ratio<< std::endl;
														   continue;
													   }
											   }
							  Vec3_t Fne = Fratio*Fc*n; //Modification 1   attention: the force direction
							   Vec3_t Fnd = Gn*dot(n,vrel)*n;
												  
					
       
		     
							  
							   Vec3_t F = Fne + Fnd ; 
							  F1  = -F;
							  F2  =  F;

					 
							
							  
							 
									
									
							 
		 Vec3_t T, Tt;
		 
        Tt = cross (x1,F);  //
        Quaternion_t q;
        Conjugate (P1->Q,q);
        Rotation  (Tt,q,T);
        T1 = -T;
		
        Tt = cross (x2,F) ;  // 
        Conjugate (P2->Q,q);
        Rotation  (Tt,q,T);
        T2 = T;
				}
				}
		if(delta>0.0){

		    First=false;
			//double a0	= pow((9*M_PI*Gama*R*R/E),(1.0/3.0)); //equilibrium contact radius
			//double delta_c= a0*a0/(2.0*R*pow(6.0,1.0/3.0)); // critical overlap
			double a	= pow(delta*R,0.5); //contact radius	

	
			
			
			double e_alpha =0.8;
			double alpha=1.2728-4.2783*e_alpha+11.087*pow(e_alpha,2)-22.348*pow(e_alpha,3)+27.467*pow(e_alpha,4)-18.022*pow(e_alpha,5)+4.8218*pow(e_alpha,6);//#########################################################
			
			 Kt 	  = 8.0*G*a; //tangential stiffness coefficient
			double kn	  = 4.0*E*a/3.0;
			Gn	   = alpha*pow(P1->Props.m*kn,0.5); //normal friction/damping coefficient
			Gt	   = Gn;  //tangential friction/damping coefficient

			 
			//Vec3_t Fc	= 3*M_PI*Gama*R*n; // critical force
			double Fc	= 3*M_PI*Gama*R; // critical force
			 double ratio = delta/delta_c;
		     double Fratio=0.0;
							   
							   //std::cout << "delta = " << delta<< "  delta_c =" << delta_c<< std::endl;
						
							   //std::cout << "R = " << R<< "   ratio =" << ratio<< std::endl;
							   
						   //  std::cout << "a = " << a<< "Gama= " << Gama<< "Kt= " << Kt<< "	kn=" << kn<< "		   Gn= " << kn<< "	  E=	 " <<E<< "	  G=	 " <<G<< std::endl;
						
						
							    double Ftable[461][2]={ 
							   -1.0,      -0.95 ,
							   -0.7233934943 , -0.956312873 ,
							   -0.7120345372 , -0.960773949,
							   -0.6111813145 ,	-0.988620724,
							   -0.5622965558 ,	-0.995718097,
							   -0.5124456352 ,	-0.999375691,
							   -0.4785866453 ,-0.999997157,
							   -0.4208657037 ,	-0.997897821,
							   -0.3718303641 ,	-0.993239474,
							   -0.3209987346 ,-0.985851105,
							   -0.2862049717 ,-0.979388125,
							   -0.2362498370 ,-0.968229621,
							   -0.1863787347 ,-0.95499984,
							   -0.1244077761 ,-0.935824917,
							   -0.0854958065 ,-0.922316149,
							   -0.0440067886 ,-0.906721951,
							   -0.0410077963 ,-0.905548342,
							   -0.0351212136 ,-0.903226807,
							   -0.0311653133 ,-0.901653405,
							   -0.0255545350 ,-0.899403571,
							   -0.0222252979 ,-0.898058529,
							   -0.0194314907 ,-0.896924039,
							   -0.0168916310 ,-0.895888113,
							   -0.0125633717 ,-0.894112789,
							   -0.0094744550 ,-0.892838144,
							   -0.0064534441 ,-0.891585367,
							   -0.0033922364 ,-0.89030973,
							   -0.0030641204 ,-0.890172631,
							   -0.0026476899 ,-0.889998529,
							   0.000653961 ,   -0.888614101    ,
							   0.001070147 ,   -0.888439076    ,
							   0.001486415 ,   -0.888263903    ,
							   0.001902765 ,   -0.888088581    ,
							   0.002319197 ,   -0.88791311 ,
							   0.00273571  ,   -0.88773749 ,
							   0.003152305 ,   -0.887561721    ,
							   0.003568982 ,   -0.887385802    ,
							   0.00398574  ,   -0.887209735    ,
							   0.00440258  ,   -0.887033519    ,
							   0.004819502 ,   -0.886857154    ,
							   0.005236506 ,   -0.886680639    ,
							   0.005653591 ,   -0.886503976    ,
							   0.006070759 ,   -0.886327163    ,
							   0.006488008 ,   -0.886150201    ,
							   0.006905338 ,   -0.88597309 ,
							   0.007322751 ,   -0.885795829    ,
							   0.007740245 ,   -0.88561842 ,
							   0.008157821 ,   -0.885440861    ,
							   0.008575479 ,   -0.885263152    ,
							   0.008993218 ,   -0.885085295    ,
							   0.009411039 ,   -0.884907288    ,
							   0.009828942 ,   -0.884729131    ,
							   0.010246927 ,   -0.884550825    ,
							   0.010664993 ,   -0.88437237 ,
							   0.011083141 ,   -0.884193765    ,
							   0.011501371 ,   -0.884015011    ,
							   0.011919682 ,   -0.883836107    ,
							   0.012338076 ,   -0.883657054    ,
							   0.012756551 ,   -0.883477851    ,
							   0.013175107 ,   -0.883298498    ,
							   0.013593746 ,   -0.883118996    ,
							   0.014012466 ,   -0.882939344    ,
							   0.014431267 ,   -0.882759542    ,
							   0.014850151 ,   -0.882579591    ,
							   0.015269116 ,   -0.88239949 ,
							   0.015688163 ,   -0.882219239    ,
							   0.016107292 ,   -0.882038838    ,
							   0.016526502 ,   -0.881858288    ,
							   0.016945794 ,   -0.881677587    ,
							   0.017365168 ,   -0.881496737    ,
							   0.017784623 ,   -0.881315737    ,
							   0.01820416  ,   -0.881134587    ,
							   0.018623779 ,   -0.880953287    ,
							   0.01904348  ,   -0.880771837    ,
							   0.019463262 ,   -0.880590237    ,
							   0.019883126 ,   -0.880408487    ,
							   0.020303071 ,   -0.880226587    ,
							   0.020723098 ,   -0.880044536    ,
							   0.021143207 ,   -0.879862336    ,
							   0.021563398 ,   -0.879679986    ,
							   0.02198367  ,   -0.879497485    ,
							   0.022404024 ,   -0.879314834    ,
							   0.02282446  ,   -0.879132033    ,
							   0.023244977 ,   -0.878949082    ,
							   0.023665576 ,   -0.87876598 ,
							   0.024086257 ,   -0.878582729    ,
							   0.024507019 ,   -0.878399326    ,
							   0.024927863 ,   -0.878215774    ,
							   0.025348789 ,   -0.878032071    ,
							   0.025769796 ,   -0.877848218    ,
							   0.026190885 ,   -0.877664214    ,
							   0.026612056 ,   -0.87748006 ,
							   0.027033308 ,   -0.877295755    ,
							   0.027454642 ,   -0.8771113  ,
							   0.027876058 ,   -0.876926694    ,
							   0.028297555 ,   -0.876741938    ,
							   0.028719134 ,   -0.876557031    ,
							   0.029140794 ,   -0.876371974    ,
							   0.029562536 ,   -0.876186766    ,
							   0.02998436  ,   -0.876001407    ,
							   0.030828253 ,   -0.875630237    ,
							   0.031672472 ,   -0.875258465    ,
							   0.032517018 ,   -0.874886089    ,
							   0.03336189  ,   -0.87451311 ,
							   0.034207089 ,   -0.874139527    ,
							   0.035052614 ,   -0.873765341    ,
							   0.035898466 ,   -0.87339055 ,
							   0.036744644 ,   -0.873015155    ,
							   0.037591149 ,   -0.872639156    ,
							   0.03843798  ,   -0.872262552    ,
							   0.039285138 ,   -0.871885343    ,
							   0.040132622 ,   -0.87150753 ,
							   0.040980432 ,   -0.871129111    ,
							   0.041828569 ,   -0.870750086    ,
							   0.042677033 ,   -0.870370456    ,
							   0.043525822 ,   -0.869990219    ,
							   0.044374938 ,   -0.869609377    ,
							   0.04522438  ,   -0.869227929    ,
							   0.046074149 ,   -0.868845873    ,
							   0.046924244 ,   -0.868463211    ,
							   0.047774665 ,   -0.868079943    ,
							   0.048625413 ,   -0.867696067    ,
							   0.049476487 ,   -0.867311583    ,
							   0.050327887 ,   -0.866926492    ,
							   0.051179613 ,   -0.866540793    ,
							   0.052031666 ,   -0.866154487    ,
							   0.052884044 ,   -0.865767572    ,
							   0.053736749 ,   -0.865380048    ,
							   0.05458978  ,   -0.864991916    ,
							   0.055443138 ,   -0.864603175    ,
							   0.056296821 ,   -0.864213825    ,
							   0.057150831 ,   -0.863823866    ,
							   0.058005167 ,   -0.863433297    ,
							   0.058859828 ,   -0.863042119    ,
							   0.059714816 ,   -0.862650331    ,
							   0.061425771 ,   -0.861864924    ,
							   0.063138029 ,   -0.861077075    ,
							   0.064851592 ,   -0.860286782    ,
							   0.066566459 ,   -0.859494043    ,
							   0.068282629 ,   -0.858698858    ,
							   0.070000104 ,   -0.857901224    ,
							   0.071718882 ,   -0.857101139    ,
							   0.073438964 ,   -0.856298603    ,
							   0.075160349 ,   -0.855493612    ,
							   0.076883038 ,   -0.854686166    ,
							   0.07860703  ,   -0.853876263    ,
							   0.080332325 ,   -0.853063901    ,
							   0.082058923 ,   -0.852249079    ,
							   0.083786825 ,   -0.851431794    ,
							   0.085516029 ,   -0.850612045    ,
							   0.087246535 ,   -0.849789831    ,
							   0.088978345 ,   -0.84896515	   ,
							   0.090711457 ,   -0.848138	   ,
							   0.092445871 ,   -0.847308379    ,
							   0.094181588 ,   -0.846476287    ,
							   0.095918607 ,   -0.84564172	   ,
							   0.097656928 ,   -0.844804677    ,
							   0.099396551 ,   -0.843965158    ,
							   0.103751305 ,   -0.84185551	   ,
							   0.108114193 ,   -0.839730343    ,
							   0.112485215 ,   -0.83758963	   ,
							   0.11686437  ,   -0.835433346    ,
							   0.121251654 ,   -0.833261463    ,
							   0.125647067 ,   -0.831073957    ,
							   0.130050607 ,   -0.828870801    ,
							   0.134462272 ,   -0.826651968    ,
							   0.138882062 ,   -0.824417434    ,
							   0.143309973 ,   -0.822167171    ,
							   0.147746005 ,   -0.819901153    ,
							   0.152190156 ,   -0.817619355    ,
							   0.156642424 ,   -0.815321749    ,
							   0.179025471 ,   -0.803595709    ,
							   0.201611215 ,   -0.791470587    ,
							   0.224399463 ,   -0.778943124    ,
							   0.247390023 ,   -0.766010064    ,
							   0.270582706 ,   -0.752668154    ,
							   0.293977329 ,   -0.73891414	   ,
							   0.317573707 ,   -0.724744775    ,
							   0.341371662 ,   -0.710156811    ,
							   0.365371015 ,   -0.695147004    ,
							   0.389571593 ,   -0.679712109    ,
							   0.413973223 ,   -0.663848888    ,
							   0.438575736 ,   -0.647554102    ,
							   0.463378964 ,   -0.630824514    ,
							   0.488382744 ,   -0.613656889    ,
							   0.513586913 ,   -0.596047997    ,
							   0.53899131  ,   -0.577994606    ,
							   0.564595779 ,   -0.559493488    ,
							   0.590400164 ,   -0.540541417    ,
							   0.616404311 ,   -0.521135168    ,
							   0.642608069 ,   -0.501271519    ,
							   0.66901129  ,   -0.480947249    ,
							   0.695613826 ,   -0.460159139    ,
							   0.722415532 ,   -0.438903973    ,
							   0.749416266 ,   -0.417178534    ,
							   0.776615886 ,   -0.394979609    ,
							   0.804014254 ,   -0.372303987    ,
							   0.831611231 ,   -0.349148458    ,
							   0.859406683 ,   -0.325509812    ,
							   0.887400476 ,   -0.301384844    ,
							   0.915592478 ,   -0.276770348    ,
							   0.943982558 ,   -0.251663121    ,
							   0.972570589 ,   -0.22605996	   ,
							   1.001356444 ,   -0.199957667    ,
							   1.030339997 ,   -0.173353041    ,
							   1.059521126 ,   -0.146242886    ,
							   1.088899707 ,   -0.118624006    ,
							   1.118475621 ,   -0.090493207    ,
							   1.148248749 ,   -0.061847296    ,
							   1.178218973 ,   -0.032683082    ,
							   1.208386179 ,   -0.002997376    ,
							   1.211413729 ,   0			   ,
							   1.24179748  ,   0.030263031 ,
							   1.272377974 ,   0.061054249 ,
							   1.303155098 ,   0.092376839 ,
							   1.334128743 ,   0.124233985 ,
							   1.365298799 ,   0.15662887  ,
							   1.39666516  ,   0.189564675 ,
							   1.428227718 ,   0.22304458  ,
							   1.45998637  ,   0.257071765 ,
							   1.491941012 ,   0.291649405 ,
							   1.52409154  ,   0.326780678 ,
							   1.556437856 ,   0.362468759 ,
							   1.588979857 ,   0.39871682  ,
							   1.621717446 ,   0.435528035 ,
							   1.654650524 ,   0.472905575 ,
							   1.687778997 ,   0.510852609 ,
							   1.721102767 ,   0.549372307 ,
							   1.754621741 ,   0.588467836 ,
							   1.788335826 ,   0.628142362 ,
							   1.822244929 ,   0.668399052 ,
							   1.85634896  ,   0.709241068 ,
							   1.890647829 ,   0.750671575 ,
							   1.925141446 ,   0.792693734 ,
							   1.959829724 ,   0.835310706 ,
							   1.994712575 ,   0.878525651 ,
							   2.029789914 ,   0.922341727 ,
							   2.065061655 ,   0.966762093 ,
							   2.100527715 ,   1.011789904 ,
							   2.13618801  ,   1.057428317 ,
							   2.172042458 ,   1.103680486 ,
							   2.208090977 ,   1.150549564 ,
							   2.244333487 ,   1.198038705 ,
							   2.280769908 ,   1.246151059 ,
							   2.317400162 ,   1.294889777 ,
							   2.354224171 ,   1.344258009 ,
							   2.391241856 ,   1.394258904 ,
							   2.428453143 ,   1.444895608 ,
							   2.465857956 ,   1.496171269 ,
							   2.50345622  ,   1.548089033 ,
							   2.54124786  ,   1.600652045 ,
							   2.579232805 ,   1.653863448 ,
							   2.617410981 ,   1.707726386 ,
							   2.655782317 ,   1.762244    ,
							   2.694346742 ,   1.817419433 ,
							   2.733104186 ,   1.873255824 ,
							   2.77205458  ,   1.929756313 ,
							   2.811197854 ,   1.986924039 ,
							   2.85053394  ,   2.04476214  ,
							   2.890062772 ,   2.103273752 ,
							   2.929784283 ,   2.162462012 ,
							   2.969698407 ,   2.222330056 ,
							   3.009805077 ,   2.282881018 ,
							   3.050104231 ,   2.344118031 ,
							   3.090595803 ,   2.406044229 ,
							   3.13127973  ,   2.468662744 ,
							   3.17215595  ,   2.531976707 ,
							   3.213224399 ,   2.595989249 ,
							   3.254485018 ,   2.660703499 ,
							   3.295937744 ,   2.726122587 ,
							   3.337582517 ,   2.792249641 ,
							   3.379419278 ,   2.859087789 ,
							   3.421447967 ,   2.926640158 ,
							   3.463668526 ,   2.994909873 ,
							   3.506080896 ,   3.063900061 ,
							   3.54868502  ,   3.133613845 ,
							   3.591480841 ,   3.20405435  ,
							   3.634468303 ,   3.2752247   ,
							   3.677647349 ,   3.347128016 ,
							   3.721017924 ,   3.41976742  ,
							   3.764579974 ,   3.493146034 ,
							   3.808333443 ,   3.567266979 ,
							   3.852278279 ,   3.642133374 ,
							   3.896414427 ,   3.717748338 ,
							   3.940741834 ,   3.794114991 ,
							   3.985260449 ,   3.871236449 ,
							   4.02997022  ,   3.94911583  ,
							   4.074871094 ,   4.027756251 ,
							   4.119963021 ,   4.107160828 ,
							   4.16524595  ,   4.187332676 ,
							   4.210719832 ,   4.268274909 ,
							   4.256384617 ,   4.349990643 ,
							   4.302240254 ,   4.43248299  ,
							   4.348286697 ,   4.515755063 ,
							   4.394523896 ,   4.599809975 ,
							   4.440951803 ,   4.684650837 ,
							   4.487570371 ,   4.77028076  ,
							   4.534379554 ,   4.856702855 ,
							   4.581379303 ,   4.943920232 ,
							   4.628569574 ,   5.031936    ,
							   4.67595032  ,   5.120753268 ,
							   4.723521496 ,   5.210375144 ,
							   4.771283056 ,   5.300804736 ,
							   4.819234956 ,   5.392045151 ,
							   4.867377152 ,   5.484099496 ,
							   4.915709601 ,   5.576970876 ,
							   4.964232257 ,   5.670662397 ,
							   5.012945078 ,   5.765177164 ,
							   5.061848022 ,   5.860518281 ,
							   5.110941046 ,   5.956688853 ,
							   5.160224108 ,   6.053691982 ,
							   5.209697165 ,   6.151530772 ,
							   5.259360178 ,   6.250208324 ,
							   5.309213104 ,   6.349727741 ,
							   5.359255903 ,   6.450092123 ,
							   5.409488535 ,   6.551304572 ,
							   5.459910959 ,   6.653368189 ,
							   5.510523136 ,   6.756286071 ,
							   5.561325027 ,   6.86006132  ,
							   5.612316592 ,   6.964697034 ,
							   5.663497792 ,   7.070196311 ,
							   5.71486859  ,   7.176562249 ,
							   5.766428947 ,   7.283797946 ,
							   5.818178825 ,   7.391906498 ,
							   5.870118187 ,   7.500891002 ,
							   5.922246996 ,   7.610754554 ,
							   5.974565214 ,   7.721500249 ,
							   6.027072806 ,   7.833131183 ,
							   6.079769734 ,   7.94565045  ,
							   6.132655962 ,   8.059061145 ,
							   6.185731456 ,   8.17336636  ,
							   6.238996178 ,   8.28856919  ,
							   6.292450095 ,   8.404672727 ,
							   6.346093172 ,   8.521680064 ,
							   6.399925372 ,   8.639594292 ,
							   6.453946663 ,   8.758418504 ,
							   6.50815701  ,   8.87815579  ,
							   6.56255638  ,   8.998809242 ,
							   6.617144737 ,   9.120381948 ,
							   6.67192205  ,   9.242877    ,
							   6.726888285 ,   9.366297487 ,
							   6.78204341  ,   9.490646498 ,
							   6.83738739  ,   9.615927121 ,
							   6.892920196 ,   9.742142445 ,
							   6.948641793 ,   9.869295558 ,
							   7.004552151 ,   9.997389547 ,
							   7.060651237 ,   10.1264275  ,
							   7.116939021 ,   10.2564125  ,
							   7.173415471 ,   10.38734764 ,
							   7.230080556 ,   10.519236   ,
							   7.286934246 ,   10.65208067 ,
							   7.343976509 ,   10.78588473 ,
							   7.401207317 ,   10.92065127 ,
							   7.458626638 ,   11.05638337 ,
							   7.516234443 ,   11.19308411 ,
							   7.574030703 ,   11.33075659 ,
							   7.632015387 ,   11.46940388 ,
							   7.690188467 ,   11.60902906 ,
							   7.748549915 ,   11.74963523 ,
							   7.8070997   ,   11.89122545 ,
							   7.865837795 ,   12.03380282 ,
							   7.924764171 ,   12.17737041 ,
							   7.9838788   ,   12.32193131 ,
							   8.043181654 ,   12.46748859 ,
							   8.102672706 ,   12.61404534 ,
							   8.162351927 ,   12.76160464 ,
							   8.222219291 ,   12.91016956 ,
							   8.282274771 ,   13.05974319 ,
							   8.342518339 ,   13.2103286  ,
							   8.402949968 ,   13.36192888 ,
							   8.463569633 ,   13.5145471  ,
							   8.524377306 ,   13.66818634 ,
							   8.585372962 ,   13.82284968 ,
							   8.646556574 ,   13.97854019 ,
							   8.707928117 ,   14.13526096 ,
							   8.769487565 ,   14.29301505 ,
							   8.831234892 ,   14.45180556 ,
							   8.893170073 ,   14.61163554 ,
							   8.955293083 ,   14.77250809 ,
							   9.017603897 ,   14.93442626 ,
							   9.08010249  ,   15.09739315 ,
							   9.142788837 ,   15.26141182 ,
							   9.205662915 ,   15.42648535 ,
							   9.268724698 ,   15.59261682 ,
							   9.331974162 ,   15.75980929 ,
							   9.395411285 ,   15.92806584 ,
							   9.45903604  ,   16.09738955 ,
							   9.522848406 ,   16.26778348 ,
							   9.586848358 ,   16.43925072 ,
							   9.651035874 ,   16.61179432 ,
							   9.715410929 ,   16.78541738 ,
							   9.779973502 ,   16.96012295 ,
							   9.844723568 ,   17.13591411 ,
							   9.909661105 ,   17.31279393 ,
							   9.974786092 ,   17.49076548 ,
							   10.0400985  ,   17.66983184 ,
							   10.70352623 ,   19.52138495 ,
							   11.38567292 ,   21.48579431 ,
							   12.08651837 ,   23.5661265  ,
							   12.80604356 ,   25.76544563 ,
							   13.54423061 ,   28.08681356 ,
							   14.30106265 ,   30.53328999 ,
							   15.07652377 ,   33.10793259 ,
							   15.87059892 ,   35.81379712 ,
							   16.68327386 ,   38.65393751 ,
							   17.51453509 ,   41.63140599 ,
							   18.36436981 ,   44.74925316 ,
							   19.23276587 ,   48.01052807 ,
							   20.11971172 ,   51.41827831 ,
							   21.02519635 ,   54.97555007 ,
							   21.9492093  ,   58.68538821 ,
							   22.89174058 ,   62.55083632 ,
							   23.85278066 ,   66.57493679 ,
							   24.83232044 ,   70.76073086 ,
							   25.83035121 ,   75.11125865 ,
							   26.84686467 ,   79.62955925 ,
							   27.88185284 ,   84.31867073 ,
							   28.93530809 ,   89.18163019 ,
							   30.0072231  ,   94.22147381 ,
							   31.09759085 ,   99.44123688 ,
							   32.2064046  ,   104.8439538 ,
							   33.33365788 ,   110.4326583 ,
							   34.47934444 ,   116.2103832 ,
							   35.64345832 ,   122.1801604 ,
							   36.82599374 ,   128.3450215 ,
							   38.02694515 ,   134.7079969 ,
							   39.24630721 ,   141.2721168 ,
							   40.48407475 ,   148.0404104 ,
							   41.74024281 ,   155.0159065 ,
							   43.01480657 ,   162.2016331 ,
							   44.30776141 ,   169.6006177 ,
							   45.61910285 ,   177.2158873 ,
							   46.94882656 ,   185.0504683 ,
							   48.29692835 ,   193.1073866 ,
							   49.66340418 ,   201.3896675 ,
							   52.45146239 ,   218.6424162 ,
							   55.31297131 ,   236.8329076 ,
							   58.24790288 ,   255.9853281 ,
							   61.25623071 ,   276.123857  ,
							   64.33792995 ,   297.2726674 ,
							   67.49297716 ,   319.4559264 ,
							   70.72135021 ,   342.6977954 ,
							   74.02302817 ,   367.0224306 ,
							   77.39799123 ,   392.453983  ,
							   80.84622061 ,   419.016599  ,
							   84.36769848 ,   446.7344203 ,
							   87.96240792 ,   475.6315845 ,
							   91.63033283 ,   505.7322249 ,
							   95.37145788 ,   537.060471  ,
							   99.18576848 ,   569.6404485 ,
							   103.0732507 ,   603.4962799 , 

							   125.0000639 ,   806.2004    ,							   
							   150.0000044 ,   1059.989    ,
							   174.9998873 ,   1335.913    ,
						       199.9998    ,   1632.321    ,
						       229.9999886 ,   2013.199    ,
						       255.9999928 ,   2364.157    ,
						       279.9999972 ,   2704.387    ,
						       300.9999992 ,   3014.343    ,
						       399.9999606 ,   4618.133    ,
						       499.9991627 ,   6454.288    ,


							   };
									   for (size_t i=0; i<461; i++)
									   {
										   if (ratio > Ftable[i][0])
										   {
											   Fratio=Ftable[i][1]+(Ftable[i+1][1]-Ftable[i][1])/(Ftable[i+1][0]-Ftable[i][0])*(ratio-Ftable[i][0]);
											   //std::cout << " i=	" <<i<< " Ftable=  " << Ftable[i][1]<< " ratio=    " << ratio<< std::endl;
											   continue;
										   }
								   }
				  Vec3_t Fne = Fratio*Fc*n; //Modification 1	  attention: the force direction
				  Vec3_t Fnd = Gn*dot(n,vrel)*n;
				  Vec3_t Fn = Fne+Fnd;


		        Vec3_t ts = vt/norm(vt);
	        x_t     +=(vt*dt);
			x_t     -= dot(x_t,n)*n;
			Vec3_t Ft =Kt* x_t+Gt*vt;
	        if (norm(Ft)>Mu*(norm(Fn)+2*Fc))  //Mu=0.3
	        {
	            
	            Ft = Mu*(norm(Fn)+2*Fc)*ts;
	        }       

		         Vec3_t F = Fn + Ft;
				 F1	= -F;
				 F2	=  F;
			
/*
		   SFr 	 += dt*vt;
		   SFr		-= dot(SFr,n)*n;
		   Vec3_t tan= SFr;
		   if(norm(tan)>0.0) tan/=norm(tan);
		   if(norm(SFr)>Mu*(norm(Fn)+2*Fc)/Kt)
		   {
			   SFr = Mu*(norm(Fn)+2*Fc)/Kt*tan;
		   }

         Vec3_t F = Fne + Fnd + Kt*SFr + Gt*vt;
        

 */

		
		 
		 //Rolling resistance
		 double miu_R=2/3*Mu;
		 double yita_R=miu_R*norm(Fn);
		 double KR = 4*Fc*pow((a/a0),1.5);
		 
		 Vec3_t Vr = P1->Props.R*P2->Props.R*cross(Vec3_t(t1 - t2),n)/(P1->Props.R+P2->Props.R);
		 Fdr += Vr*dt;
		 Fdr -= dot(Fdr,n)*n;
		 Vec3_t tanr = Fdr;
		 if (norm(tanr)>0.0) tanr/=norm(tanr);
		 Vec3_t Mr = KR*Fdr+yita_R*Vr;
	
		  if (norm(Mr)>KR*0.03*2*M_PI*R)  
		 {
			 Mr = KR*0.03*2*M_PI*R*tanr;
			
		 }
   

		 

		 /*******************************      twisting torque           *******************************/
			//Vec3_t omiga_n = dot(n,(P1->w-P2->w))*n;
			
        Vec3_t omiga_n = dot(n,(t1-t2))*n;
		 thei_tw+=(omiga_n*dt);
        Vec3_t Mt = ((Kt*a*a)/2)*(thei_tw)  +  ((Gn*a*a)/2)*omiga_n;
		 if (norm(Mt)>2/3*a*Mu*(norm(Fne)+2*Fc))	//Mu=0.3
				{
					Mt = 2/3*a*Mu*(norm(Fne)+2*Fc)*omiga_n/norm(omiga_n);
				}
				/*****************************************************/

		                               


if(Util::IsNan(norm(F))) 
								   {
							
									std::cout << "Fnesphere = " << Fne		<< std::endl;

									std::cout << "P1->v = " << P1->v<< "P2->v = " << P2->v	<< std::endl;

									std::cout << "t1 = " << t1<< "t2 = " << t2		<< std::endl;

									std::cout << "vrel = " << vrel		<< std::endl;

									std::cout << "n = " << n		<< std::endl;
									
									std::cout << "Fnd = " << Fnd		<< std::endl;
									
									std::cout << "Ft = " << Ft		<< std::endl;
													  
								   }



		  /*******************************      rolling torque           *******************************/
/*
		 double miu_R=2/3*Mu;
		 double yita_R=miu_R*norm(Fne);
		 
         double KR = 4*Fc*pow((a/a0),1.5);

		Vec3_t omiga = t1-t2;
		Vec3_t Vl = -R*cross(omiga,n);	//??????????//negative ???-0.5*((P2->Props.R-P1->Props.R)/(P2->Props.R+P1->Props.R))*vt
		Vec3_t tR =  Vl/norm(Vl);
		 
		  thei_t+=(Vl*dt);
		 Vec3_t Mr = KR*thei_t;       //+  yita_R*Vl
		 
		 //19–64 milliradians=(19–64)*0.0573
		 //double angle = 60*0.0573*M_PI/180.0;
        if (norm(thei_t)>0.01*2*M_PI*P1->Props.R)	//Mu=0.3      20*10-3rad
				{
		    	 Mr = (KR*0.01*2*M_PI*P1->Props.R)*tR;		
				}
				
				/*****************************************************/


        Vec3_t T, Tt;
		 
        Tt = cross (x1,F) + Mt + P1->Props.R*cross(n,Mr); //+ cross(Mr,n)
        Quaternion_t q;
        Conjugate (P1->Q,q);
        Rotation  (Tt,q,T);
        T1 = -T;
		
        Tt = cross (x2,F) + Mt+ P2->Props.R*cross(n,Mr); // +cross(Mr,n)
        Conjugate (P2->Q,q);
        Rotation  (Tt,q,T);
        T2 = T;
			}//##########################
        }

   
    return false;
}

inline bool CInteractonSphere::UpdateContacts (double alpha)
{
	//std::cout << "**************************9*******************************" << std::endl;


		
		


	if(P1->Props.fixedparticle==true)
		if(P2->Props.fixedparticle==true)  return false;

	
    if (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*alpha) 
    	{
	        P1->Props.interact=true;//2022.03.26  by qian
			P2->Props.interact=true;//2022.03.26  by qian

			
			if(P1->Props.fixedparticle==true||P2->Props.fixedparticle==true)//2022.04.01
						{
						P1->Props.ready_fixedparticle=true;
							P2->Props.ready_fixedparticle=true;
						}

			
	        return true;
    	}
    else 
		{
	P1->Props.interact=false;//2022.03.26  by qian
    P2->Props.interact=false;//2022.03.26  by qian

	P1->Props.ready_fixedparticle=false;//2022.04.01
	P2->Props.ready_fixedparticle=false;//2022.04.01
				
	return false;
	}
}

inline void CInteractonSphere::UpdateParameters ()
{
	//std::cout << "**************************10*******************************" << std::endl;

    Kn   = 2*ReducedValue(P1->Props.Kn,P2->Props.Kn);
    Kt   = 2*ReducedValue(P1->Props.Kt,P2->Props.Kt);
    double me       = ReducedValue(P1->Props.m ,P2->Props.m );
    Gn              = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn);
    Gt              = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt);
    if (Gn < -0.001)
    {
        if (fabs(Gn)>1.0) throw new Fatal("CInteracton the restitution coefficient is greater than 1");
        Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
        Gt = 0.0;
    }
    Gn *= me;
    Gt *= me;
    //Mu              = 2*ReducedValue(P1->Props.Mu,P2->Props.Mu);
    if (P1->Props.Mu>1.0e-12&&P2->Props.Mu>1.0e-12)
    {
        if (!P1->IsFree())      Mu = P1->Props.Mu;
        else if (!P2->IsFree()) Mu = P2->Props.Mu;
        else                    Mu = std::max(P1->Props.Mu,P2->Props.Mu);
    }
    else 
    {
        Mu          = 0.0;
    }
    beta = 2*ReducedValue(P1->Props.Beta,P2->Props.Beta);
    eta  = 2*ReducedValue(P1->Props.Eta,P2->Props.Eta);
    Bn   = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn);
    Bt   = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt);
    eps  = 2*ReducedValue(P1->Props.eps,P2->Props.eps);
}

//Cohesion interacton

inline BInteracton::BInteracton (Particle * Pt1, Particle * Pt2, size_t Fi1, size_t Fi2)
{
	//std::cout << "**************************11*******************************" << std::endl;

    P1              = Pt1;
    P2              = Pt2;
    IF1             = Fi1;
    if (Fi2>=P2->Faces.Size())
    {
        IF2 = P2->Faces.Size()-1;
        Area            = P1->Faces[IF1]->Area();
    }
    else
    {
        IF2 = Fi2;
        Area            = 0.5*(P1->Faces[IF1]->Area()+P2->Faces[IF2]->Area());
    }
    I1              = P1->Index;
    I2              = P2->Index;
    Bn              = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn)*Area;
    Bt              = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt)*Area;
    Bm              = 2*ReducedValue(P1->Props.Bm,P2->Props.Bm)*Area;
    double me       = ReducedValue(Pt1->Props.m ,Pt2->Props.m );
    Gn              = 2*ReducedValue(Pt1->Props.Gn,Pt2->Props.Gn)*me;
    Gt              = 2*ReducedValue(Pt1->Props.Gt,Pt2->Props.Gt)*me;
    eps             = 2*ReducedValue(P1->Props.eps,P2->Props.eps);
    eta             = 2*ReducedValue(P1->Props.Eta,P2->Props.Eta);

    Vec3_t n1,n2;
    P1->Faces[IF1]->Normal(n1);
    P2->Faces[IF2]->Normal(n2);
    Vec3_t c1,c2;
    An         = 0.0;
    valid      = true;
    P1->Faces[IF1]->Centroid(c1);
    P2->Faces[IF2]->Centroid(c2);
    L0         = dot(n1,c2-c1);
    if (Fi2>=P2->Faces.Size()) c2 = c1;
    Vec3_t V   = 0.5*(c1+c2);

    Face * F   = P1->Faces[IF1];
    double a   = dot(*F->Edges[0]->X0-V , F->Edges[0]->dL);
    double b   = dot(F->Edges[0]->dL    , F->Edges[0]->dL);
    double c   = dot(F->Edges[0]->dL    , F->Edges[1]->dL);
    double d   = dot(*F->Edges[0]->X0-V , F->Edges[1]->dL);
    double f   = dot(F->Edges[1]->dL    , F->Edges[1]->dL);
    s1         = (c*d-a*f)/(b*f-c*c);
    t1         = (a*c-b*d)/(b*f-c*c);
    
    Vec3_t p1  = -s1*F->Edges[0]->dL - t1*F->Edges[1]->dL;

    F          = P2->Faces[IF2];
    a          = dot(*F->Edges[0]->X0-V , F->Edges[0]->dL);
    b          = dot(F->Edges[0]->dL    , F->Edges[0]->dL);
    c          = dot(F->Edges[0]->dL    , F->Edges[1]->dL);
    d          = dot(*F->Edges[0]->X0-V , F->Edges[1]->dL);
    f          = dot(F->Edges[1]->dL    , F->Edges[1]->dL);
    s2         = (c*d-a*f)/(b*f-c*c);
    t2         = (a*c-b*d)/(b*f-c*c);

    Vec3_t p2  = -s2*F->Edges[0]->dL - t2*F->Edges[1]->dL;
    
    double arg = dot(p1,p2)/(norm(p1)*norm(p2));
    if (arg> 1.0) arg= 1.0;
    if (arg<-1.0) arg=-1.0;
    An         = acos(arg);
    if (dot(cross(p1,p2),n1)<0) An = -An;
#ifdef USE_THREAD
    pthread_mutex_init(&lck,NULL);
#endif
}

inline bool BInteracton::UpdateContacts (double alpha)
{
    return valid;
}

inline bool BInteracton::CalcForce(double dt)
{

	//std::cout << "**************************12*******************************" << std::endl;

    F1 = OrthoSys::O;
    F2 = OrthoSys::O;
    T1 = OrthoSys::O;
    T2 = OrthoSys::O;
    if (valid)
    {
        // Calculate the normal vector and centroid of the contact face
        Vec3_t n1,n2,n;
        P1->Faces[IF1]->Normal(n1);
        P2->Faces[IF2]->Normal(n2);
        n = 0.5*(n1-n2);
        n/=norm(n);

        Face * F    = P1->Faces[IF1];
        Vec3_t pro1 = *F->Edges[0]->X0 + s1*F->Edges[0]->dL + t1*F->Edges[1]->dL;
        Vec3_t p1   = -s1*F->Edges[0]->dL - t1*F->Edges[1]->dL;
        F           = P2->Faces[IF2];
        Vec3_t pro2 = *F->Edges[0]->X0 + s2*F->Edges[0]->dL + t2*F->Edges[1]->dL;
        Vec3_t p2   = -s2*F->Edges[0]->dL - t2*F->Edges[1]->dL;

        // force
        double delta = (dot(pro2-pro1,n)-L0)/L0;
        if (delta<0.0) delta = 0.0;

        xnet        = 0.5*(pro1+pro2);
        Vec3_t t1,t2,x1,x2;
        Rotation(P1->w,P1->Q,t1);
        Rotation(P2->w,P2->Q,t2);
        x1 = xnet - P1->x;
        x2 = xnet - P2->x;
        Vec3_t td   = pro2-pro1-dot(pro2-pro1,n)*n;
        //Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
        //Vec3_t vt   = vrel - dot(n,vrel)*n;
        //Fnet        = -Bn*(delta)*n - Gn*dot(n,vrel)*n;
        //Ftnet       = -Bt*td/L0 - Gt*vt;
        Fnet        = -Bn*(delta)*n;
        Ftnet       = -Bt*td/L0;

#ifdef USE_THREAD
        //pthread_mutex_lock(&lck);
        //pthread_mutex_lock(&P1->lck);
        //pthread_mutex_lock(&P2->lck);
        //std::lock_guard<std::mutex> lk1(P1->mtex);
        //std::lock_guard<std::mutex> lk2(P2->mtex);
#endif
        //Adding forces and torques
        //P1->F -= Fnet;
        //P2->F += Fnet;
        Vec3_t FT = Fnet + Ftnet;
        F1 -= FT;
        F2 += FT;
        Vec3_t T, Tt;
        Tt = cross (x1,FT);
        Quaternion_t q1;
        Conjugate (P1->Q,q1);
        Rotation  (Tt,q1,T);
        //P1->T -= T;
        T1 -= T;
        Tt = cross (x2,FT);
        Quaternion_t q2;
        Conjugate (P2->Q,q2);
        Rotation  (Tt,q2,T);
        //P2->T += T;
        T2 += T;

        //Torque
        double arg = dot(p1,p2)/(norm(p1)*norm(p2));
        if (arg> 1.0) arg= 1.0;
        if (arg<-1.0) arg=-1.0;
        double Ant = acos(arg);
        if (dot(cross(p1,p2),n)<0) Ant = -Ant;
        Tt         = -Bm*(Ant-An)*n/L0;
        Rotation  (Tt,q1,T);
        //P1->T -= T;
        T1 -= T;
        Rotation  (Tt,q2,T);
        //P2->T += T;
        T2 += T;
#ifdef USE_THREAD
        //pthread_mutex_unlock(&lck);
        //pthread_mutex_unlock(&P1->lck);
        //pthread_mutex_unlock(&P2->lck);
#endif

        //Breaking point
        if ((fabs(delta)/(eta*eps)+norm(td)/(L0*eps)+0.0*fabs(Ant-An)*L0>1.0)&&(eps>=0.0))
        //double en = fabs(delta)/(eta*eps);
        //double et = norm(td)/(L0*eps);
        //if ((en*en + et*et > 1.0)&&(eps>=0.0))
        {
            valid = false;
        }
    }
    return false;
}

inline void BInteracton::UpdateParameters ()
{
    Bn              = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn)*Area;
    Bt              = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt)*Area;
    Bm              = 2*ReducedValue(P1->Props.Bm,P2->Props.Bm)*Area;
    double me       = ReducedValue(P1->Props.m ,P2->Props.m );
    Gn              = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn)*me;
    Gt              = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt)*me;
    eps             = 2*ReducedValue(P1->Props.eps,P2->Props.eps);
    eta             = 2*ReducedValue(P1->Props.Eta,P2->Props.Eta);
    I1              = P1->Index;
    I2              = P2->Index;
}

}
#endif //  MECHSYS_DEM_INTERACTON_H
