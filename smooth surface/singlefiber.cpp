/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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
// Stokes law


//STD
#include<iostream>

// MechSys
#include "Domain_qian.h"
#include "InputData.h"      ///<data transfer by **************************qian****************************************

using namespace std;

InputData indat;               ///<data transfer by **************************qian****************************************

struct UserData
{
    std::ofstream      oss_ss;       ///< file for particle data
    Vec3_t                acc;
    double                 nu;
    double                  R;
    Array<Cell *>        xmin;
    Array<Cell *>        xmax;
	Vec3_t g;
};

  
   void AddCub3 (LBM::Domain & dom,int Tag, Vec3_t const & X, double R, double Lx, double Ly, double Lz,double rho, double Angle, Vec3_t * Axis)
  {
	  // vertices
	  Array<Vec3_t> V(8);
	  double lx = Lx/2.0;
	  double ly = Ly/2.0;
	  double lz = Lz/2.0;
	  
	  V[0] = -lx-1, -ly, -lz-1;
	  V[1] =  lx+1, -ly, -lz-1;
	  V[2] =  lx,  ly, -lz;
	  V[3] = -lx,  ly+1.5, -lz;
	  V[4] = -lx-1, -ly,  lz;
	  V[5] =  lx+1, -ly,  lz;
	  V[6] =  lx,  ly,	lz;
	  V[7] = -lx,  ly+1.5,	lz;
  
	  // edges
	  Array<Array <int> > E(12);
	  for (size_t i=0; i<12; ++i) E[i].Resize(2);
	  E[ 0] = 0, 1;
	  E[ 1] = 1, 2;
	  E[ 2] = 2, 3;
	  E[ 3] = 3, 0;
	  E[ 4] = 4, 5;
	  E[ 5] = 5, 6;
	  E[ 6] = 6, 7;
	  E[ 7] = 7, 4;
	  E[ 8] = 0, 4;
	  E[ 9] = 1, 5;
	  E[10] = 2, 6;
	  E[11] = 3, 7;
  
	  // faces
	  Array<Array <int> > F(6);
	  for (size_t i=0; i<6; i++) F[i].Resize(4);
	  F[0] = 4, 7, 3, 0;
	  F[1] = 1, 2, 6, 5;
	  F[2] = 0, 1, 5, 4;
	  F[3] = 2, 3, 7, 6;
	  F[4] = 0, 3, 2, 1;
	  F[5] = 4, 5, 6, 7;
  
	    double vol; // volume of the polyhedron
		  Vec3_t CM;  // Center of mass of the polyhedron
		  Mat3_t It;  // Inertia tensor of the polyhedron
		  DEM::PolyhedraMP(V,F,vol,CM,It);	 
		 DEM::Erosion(V,E,F,R);
 
			  // calculate the rotation
		  bool ThereisanAxis = true;
		  if (Axis==NULL)
		  {
			  Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
			  Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
			  ThereisanAxis = false;
		  }
		  Quaternion_t q;
		  NormalizeRotation (Angle,(*Axis),q);
		  for (size_t i=0; i<V.Size(); i++)
		  {
			  Vec3_t t;
			  Rotation (V[i],q,t);
			  V[i] = t;
		  }
 
	  
		  dom.Particles.Push(new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
		  dom.Particles[dom.Particles.Size()-1]->x			= CM;
		  dom.Particles[dom.Particles.Size()-1]->Props.V	= vol;
		  std::cout << "V="<< dom.Particles[dom.Particles.Size()-1]->Props.V<<std::endl;
		  dom.Particles[dom.Particles.Size()-1]->Props.m	= rho*dom.Particles[dom.Particles.Size()-1]->Props.V;
 
		   
		  Vec3_t I;
		  Quaternion_t Q;
		  Vec3_t xp,yp,zp;
		  
		  Eig(It,I,xp,yp,zp);
		  
		  CheckDestroGiro(xp,yp,zp);
		  
		  I *= rho;
		  Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
		  Q(1) = (yp(2)-zp(1))/(4*Q(0));
		  Q(2) = (zp(0)-xp(2))/(4*Q(0));
		  Q(3) = (xp(1)-yp(0))/(4*Q(0));
		  Q = Q/norm(Q);
  
		  double Dmax = DEM::Distance(CM,V[0])+R;
		  for (size_t i=1; i<V.Size(); ++i)
		   {
			  if (DEM::Distance(CM,V[i])+R > Dmax) Dmax = DEM::Distance(CM,V[i])+R;
			  }
		  
		  dom.Particles[dom.Particles.Size()-1]->Q			= Q;
			  
		  dom.Particles[dom.Particles.Size()-1]->I			= I;
  
		  dom.Particles[dom.Particles.Size()-1]->Ekin		= 0.0;
  
		  dom.Particles[dom.Particles.Size()-1]->Erot	 = 0.0;
  
		  dom.Particles[dom.Particles.Size()-1]->Dmax		= Dmax;
  
		  dom.Particles[dom.Particles.Size()-1]->PropsReady = true;
		  
		  dom.Particles[dom.Particles.Size()-1]->Index		= dom.Particles.Size()-1;
		  
		  dom.Particles[dom.Particles.Size()-1]->Position(X);
	  
  
  }


  void AddCub2 (LBM::Domain & dom,int Tag, Vec3_t const & X, double R, double Lx, double Ly, double Lz,double rho, double Angle, Vec3_t * Axis)
 {
	 // vertices
	 Array<Vec3_t> V(8);
	 double lx = Lx/2.0;
	 double ly = Ly/2.0;
	 double lz = Lz/2.0;
	 
	 V[0] = -lx, -ly, -lz;
	 V[1] =  lx, -ly, -lz;
	 V[2] =  lx,  ly, -lz;
	 V[3] = -lx,  ly, -lz;
	 V[4] = -lx, -ly,  lz;
	 V[5] =  lx, -ly,  lz;
	 V[6] =  lx,  ly,  lz;
	 V[7] = -lx,  ly,  lz;
 
	 // edges
	 Array<Array <int> > E(12);
	 for (size_t i=0; i<12; ++i) E[i].Resize(2);
	 E[ 0] = 0, 1;
	 E[ 1] = 1, 2;
	 E[ 2] = 2, 3;
	 E[ 3] = 3, 0;
	 E[ 4] = 4, 5;
	 E[ 5] = 5, 6;
	 E[ 6] = 6, 7;
	 E[ 7] = 7, 4;
	 E[ 8] = 0, 4;
	 E[ 9] = 1, 5;
	 E[10] = 2, 6;
	 E[11] = 3, 7;
 
	 // faces
	 Array<Array <int> > F(6);
	 for (size_t i=0; i<6; i++) F[i].Resize(4);
	 F[0] = 4, 7, 3, 0;
	 F[1] = 1, 2, 6, 5;
	 F[2] = 0, 1, 5, 4;
	 F[3] = 2, 3, 7, 6;
	 F[4] = 0, 3, 2, 1;
	 F[5] = 4, 5, 6, 7;
 
	 // calculate the rotation
	 bool ThereisanAxis = true;
	 if (Axis==NULL)
	 {
		 Angle	 = (1.0*rand())/RAND_MAX*2*M_PI;
		 Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
		 ThereisanAxis = false;
	 }
	 Quaternion_t q;
	 NormalizeRotation (Angle,(*Axis),q);
	 for (size_t i=0; i<V.Size(); i++)
	 {
		 Vec3_t t;
		 Rotation (V[i],q,t);
		 V[i] = t+X;
	 }
 
	 // add particle
	  dom.Particles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
 
	 // clean up
	 if (!ThereisanAxis) delete Axis;
	 q(0) = 1.0;
	 q(1) = 0.0;
	 q(2) = 0.0;
	 q(3) = 0.0;
	 q = q/norm(q);
 
	  dom.Particles[ dom.Particles.Size()-1]->Q		   = q;
	 dom. Particles[ dom.Particles.Size()-1]->Props.V    = Lx*Ly*Lz;
	 dom. Particles[ dom.Particles.Size()-1]->Props.m    = rho*Lx*Ly*Lz;
	 dom. Particles[ dom.Particles.Size()-1]->I		   = Lx*Lx, Ly*Ly, Lz*Lz;
	 dom. Particles[ dom.Particles.Size()-1]->I		  *=  dom.Particles[ dom.Particles.Size()-1]->Props.m/6.0;
	 dom. Particles[ dom.Particles.Size()-1]->x		   = X;
	 dom. Particles[ dom.Particles.Size()-1]->Ekin	   = 0.0;
	 dom. Particles[ dom.Particles.Size()-1]->Erot	   = 0.0;
	 dom. Particles[ dom.Particles.Size()-1]->Dmax	   = sqrt(Lx*Lx+Ly*Ly+Lz*Lz)/2+R;
	 dom. Particles[ dom.Particles.Size()-1]->PropsReady = true;
	 dom. Particles[ dom.Particles.Size()-1]->Index	   =  dom.Particles.Size()-1;
	 
 
 }

 
 void AddRiceR (LBM::Domain & dom, int Tag, const Vec3_t & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
 {
	 // vertices
	 Array<Vec3_t> V(2);
	 V[0] = 0.0, 0.0,  L/2;
	 V[1] = 0.0, 0.0, -L/2;
 
	 // edges
	 Array<Array <int> > E(1);
	 E[0].Resize(2);
	 E[0] = 0, 1;
 
	 // faces
	 Array<Array <int> > F(0); // no faces
 
	 // calculate the rotation
	 bool ThereisanAxis = true;
	 if (Axis==NULL)
	 {
		 Angle	 = (1.0*rand())/RAND_MAX*2*M_PI;
		 Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
		 ThereisanAxis = false;
	 }
	 Quaternion_t q;
	 NormalizeRotation (Angle,(*Axis),q);
	 for (size_t i=0; i<V.Size(); i++)
	 {
		 Vec3_t t;
		 Rotation (V[i],q,t);
		 V[i] = t+X;
	 }
 
	 // add particle
	 dom.Particles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
 
	 dom.Particles[dom.Particles.Size()-1]->Q		   = q;
	 dom.Particles[dom.Particles.Size()-1]->Props.V    = (4.0/3.0)*M_PI*R*R*R + M_PI*L*R*R;
	 dom.Particles[dom.Particles.Size()-1]->Props.m    = rho*dom.Particles[dom.Particles.Size()-1]->Props.V;
	 dom.Particles[dom.Particles.Size()-1]->I		   = (1.0/3.0)*M_PI*rho*R*R*R*L*L + (1.0/12.0)*M_PI*rho*R*R*L*L*L + (3.0/4.0)*M_PI*rho*R*R*R*R*L + (8.0/15.0)*M_PI*rho*R*R*R*R*R,
												 (1.0/3.0)*M_PI*rho*R*R*R*L*L + (1.0/12.0)*M_PI*rho*R*R*L*L*L + (3.0/4.0)*M_PI*rho*R*R*R*R*L + (8.0/15.0)*M_PI*rho*R*R*R*R*R,
												 0.5*M_PI*rho*R*R*R*R*L + (8.0/15.0)*M_PI*rho*R*R*R*R*R;
	 dom.Particles[dom.Particles.Size()-1]->x		   = X;
	 dom.Particles[dom.Particles.Size()-1]->Ekin	   = 0.0;
	 dom.Particles[dom.Particles.Size()-1]->Erot	   = 0.0;
	 dom.Particles[dom.Particles.Size()-1]->Dmax	   = 0.5*L + R;
	 dom.Particles[dom.Particles.Size()-1]->PropsReady = true;
	 dom.Particles[dom.Particles.Size()-1]->Index	   = dom.Particles.Size()-1;
 
	 // clean up
	 if (!ThereisanAxis) delete Axis;
 }



 void AddicosahedronParticle(LBM::Domain & dom, int Tag, Vec3_t const & X, double size,double R,   double rho, double Angle, Vec3_t * Axis)
  {
	  
	 
	  // vertices
		  Array<Vec3_t> V(12);
 
		  double m,n;
		  m=sqrt(50-10*sqrt(5))/10*size;
		  n=sqrt(50+10*sqrt(5))/10*size;
 
 
		  V[0]=m,0,n;
		  V[1]=m,0,-n;
		  V[2]=-m,0,n;
		  V[3]=-m,0,-n;
		  
		  V[4]=0,n,m;
		  V[5]=0,-n,m;
		  V[6]=0,n,-m;
		  V[7]=0,-n,-m;
		  
		  V[8]=n,m,0;
		  V[9]=-n,m,0;
		  V[10]=n,-m,0;
		  V[11]=-n,-m,0;
 
		 
		  // edges
		  Array<Array <int> > E(30);
		  
		  for (size_t i=0; i<30; ++i) E[i].Resize(2);
 
		  E[0] = 6, 4;//
		  E[1] = 6, 8;//
		  E[2] = 4, 8;//
		  E[3] = 9, 6;//
		  E[4] = 6, 1;//
		  E[5] = 6, 3;//
 
		  E[6] = 1, 3;//
		  E[7] = 8, 1;//
		  E[8] = 10, 1;//
		  E[9] = 8, 0;//
		  E[10] = 8, 10;//
 
		  E[11] = 0, 10;//
		  E[12] = 4, 0;//
		  E[13] = 2, 0;//
		  E[14] = 4, 9;//
		  E[15] = 4, 2;//
 
		  E[16] = 9, 2;//
		  E[17] = 9, 3;//
		  E[18] = 9, 11;//
		  E[19] = 3, 11;//
		  E[20] = 3, 7;////////
 
		  E[21] = 1, 7;////////
		  E[22] = 0, 5;//
		  E[23] = 2, 11;//
		  E[24] = 2, 5;//
		  E[25] = 11, 5;//
 
		  E[26] = 7, 11;//
		  E[27] = 10, 5;//
		  E[28] = 10, 7;//
		  E[29] = 5, 7;//////

		  
 
	  
		  // faces
		  Array<Array <int> > F(20);
		  
		  for (size_t i=0; i<20; i++) F[i].Resize(3);
 
		  F[0]=4,6,8;
		  F[1]=4,9,6;
		  F[2]=6,9,3;
		  F[3]=6,1,3;
		  F[4]=6,3,1;
		  
		  F[5]=8,1,10;
		  F[6]=8,10,1;
		  F[7]=8,0,4;
		  F[8]=4,0,2;
		  F[9]=4,2,9;
		  
		  F[10]=9,2,11;
		  F[11]=9,11,3;
		  F[12]=3,7,1;
		  F[13]=1,7,10;
		  F[14]=10,5,0;
 
		  F[15]=0,5,2;
		  F[16]=2,5,11;
		  F[17]=3,11,7;
		  F[18]=5,7,11;
		  F[19]=10,7,5;
		 
	 
		  
		  double vol; // volume of the polyhedron
		  Vec3_t CM;  // Center of mass of the polyhedron
		  Mat3_t It;  // Inertia tensor of the polyhedron
		  DEM::PolyhedraMP(V,F,vol,CM,It);	 
		 DEM::Erosion(V,E,F,R);
 
			  // calculate the rotation
		  bool ThereisanAxis = true;
		  if (Axis==NULL)
		  {
			  Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
			  Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
			  ThereisanAxis = false;
		  }
		  Quaternion_t q;
		  NormalizeRotation (Angle,(*Axis),q);
		  for (size_t i=0; i<V.Size(); i++)
		  {
			  Vec3_t t;
			  Rotation (V[i],q,t);
			  V[i] = t;
		  }
 
	  
		  dom.Particles.Push(new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
		  dom.Particles[dom.Particles.Size()-1]->x			= CM;
		  dom.Particles[dom.Particles.Size()-1]->Props.V	= vol;
		  std::cout << "V="<< dom.Particles[dom.Particles.Size()-1]->Props.V<<std::endl;
		  dom.Particles[dom.Particles.Size()-1]->Props.m	= rho*dom.Particles[dom.Particles.Size()-1]->Props.V;
 
		   
		  Vec3_t I;
		  Quaternion_t Q;
		  Vec3_t xp,yp,zp;
		  
		  Eig(It,I,xp,yp,zp);
		  
		  CheckDestroGiro(xp,yp,zp);
		  
		  I *= rho;
		  Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
		  Q(1) = (yp(2)-zp(1))/(4*Q(0));
		  Q(2) = (zp(0)-xp(2))/(4*Q(0));
		  Q(3) = (xp(1)-yp(0))/(4*Q(0));
		  Q = Q/norm(Q);
  
		  double Dmax = DEM::Distance(CM,V[0])+R;
		  for (size_t i=1; i<V.Size(); ++i)
		   {
			  if (DEM::Distance(CM,V[i])+R > Dmax) Dmax = DEM::Distance(CM,V[i])+R;
			  }
		  
		  dom.Particles[dom.Particles.Size()-1]->Q			= Q;
			  
		  dom.Particles[dom.Particles.Size()-1]->I			= I;
  
		  dom.Particles[dom.Particles.Size()-1]->Ekin		= 0.0;
  
		  dom.Particles[dom.Particles.Size()-1]->Erot	 = 0.0;
  
		  dom.Particles[dom.Particles.Size()-1]->Dmax		= Dmax;
  
		  dom.Particles[dom.Particles.Size()-1]->PropsReady = true;
		  
		  dom.Particles[dom.Particles.Size()-1]->Index		= dom.Particles.Size()-1;
		  
		  dom.Particles[dom.Particles.Size()-1]->Position(X);
  
  }



 void AddIrregularParticle(LBM::Domain & dom, int Tag, Vec3_t const & X, double Cylinder_radius,double R, double L, double nodes_number, double rho, double Angle, Vec3_t * Axis)
 {
	 
	
	 // vertices
		 Array<Vec3_t> V(2*nodes_number);


		 for(int i=0;i<nodes_number;i++)
		 	{
			V[i] = Cylinder_radius*cos((360.0/nodes_number)*(M_PI/180)*i), Cylinder_radius*sin((360.0/nodes_number)*(M_PI/180)*i), L/2;
		 }

		 for(int i=nodes_number;i<2*nodes_number;i++)
					 {
					 V[i] = Cylinder_radius*cos((360.0/nodes_number)*(M_PI/180)*i), Cylinder_radius*sin((360.0/nodes_number)*(M_PI/180)*i), -L/2;
				  }		 
		 
		
		 // edges
		 Array<Array <int> > E(3*nodes_number);
		 for (size_t i=0; i<3*nodes_number; ++i) E[i].Resize(2);

		 for(int i=0;i<nodes_number-1;i++)   E[ i ] = i, i+1;
		 E[ nodes_number-1 ] = nodes_number-1, 0;
		 
		 for(int i=nodes_number;i<2*nodes_number;i++)   E[ i ] = i-nodes_number, i;
		 
		 for(int i=2*nodes_number;i<3*nodes_number-1;i++)  
		 	{

			int j=i-nodes_number;				
		    E[ i ] = j, j+1;
      
		 	}
					 E[ 3*nodes_number-1 ] = 2*nodes_number-1, nodes_number;	

	
	 
		 // faces
		 Array<Array <int> > F(nodes_number+2);
		 
		 for (size_t i=0; i<nodes_number; i++) F[i].Resize(4);
		 for(int i=0;i<nodes_number-1;i++)   
		 	{
			 F[i] = i, i+nodes_number, i+1+nodes_number, i+1;			
		 }
		 
		 F[nodes_number-1] = nodes_number-1, nodes_number-1+nodes_number, nodes_number, 0;

		
		 for (size_t i=nodes_number; i<nodes_number+2; i++) F[i].Resize(nodes_number);
		 for(int j=0;j<nodes_number;j++)  
		 	{

		 F[nodes_number][j] = j;
		 F[nodes_number+1][j] = 2*nodes_number-j-1;	 
		 
		 	}    
	
		 
		 double vol; // volume of the polyhedron
		 Vec3_t CM;  // Center of mass of the polyhedron
		 Mat3_t It;  // Inertia tensor of the polyhedron
		 DEM::PolyhedraMP(V,F,vol,CM,It);	 
		DEM::Erosion(V,E,F,R);

		
			 // calculate the rotation
		 bool ThereisanAxis = true;
		 if (Axis==NULL)
		 {
			 Angle	 = (1.0*rand())/RAND_MAX*2*M_PI;
			 Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
			 ThereisanAxis = false;
		 }
		 Quaternion_t q;
		 NormalizeRotation (Angle,(*Axis),q);
		 for (size_t i=0; i<V.Size(); i++)
		 {
			 Vec3_t t;
			 Rotation (V[i],q,t);
			 V[i] = t;
		 }


	 
		 dom.Particles.Push(new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
		 dom.Particles[dom.Particles.Size()-1]->x		   = CM;
		 dom.Particles[dom.Particles.Size()-1]->Props.V    = vol;
		 dom.Particles[dom.Particles.Size()-1]->Props.m    = rho*dom.Particles[dom.Particles.Size()-1]->Props.V;

		
		  
		 Vec3_t I;
		 Quaternion_t Q;
		 Vec3_t xp,yp,zp;
		 
		 Eig(It,I,xp,yp,zp);
		 
		 CheckDestroGiro(xp,yp,zp);
		 
		 I *= rho;
		 
		 Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
		 Q(1) = (yp(2)-zp(1))/(4*Q(0));
		 Q(2) = (zp(0)-xp(2))/(4*Q(0));
		 Q(3) = (xp(1)-yp(0))/(4*Q(0));
		 Q = Q/norm(Q);
 
		 double Dmax = DEM::Distance(CM,V[0])+R;
		 for (size_t i=1; i<V.Size(); ++i)
		  {
			 if (DEM::Distance(CM,V[i])+R > Dmax) Dmax = DEM::Distance(CM,V[i])+R;
			 }
		 
		 dom.Particles[dom.Particles.Size()-1]->Q		   = Q;
			 
		 dom.Particles[dom.Particles.Size()-1]->I		   = I;
 
		 dom.Particles[dom.Particles.Size()-1]->Ekin	   = 0.0;
 
		 dom.Particles[dom.Particles.Size()-1]->Erot	= 0.0;
 
		 dom.Particles[dom.Particles.Size()-1]->Dmax	   = Dmax;
 
		 dom.Particles[dom.Particles.Size()-1]->PropsReady = true;
		 
		 dom.Particles[dom.Particles.Size()-1]->Index	   = dom.Particles.Size()-1;
		 
		 dom.Particles[dom.Particles.Size()-1]->Position(X);
 
	 
 }


 
void Setup(LBM::Domain & dom, void * UD)
{


// Boundary conditon
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
	for (size_t i=0; i<indat.Left.Size(); i++) // velocity inlet,Non-equilibrium extrapolation
	{


	//std::cout<<"left size"<<indat.Left.Size()<<std::endl;
		Cell * c = indat.Left[i];
        Cell * cR= indat.LeftR[i];
		if (c->IsSolid) continue;
		c->Vel = indat.uf_in;
		c->Rho=cR->Rho;
        for (size_t k=0; k< c->Nneigh; k++)
        {
            c->F[k] = c->Feq(k,indat.uf_in,c->Rho)+cR->F[k]-cR->Feq(k,cR->Vel,cR->Rho);
	    }
    }
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
    for (size_t i=0; i<indat.Right.Size(); i++) // fully developed,Non-equilibrium extrapolation
	{
		Cell * c = indat.Right[i];
        Cell * cL= indat.RightL[i];
		if (c->IsSolid) continue;
		c->Vel 	 = cL->Vel;
		 c->Rho   =cL->Rho;
        for (size_t k=0; k<c->Nneigh;k++)
        {
           c->F[k] = c->Feq(k,c->Vel,c->Rho)+ cL->F[k]-cL->Feq(k,cL->Vel, cL->Rho);
	    }
    }


	

}

void Report(LBM::Domain & dom, void * UD)
{


}

int main(int argc, char **argv) try
{
  
    size_t Nproc = 16; 
    if (argc==3) Nproc=atoi(argv[2]);
	
    bool   Render = true;
     size_t nx =180;
    size_t ny =180;
    size_t nz =220;
	
    double nu = 0.14855;
    double dx = 1.0;
    double dt = 1.0;//2.5e-10
    double Dp = 0.0;
    double R  = 4.5;

	double rp    = 4.5;   //jubenyixikeli
    double rp2    = 3.0;    //eryanghuagui
	double rp3    = 2.5;    //yanghualv

	double        E=2.153e9;   //PA6
	double        Sig=0.33;
	double        G=E/(2*(1+Sig));
	double        Gama=400;

    double        E2=2.4e9;// YAKELIBAN
	double        Sig2=0.33;
	double        G2=E2/(2*(1+Sig2));
	double        Gama2=700;

    double		  E3=2.5e8;    //yanghualv
	double		  Sig3=0.25;
	double		  G3=E3/(2*(1+Sig3));
	double		  Gama3=625;
	
    double w  = 0.00;
    double ang= 0.0;
    double Tf = 5.0e4;
	double rho_g = 1.205;
	Vec3_t uf_in = {0.000755767, 0.0, -0.47236};
	Vec3_t uf_in2 = {0.0, 0.0, 0.0};

	size_t inter_time=5; 
  



   indat.nx 	= nx;
   indat.ny 	= ny;
   indat.nz 	= nz;

   indat.dt	 = dt;
   indat.dx	 = dx;
	
   indat.Length  = dx*nx;
   indat.Width	= dx*ny;
   indat.Depth	= dx*nz;
   
   indat.nu 	= nu;
   indat.rp 	= rp;		 //sphere particle radius
   indat.rp2 	= rp2;		 //sphere particle radius
   indat.rp3 	= rp3;		 //sphere particle radius
   indat.dp 	= R;

   indat.rho_p	  = 1156.; 
	 indat.rho_p2 = 1190.; 
	 indat.rho_p3 = 3940.; 

   indat.inter_time  =inter_time;	//sphere inter every inter_sphere time
   indat.inter_number  =1;

   indat.E    =E;
   indat.Sig  =Sig;
   indat.G    =G;
   indat.Gama =Gama;
   
   indat.E2    =E2;
   indat.Sig2  =Sig2;
   indat.G2    =G2;
   indat.Gama2 =Gama2;

   indat.E3    =E3;
   indat.Sig3  =Sig3;
   indat.G3    =G3;
   indat.Gama3 =Gama3;


   indat.rho_g	= rho_g;
   indat.uf_in=0.0;
   indat.uf_in_particle=0.0;
 
	indat.Np_tol = 0;
	indat.Np_tol_in=0;
        indat.Np_tol2 = 0;
	indat.Np_tol_in2=0;
	indat.Np_tol_in3=0;
	indat.Np_cap = 0;
	indat.Np_ecp = 0;



    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);


   Dom.initInputData (&indat);	  ///<data transfer by **************************qian****************************************
    Dom.Step     = 1;
    Dom.Sc       = 0.0;
	
    UserData dat;
    Dom.UserData = &dat;
    dat.acc      = Vec3_t(0.0,0.0,0.0);
    dat.R        = R;
    dat.nu       = nu;
	dat.g		 = 0.0,0.0,-9.8e-4;



	
	//Assigning boundaries at left and right (inlet and outlet)
	  for (size_t i=0;i<ny;i++)
	  	  for (size_t j=0;j<nz;j++)
	  {
		  indat.Left  .Push(Dom.Lat[0].GetCell(iVec3_t(0   ,i,j)));
		  indat.LeftR .Push(Dom.Lat[0].GetCell(iVec3_t(1   ,i,j)));
		  
		  indat.RightL.Push(Dom.Lat[0].GetCell(iVec3_t(nx-2,i,j)));
		  indat.Right .Push(Dom.Lat[0].GetCell(iVec3_t(nx-1,i,j)));
		  
		  indat.Left[i]->RhoBC = indat.rho_g;
		  indat.Right[i]->RhoBC = indat.rho_g;
	  }



    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        Vec3_t v0(0.0,0.0,0.0);
        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(indat.rho_g,v0);
    }



/*               
	//Assigning solid boundaries at top and bottom
	   for (size_t i=0;i<nx;i++)
	   	  for (size_t j=0;j<ny;j++)
	   {
		   Dom.Lat[0].GetCell(iVec3_t(i,j,0))->IsSolid = true;
		   Dom.Lat[0].GetCell(iVec3_t(i,j,nz-1))->IsSolid = true;
		   
	   }
	   
	   //Assigning solid boundaries at left and right

		   for (size_t i=0;i<nx;i++)
		   	  for (size_t k=0;k<nz;k++)
		   {
			   Dom.Lat[0].GetCell(iVec3_t(i,0,k))->IsSolid = true;
			   Dom.Lat[0].GetCell(iVec3_t(i,ny-1,k))->IsSolid = true;
		   }
*/




			  
			

				//AddRice (LBM::Domain & dom, int Tag, const Vec3_t & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
/*
	         Dom.AddSphere(-1,Vec3_t(0.3* nx*dx,0.5* ny*dx,0.5* nz*dx),3,indat.rho_p);
			 Dom.Particles[0]->v =uf_in ;	  //*************************************************************************fiber velocity add by qian
			 Dom.Particles[0]->Props.particle_radius=rp ;
			Dom.Particles[0]->Props.E =E ;
			Dom.Particles[0]->Props.Sig =Sig ;
			Dom.Particles[0]->Props.G =G ;
			Dom.Particles[0]->Props.Gama =Gama ;
                       */

	
			// AddicosahedronParticle(Dom, -1, Vec3_t(0.5* nx*dx,0.5* ny*dx,0.8* nz*dx), 5.3, 0.5,	indat.rho_p,0.0,&OrthoSys::e1); 	 
			 AddCub3 (Dom,-1, Vec3_t(0.5* nx*dx,0.5* ny*dx, 55.0094), 0.5, 16.6, 16.9, 19.7, indat.rho_p, 6.0*M_PI/180.0,&OrthoSys::e0);	

	       //  Dom.AddSphere(-1,Vec3_t(0.5* nx*dx,0.5* ny*dx,0.8* nz*dx),rp,indat.rho_p);
			 Dom.Particles[0]->v =uf_in ;	  //*************************************************************************fiber velocity add by qian
			Dom.Particles[0]->Props.particle_radius= Dom. Particles[0]->Dmax;   
			Dom.Particles[0]->Props.E =E ;
			Dom.Particles[0]->Props.Sig =Sig ;
			Dom.Particles[0]->Props.G =G ;
			Dom.Particles[0]->Props.Gama =Gama ;

			
            Dom.Particles[0]->Ff = Dom.Particles[0]->Props.m*dat.g;    

				
				//void AddPlane		 (int Tag, Vec3_t const & X, double R, double Lx,double Ly , double rho, double Angle=0, Vec3_t * Axis=NULL);					 ///< Add a cube at position X with spheroradius R, side of length L and density rho
			//	Dom.AddPlane		 (-1, Vec3_t(0.5* nx*dx,0.5* ny*dx,0.3* nz*dx), 5, 60,60 , indat.rho_p, 0.0, &OrthoSys::e1);					 ///< Add a cube at position X with spheroradius R, side of length L and density rho


		   AddCub2 (Dom,-1, Vec3_t(0.5* nx*dx,0.5* ny*dx,20), 0.5, nx, ny, 20,indat.rho_p2, 0.0, &OrthoSys::e1);
		   
		   
		   Dom.Particles[1]->Props.particle_radius=Dom. Particles[1]->Dmax;   
		   Dom.Particles[1]->Props.E =E2 ;
		   Dom.Particles[1]->Props.Sig =Sig2 ;
		   Dom.Particles[1]->Props.G =G2 ;
		   Dom.Particles[1]->Props.Gama =Gama2 ;
		   
		   Dom.Particles[1]->FixVeloc();


	

		
	  /* for (size_t i=0;i<Dom.Particles.Size();i++)
	   {
		   Dom.Particles[i]->Ff = Dom.Particles[i]->Props.m*dat.g;
		  
	   }*/
	

    //Solving
    Dom.Solve(Tf,10,Setup,Report,"singerfiber",3,Nproc);
}
MECHSYS_CATCH
