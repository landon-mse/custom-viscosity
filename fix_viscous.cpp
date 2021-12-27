/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_viscous.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"	
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixViscous::FixViscous(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  gamma(NULL)
{
  dynamic_group_allow = 1;

  if (narg < 4) error->all(FLERR,"Illegal fix viscous command");

  double gamma_one = force->numeric(FLERR,arg[3]);
  gamma = new double[atom->ntypes+1];
  for (int i = 1; i <= atom->ntypes; i++) gamma[i] = gamma_one;

  // optional args

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix viscous command");
      int itype = force->inumeric(FLERR,arg[iarg+1]);
      double scale = force->numeric(FLERR,arg[iarg+2]);
      if (itype <= 0 || itype > atom->ntypes)
        error->all(FLERR,"Illegal fix viscous command");
      gamma[itype] = gamma_one * scale;
      iarg += 3;
    } else error->all(FLERR,"Illegal fix viscous command");
  }

  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

FixViscous::~FixViscous()
{
  delete [] gamma;
}

/* ---------------------------------------------------------------------- */

int FixViscous::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixViscous::init()
{
  int max_respa = 0;

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = max_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,max_respa);
  }
  
  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 0;
  
}

/* ---------------------------------------------------------------------- */

void FixViscous::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixViscous::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixViscous::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscous::post_force(int /*vflag*/)
{
  // apply drag force to atoms in group
  // direction is opposed to velocity vector
  // magnitude depends on atom type
	
  double rpsq = 25;  //need to add as fix input eventually
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  double xreg[100], yreg[100], zreg[100];
  double xx, yy, zz, xy, xz, yz, xsum, ysum, zsum, vdotn;
  double A[3][3], b[3], n[3], nmag;
  int inum, jnum, i, ii, j, jj, nreg, k;
	
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  double shift = 1.03040961; //arbitrary shift to prevent ill-conditioning

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  // invoke full neighbor list (will copy or build if necessary)

  //neighbor->build_one(list);
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double drag;

	
  // determine the local normal direction
  // based on all neighbors within rp of
  // the atom, then impose the drag force in
  // that direction
  
  for (ii = 0; ii < inum; ii++) {
	  i = ilist[ii];
	  nreg = 0;
	
    if (mask[i] & groupbit) {
	    xtmp = x[i][0];
	    ytmp = x[i][1];
	    ztmp = x[i][2];
	    jlist = firstneigh[i];
	    jnum = numneigh[i];
		
		//add point to regression array
		xreg[nreg] = xtmp + shift;
		yreg[nreg] = ytmp + shift;
		zreg[nreg] = ztmp + shift;
		nreg += 1;
		
		
		
		 for (jj = 0; jj < jnum; jj++) {
	    	j = ilist[jj];
	    	if (mask[j] & groupbit) {
		      delx = xtmp - x[j][0];
		      dely = ytmp - x[j][1];
		      delz = ztmp - x[j][2];
		      domain -> minimum_image(delx,dely,delz);
			  rsq = delx*delx + dely*dely + delz*delz;

	    		if (rsq < rpsq) {
					xreg[nreg] = x[j][0] + shift;
					yreg[nreg] = x[j][1] + shift;
					zreg[nreg] = x[j][2] + shift;
					nreg += 1;
				}
	        }  
	    }
	  }
	  
	  
	  //perform the regression to identify the local normal n
	  xx = 0.0; yy = 0.0; zz = 0.0;
	  xy = 0.0; xz = 0.0; yz = 0.0;
	  xsum = 0.0; ysum = 0.0; zsum = 0.0;
	  for (k = 0; k < nreg; k++) {
		  xx += xreg[k]*xreg[k];
		  xy += xreg[k]*yreg[k];
		  xz += xreg[k]*zreg[k];
		  yy += yreg[k]*yreg[k];
		  yz += yreg[k]*zreg[k];
		  zz += zreg[k]*zreg[k];
		  xsum += xreg[k];
		  ysum += yreg[k];
		  zsum += zreg[k];
	  }
	  A[0][0] = xx; A[0][1] = xy; A[0][2] = xz;
	  A[1][0] = xy; A[1][1] = yy; A[1][2] = yz;
	  A[2][0] = xz; A[2][1] = yz; A[2][2] = zz;
	  b[0] = xsum; b[1] = ysum; b[2] = zsum;
	  
	  int ierror = MathExtra::mldivide3(A, b, n);
	  if (ierror==1)
 	  {
	    nmag=0;
	    n[0]=0;
	    n[1]=0;
	    n[2]=0;
            vdotn=0;
	  }
	  else
	  {
	  nmag = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	  n[0] /= nmag; n[1] /= nmag; n[2] /= nmag;
	  
	  vdotn = v[i][0]*n[0]+v[i][1]*n[1]+v[i][2]*n[2];
	  }
	  
      drag = gamma[type[i]];
      f[i][0] -= drag*vdotn*n[0];
      f[i][1] -= drag*vdotn*n[1];
      f[i][2] -= drag*vdotn*n[2];
  }
  
 
  //Original fix code  
  //for (int i = 0; i < nlocal; i++) {
//	if (mask[i] & groupbit) {
//	  drag = gamma[type[i]];
//	  f[i][0] -= drag*v[i][0];
//	  f[i][1] -= drag*v[i][1];
//	  f[i][2] -= drag*v[i][2];
//	}
//  }
    
	
}

/* ---------------------------------------------------------------------- */

void FixViscous::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscous::min_post_force(int vflag)
{
  post_force(vflag);
}
