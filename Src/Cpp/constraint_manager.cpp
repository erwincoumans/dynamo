/*
  DYNAMO - Dynamic Motion library
  Copyright (C) 1996-1999 Bart Barenbrug

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Library General Public
  License as published by the Free Software Foundation; either
  version 2 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Library General Public License for more details.

  You should have received a copy of the GNU Library General Public
  License along with this library; if not, write to the Free
  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  Please send remarks, questions and bug reports to bartb@win.tue.nl,
  or write to:
                  Bart Barenbrug
		  Department of Mathematics and Computing Science
		  Eindhoven University of Technology
		  P.O. Box 513, 5600 MB Eindhoven, The Netherlands
*/

//
// filename     : constraint_manager.cpp
// description	: non-inline methods of class DL_constraint_manager
// author	: Bart Barenbrug     April 1996
//                Mischa Courtin     September 1996 (sparse dCdR)
//

#include "constraint_manager.h"
#include "dyna_system.h"
#include "euler.h"
#include "NaN.h"

//#define DCDR
//#define DEBUG

// pointer to the one and only constraint manager:
DL_constraint_manager* DL_constraints=NULL;

// ************************** //
// non-inline member fuctions //
// ************************** //

void DL_constraint_manager::calc_all_errors(DL_largevector *lv) {
  static DL_largevector err;
  DL_constraint *constr=(DL_constraint *)c->getfirst();
  while (constr) {
    err.resize(constr->dim);
    constr->get_error(&err);
    lv->setsubvector(constr->index,&err);
    constr=(DL_constraint*)c->getnext(constr);
  }
}

boolean DL_constraint_manager::apply_all_restriction_changes(DL_largevector *lv) {
  static DL_largevector restr;
  DL_constraint *next_c,*constr=(DL_constraint *)c->getfirst();
  // first see if no reactionforces become too large:
  while (constr) {
    next_c=(DL_constraint*)c->getnext(constr);
    restr.resize(constr->dim);
    lv->getsubvector(constr->index,&restr);
    constr->test_restriction_changes(&restr);
    constr=next_c;
  }
  if (c_changed) {
    //one or more constraints deleted themselves: recalculate dCdR:
    // first redo the index administration:
    redo_index_administration();
    
    // then recalculate dCdR:
    calc_dCdR_full();
    return TRUE;
  }
  else { // everything ok: really apply the changes:
    constr=(DL_constraint *)c->getfirst();
    while (constr) {
      restr.resize(constr->dim);
      lv->getsubvector(constr->index,&restr);
      constr->apply_restriction_changes(&restr);
      constr=(DL_constraint*)c->getnext(constr);
    }
    return FALSE;
  }
}

void DL_constraint_manager::add(DL_constraint *constr) {
  c->addelem(constr);
  if (show_con_forces) constr->show_forces();
  else constr->hide_forces();
  c_changed=TRUE;
}

void DL_constraint_manager::del(DL_constraint *constr) {
  constr->hide_forces();
  c->remelem(constr);
  c_changed=TRUE;
}

void DL_constraint_manager::satisfy() {
  static DL_largevector dC;
  int i, nr_collisionloops=0;

  if (c_changed) redo_index_administration();
  dC.resize(totdim);

  new_frame();
  
  calc_all_errors(&dC);
  first_error=dC.norm();

  apply_first_estimates();
  
  calc_all_errors(&dC);
  error=dC.norm();
  
  if ((error>4*first_error) || NaN(error)) {
#ifdef DEBUG    
  DL_dsystem->get_companion()->Msg("resetting at %d\n", DL_dsystem->frame_number() );
#endif
#undef DEBUG 
    reset_undo_all();
    calc_all_errors(&dC);
    error=dC.norm();  
  }

  nriter=0;

  do {  //  while ((nrcollisions>0) && (nr_collisionloops<max_collisionloops));
    nr_collisionloops++;
    // here is the point to insert collision detection: all geos have
    // been moved to their new location and to all dynas inertia and the
    // first estimates of reaction forces have been applied (so getpoint
    // will return the old position, and getnewpoint will return a good
    // estimate of the position at the next frame).

    nrcollisions=0;
    if (max_collisionloops>0) {
//      DL_dsystem->update_dyna_companions();
      DL_dsystem->get_companion()->do_collision_detection();
    }
    
    if (c->length()==0) return;  // no constraints to satisfy

    if (c_changed) {  // redo index administration:
      redo_index_administration();
      dC.resize(totdim);
      calc_all_errors(&dC);
      error=dC.norm();
    }
    else {
      if (nr_collisionloops>1) {
	dC.resize(totdim);
	calc_all_errors(&dC);
	error=dC.norm();
      }
    }
    if (error>max_error) {  // recalculate dCdR
      if (c_changed) { // have to rebuild dCdR from scratch
        calc_dCdR_full();
        dCdRToGo=NrSkip;
	// calc_dCdR_full might have permutated constraints:
	dC.resize(totdim);
	calc_all_errors(&dC);
	error=dC.norm();
      }
      else {
        if (dCdRToGo==0) { // rebuild dCdR using the info from cp.
          if (analytical) calc_dCdR_analytical();
          else calc_dCdR_empirical();
          dCdRToGo=NrSkip;
	}
        else dCdRToGo--;
      }
      iterate(&dC);
    }
    // delete collision constraints here
    for (i=0;i<nrcollisions;i++) collisions[i]->post_processing();
  }
  while ((nrcollisions>0) && (nr_collisionloops<max_collisionloops));

/*
// for testing of the collision constraint (floore.lks):
static int nr_coll_changes=0;
static DL_Scalar total_change=0;
if ((nrcollisions>0) || (nr_collisionloops>0)) {
  nr_coll_changes++;
  total_change+=fabs(DL_dsystem->totenergy()-DL_dsystem->newtotenergy());
  DL_dsystem->get_companion()->Msg("energychange #%d (time=%f): %f (avg so far: %f)\n",
                                   nr_coll_changes, DL_dsystem->time(),
				   fabs(DL_dsystem->totenergy()-DL_dsystem->newtotenergy()),
				   total_change/nr_coll_changes );
}
*/
  
  do_post_processing();
}

void DL_constraint_manager::iterate(DL_largevector *dC) {
  static DL_largevector dR;
  first_error=error;
  boolean singular=dCdR->prep_for_solve();
  dR.resize(dC->get_dim());
  while ((error>max_error) && (nriter<MaxIter)) {
    nriter++;
    dC->neg(dC);
    dCdR->solve(&dR,dC);
    if (apply_all_restriction_changes(&dR)) {
      dR.resize(totdim);
      dC->resize(totdim);
      calc_all_errors(dC);
      error=dC->norm();
    }
    else {
      calc_all_errors(dC);
      error=dC->norm();
      if ((error>4*first_error) || NaN(error)) {
	// clear divergence: try a more stable solve method
	switch (dCdR->get_solve_method()) {
	case lud_bcksub:
	  dCdR->set_solve_method(conjug_grad);
	  singular=dCdR->prep_for_solve();
	  break;
	case conjug_grad:
	  dCdR->set_solve_method(svd);
	  singular=dCdR->prep_for_solve();
	  break;
	case svd:
	  // nothing we can do...
	  nriter=MaxIter;
          DL_dsystem->get_companion()->Msg("Warning: Can not solve constraints at frame %d\n", DL_dsystem->frame_number() );
	  // possibly raise an event here
	  break;
	}
	dR.neg(&dR);
	apply_all_restriction_changes(&dR);
	calc_all_errors(dC);
	error=dC->norm();
      }
    }
  }
  switch (dCdR->get_solve_method()) {
  case lud_bcksub: break;
  case conjug_grad:
    if ((0<nriter) && (nriter<MaxIter)) {
      if ((!singular) && (error<first_error) && (nr_cg>10)) {
	nr_cg=0;
	dCdR->set_solve_method(lud_bcksub);
      }
      else nr_cg++;
    }
    else if (nriter>=MaxIter) nr_cg=0;
    break;
  case svd:
    if ((0<nriter) && (nriter<MaxIter)) {
      if ((!singular) && (error<first_error) && (nr_cg>10)) {
	nr_cg=0;
	dCdR->set_solve_method(lud_bcksub);
      }
      else nr_cg++;
    }
    else if (nriter>=MaxIter) nr_cg=0;
    break;
  }
}

void DL_constraint_manager::calc_dCdR_empirical() {
  static DL_largevector org_err;
  static DL_largevector new_err;
  static DL_largevector err_dif;
  static DL_largevector test_restr;
  int i;
  
  org_err.resize(totdim);
  new_err.resize(totdim);
  err_dif.resize(totdim);

  // use an Euler-integrator (modified with the discretisation factor)
  // to do the testing:
  DL_m_integrator *save_int=DL_dsystem->get_integrator();
  DL_euler my_int;
  DL_dsystem->set_integrator(&my_int);

  // first force reintegration with the new integrator
  calc_all_errors(&org_err);
  // then establish the baseline
  begin_test();
  calc_all_errors(&org_err);
  end_test();
  
  // then traverse all restrictions of all constraints and apply
  // a testvector and observe the difference in error it yields

  dCdR->makezero();
  DL_constraint *constr=(DL_constraint *)c->getfirst();
  while (constr) {
    test_restr.resize(constr->dim);
    test_restr.makezero();
    for (i=0;i<constr->dim;i++) {
       begin_test();
       test_restr.set(i,1.0); // test with unit vector i
       constr->apply_restriction_changes(&test_restr);
       calc_all_errors(&new_err);
       new_err.minus(&org_err,&err_dif);
       dCdR->setcolumn(constr->index+i,&err_dif);
       test_restr.set(i,0.0); // make zero vector again
       end_test();            // restore the state
    }
    constr=(DL_constraint*)c->getnext(constr);
  }

  // restore the motion integrator:
  DL_dsystem->set_integrator(save_int);
  c_changed=FALSE;
                                     #ifdef DCDR
                                       DL_dsystem->get_companion()->Msg("dCdR (empirical)\n");
				       dCdR->show();
                                     #endif
}

void DL_constraint_manager::calc_dCdR_full(){
  DL_constraint_pair *cpe;
  
  // first clear the old cp:
  cp.delete_all();
  
  // then calculate dCdR analytically and build cp:
  DL_constraint* cc=(DL_constraint *)c->getfirst();
  DL_constraint* cf;
  DL_largematrix sub;
  int i=0;
  dCdR->resize(totdim,totdim);
  while (cc) {
    cf=(DL_constraint *)c->getfirst();
    while (cf) {
      sub.resize(cc->dim,cf->dim);
      if (cf->dCdRsub(cc,&sub)) {
        dCdR->setsubmatrixnonzero(cc->index,cf->index,&sub);
	cpe=new DL_constraint_pair(cc,cf);
	i++;
	cp.addelem(cpe);
      }
      else dCdR->setsubmatrixzero(cc->index,cf->index,cc->dim,cf->dim);
      cf=(DL_constraint*)c->getnext(cf);
    }
    cc=(DL_constraint*)c->getnext(cc);
  }
//  if (dCdR->get_solve_method()==lud_bcksub) // always sort: sm might change
  sort_constraints();
  dCdR->analyse_structure();
                                     #ifdef DCDR
                                       DL_dsystem->get_companion()->Msg("-=-\ndCdR (analytical):\n");
                                       dCdR->show();
                                     #endif
  // if we had to calculate dCdR empirically: do so:
  if (!analytical) calc_dCdR_empirical();
  c_changed=FALSE;
}

void DL_constraint_manager::calc_dCdR_analytical(){
// using the cp list calculated by calc_dCdR_full, calculate
// dCdR analytically
// let the constraints do all of the work...
  static DL_largematrix sub;
  DL_constraint *cc,*cf;
  DL_constraint_pair *cpe=(DL_constraint_pair*)cp.getfirst();
  dCdR->makezero();
  int i=0;
  while (cpe) {
    cc=cpe->cc; cf=cpe->cf;
    sub.resize(cc->dim,cf->dim);
    cf->dCdRsub(cc,&sub);
    dCdR->setsubmatrix(cc->index,cf->index,&sub);
    cpe=(DL_constraint_pair*)cp.getnext(cpe);
    i++;
  }
                                     #ifdef DCDR
                                       DL_dsystem->get_companion()->Msg("dCdR (analytical):\n");
				       dCdR->show();
                                     #endif
  c_changed=FALSE;
}

#include "queue.h"

void DL_constraint_manager::sort_constraints() {
  // based on the connectivity info in cp, sort the constraints so that
  // the dCdR matrix has minimal bandwidth (or an approximation of that).
  // also re-orders dCdR accordingly. This is useful for sparse LU
  // decomposition.
  // The constraints are sorted into their CutHillMcKee ordering
  // (starting with the node that is the first leaf encountered in the
  //  depth-first traversal started from the last node)

  // nr of constraints:
  int N=c->length();

  if (N<2) return;
  
  int i; DL_constraint *cc; DL_constraint_pair *cpe;
  
  // old indices of the constraints:
  int *oldindex=new int[N];
    
  // graph representation:
  int *nrneighbours=new int[N];
  for (i=0;i<N;i++) nrneighbours[i]=0;
  int **neighbour=new int*[N];
  for (i=0;i<N;i++) neighbour[i]=new int[N];

  // save the old indices and number the constraints sequentially:
  i=0; cc=(DL_constraint*)c->getfirst();
  while (cc) {
    oldindex[i]=cc->index;
    cc->index=i;
    cc=(DL_constraint*)c->getnext(cc); i++;
  }
    
  //first build up graph representation:
  cpe=(DL_constraint_pair *)cp.getfirst();
  while (cpe) {
    if (cpe->cc->index!=cpe->cf->index) {
      // add cpe->cf to the neighbours of cpe->cc
      neighbour[cpe->cc->index][nrneighbours[cpe->cc->index]]=cpe->cf->index;
      nrneighbours[cpe->cc->index]++;
    }
    cpe=(DL_constraint_pair *)cp.getnext(cpe);
  }
  
  // to use CutHillMcKee ordering, sort the neighbours of each
  // node by ascending degree here:
  for (i=0;i<N;i++) {
    int j,k,*neighbouri=neighbour[i];  // simply use bubblesort
    for (j=1;j<nrneighbours[i];j++) {
      k=j;
      while (k>0) {
	if (nrneighbours[neighbouri[k-1]]>nrneighbours[neighbouri[k]]){
	  int tmp=neighbouri[k-1];
	  neighbouri[k-1]=neighbouri[k];
	  neighbouri[k]=tmp;
	  k=0;
	}
	else k--;
      }
    }
  }
  
  int *map=new int[N];
  for (i=0;i<N;i++) map[i]=-1; // not yet mapped
  
  // do the breadth first search:
  {
    DL_Queue<int> q;
    int node,j=0,mn;
    do {
      // find an unmapped node with the minimum number of neighbours
      // to start searching from:
      node=-1; mn=N;
      for (i=0;i<N;i++) {
		  if (map[i]<0) { // node is unmapped
			  if (nrneighbours[i]<mn) { // less neighbours than current minimum
				  node=i;
				  mn=nrneighbours[i];
			  }
		  }
      }
      if (node>=0) {
	// do a breadth first search for this part of the graph,
	// starting with the node we just found:
	q.put(node);
	while (q.size()) {
	  node=q.get();
	  map[node]=j++;
	  for (i=0;i<nrneighbours[node];i++) {
	    int nb=neighbour[node][i];
	    if ((map[nb]<0) && !q.contains(nb)) q.put(nb);
	  }
	}
      } else break;
    } while (TRUE);
  }

  for (i=0;i<N;i++) delete[] neighbour[i];
  delete[] neighbour;
  delete nrneighbours;

  // now we have the map: reorder the constraints and dCdR:
  int *invmap=new int[N];
  for (i=0;i<N;i++) invmap[map[i]]=i;
  delete[] map;

  // old constraints:
  DL_constraint* *oldc=new DL_constraint*[N];
  cc=(DL_constraint*)c->getfirst(); i=0;
  while (cc) {
    oldc[i]=cc;
    cc=(DL_constraint*)c->getnext(cc); i++;
  }
  
  // new constraint order:
  DL_List *newc=new DL_List;

  for (i=0;i<N;i++) newc->addelem(oldc[invmap[i]]);
  delete[] invmap;
  
  // renumber the indices of the constraints:
  int *newindex=new int[N];
  totdim=0;
  cc=(DL_constraint *)newc->getfirst();
  while (cc) {
    newindex[cc->index]=totdim;
    totdim+=cc->dim;
    cc=(DL_constraint*)newc->getnext(cc);
  }
  
  // re-arrange dCdR:
  DL_largematrix* newdCdR=new DL_largematrix(totdim,totdim,dCdR->get_solve_method());
  newdCdR->makezero();
  {
    DL_largematrix sub;
    cpe=(DL_constraint_pair *)cp.getfirst();
    while (cpe) {
      sub.resize(cpe->cc->dim,cpe->cf->dim);
      dCdR->getsubmatrix(oldindex[cpe->cc->index],
			 oldindex[cpe->cf->index],
			 &sub);
      newdCdR->setsubmatrixnonzero(newindex[cpe->cc->index],
				   newindex[cpe->cf->index],
				   &sub);
      cpe=(DL_constraint_pair *)cp.getnext(cpe);
    }
    if (newdCdR->get_min_solve_method()!=dCdR->get_min_solve_method()) {
      newdCdR->set_min_solve_method(dCdR->get_min_solve_method());
    }
  }
//  newdCdR->show_all();

  // ok: everything re-arranged: now finalize everything:
  if (newdCdR->get_bandwidth()<dCdR->get_bandwidth()) {
    delete dCdR; dCdR=newdCdR;
    delete c; c=newc;
  
    cc=(DL_constraint *)c->getfirst();
    while (cc) {
      cc->index=newindex[cc->index];
      cc=(DL_constraint*)newc->getnext(cc);
    }
  }
  else {
    delete newc;
    delete newdCdR;
    // rebuild old c:
    delete c; c=new DL_List;
    for (i=0;i<N;i++){
      oldc[i]->index=oldindex[i];
      c->addelem(oldc[i]);
    }
  }
  
  delete[] oldc;
  delete[] oldindex;
  delete[] newindex;
}

void DL_constraint_manager::add_collision(DL_collision *col){
  if (size_collisions==nrcollisions) {
    // have to increase the size of collisions:
    DL_collision* *newcollisions=new DL_collision*[size_collisions+10];
    for(int i=0;i<size_collisions;i++) newcollisions[i]=collisions[i];
    if (size_collisions) delete[] collisions;
    size_collisions+=10;
    collisions=newcollisions;
  }
  collisions[nrcollisions]=col;
  nrcollisions++;
}

void DL_constraint_manager::show_constraint_forces() {
  if (show_con_forces) return;
  DL_constraint *cc=(DL_constraint*)c->getfirst();
  while (cc) {
    cc->show_forces();
    cc=(DL_constraint*)c->getnext(cc);
  }
  show_con_forces=TRUE;
}

void DL_constraint_manager::hide_constraint_forces() {
  if (!show_con_forces) return;
  DL_constraint *cc=(DL_constraint*)c->getfirst();
  while (cc) {
    cc->hide_forces();
    cc=(DL_constraint*)c->getnext(cc);
  }
  show_con_forces=FALSE;
}

void DL_constraint_manager::solve_using_lud(){
    dCdR->set_min_solve_method(lud_bcksub);
}

void DL_constraint_manager::solve_using_cg(){
  dCdR->set_min_solve_method(conjug_grad);
  nr_cg=0;
}

void DL_constraint_manager::solve_using_svd(){
  dCdR->set_min_solve_method(svd);
  nr_cg=0;
}

#undef DCDR
//#undef DEBUG
