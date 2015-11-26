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
// filename: containerlist.h
// description: lists that are implemented through containers with
//              pointers to the actual objects. Compared to the DL_List
//              class it is a little less efficient, but it does allow
//              for objects to be part of more than one list.
//              For iteration a user will have to have its own pointer
//              to the containers: trying to hide this will mean that
//              the list will have to maintain the current position,
//              meaning that only one client can iterate at a time
//              (no support for parallellism or nested iteration)
// author: Bart Barenbrug   Jan '98
//
 
#ifndef DL_CONTAINERLISTH
#define DL_CONTAINERLISTH

template <class T>
class DL_ListElemContainer {
  public:
    DL_ListElemContainer *next,*prev; // NULL as value for either of these
                                      // means no next/prev element
    T *elem;                          // a pointer to the element

    DL_ListElemContainer(T* e){elem=e; prev=next=NULL;};   // constructor
    ~DL_ListElemContainer(){};                             // destructor
};

template <class T>
class DL_containerlist {
  private:
    DL_ListElemContainer<T> *first,*last;
  public:
  
/*| Interface:
 *|   void addelem(T*)
 *|   void addbefore(T*,DL_ListElemContainer<T>*)
 *|   void addafter(T*,DL_ListElemContainer<T>*)
 *|   void remelem(T*)
 *|   void remelem(DL_ListElemContainer<T>*)
 *|   void rem_all();
 *|   void delete_all();
 *|   int length();
 *|   T* getfirst(DL_ListElemContainer<T>**)
 *|   T* getnext(DL_ListElemContainer<T>**)
 *|   T* getlast(DL_ListElemContainer<T>**)
 *|   T* getprev(DL_ListElemContainer<T>**)
 *|   T* element(int)
 *|   assign(DL_containerlist<T>*)
 */
    void addelem(T *e) {
      // add an element of type T to end of the list
      if (last) {     // add to the end of the list
	last->next=new DL_ListElemContainer<T>(e);
	last->next->prev=last;
	last=last->next;
      }
      else { // the list is empty sofar
	first=last=new DL_ListElemContainer<T>(e);
      }
    } // addelem

    void addbefore(T *e, DL_ListElemContainer<T> *cont) {
      // add e before the element contaioned in cont
      if (cont) {
	if (cont==first) addafter(e,NULL); // add at the head
	else { // add in the middle
	  DL_ListElemContainer<T>* tmp=new DL_ListElemContainer<T>(e);
	  tmp->prev=cont->prev;
	  cont->prev=tmp->prev->next=tmp;
	  tmp->next=cont;
	}
      }
      else addelem(e); // adding before NULL means adding at the end
    } // addbefore

    void addafter(T *e, DL_ListElemContainer<T> *cont) {
      // add an element of type T after teh element contained in cont
      if (cont) {
        if (cont==last) addelem(e); // add at the tail;
	else { //add in the middle
	    DL_ListElemContainer<T>* tmp=new DL_ListElemContainer<T>(e);
	    tmp->next=cont->next;
	    cont->next=tmp->next->prev=tmp;
	    tmp->prev=cont;
	}
      }
      else { // adding after NULL means adding at the head
	DL_ListElemContainer<T>*tmp=first;
	first=new DL_ListElemContainer<T>(e);
	first->prev=NULL;
	if (tmp) {
	  first->next=tmp;
	  tmp->prev=first;
	}
	else {
	  last=first;
	  first->next=NULL;
	}
      }
    } // addafter

    void remelem(T *e) {
      // remove all occurrences of e from the list if it is in there

      // traverse the list and remove any container that contains e
      // So this one's quite expensive!

      if (!first) { // the list is already empty
	return;
      }

      DL_ListElemContainer<T> *next;

      // first see if the head of the list needs to be removed
      while (first->elem==e) { //the first element has to be deleted
	next=first->next;
	if (next) { // this isn't the only element in the list
	  next->prev=NULL;
	  delete first;
	  first=next;
	}
	else { // this was the only element of the list
	  delete first;
	  first=last=NULL;
	  return;
	}
      }
     
      // then see if the tail of the list needs to be removed
      DL_ListElemContainer<T> *prev;       
      while (last->elem==e) { //the last element has to be deleted
	prev=last->prev;
	if (prev) { // this isn't the only element in the list
	  prev->next=NULL;
	  delete last;
	  last=prev;
	  prev=last->prev;
	}
	else { // this was the final element of the list
	  delete last;
	  first=last=NULL;
	  return;
	}
      }
     
      // now we know that any elements equal to e are not the first or
      // the last, so we don't have to change first and last anymore.
      // traverse the rest of the list to find any middle elements
      // equal to e
      DL_ListElemContainer<T> *pos=first->next;

      while (pos) {
	if (pos->elem==e) {
          next=pos->next;
	  pos->prev->next=next;
	  if (next) next->prev=pos->prev;
	  delete pos;
	  pos=next;
	}
	else pos=pos->next;
      }
    } //remelem
     
    void remelem(DL_ListElemContainer<T> *cont) {
      // this version of remove already has a pointer to the container
      // which is to be removed
      if (first==cont) {
        first=first->next;
        if (first) first->prev=NULL;
        else last=NULL;
      }
      else {
        if (last==cont) {
	  last=last->prev;
	  last->next=NULL;
        }
        else {
          cont->next->prev=cont->prev;
	  cont->prev->next=cont->next;
	}
      }
      delete cont;
    } //remelem
    
    void rem_all(){
      DL_ListElemContainer<T> *next;
      while (first) {
	next=first->next;
	delete first;
	first=next;
      }
      last=NULL;
    }; // rem_all

    void delete_all(){
      DL_ListElemContainer<T> *next;
      while (first) {
	next=first->next;
	delete first->elem;
	delete first;
	first=next;
      }
      last=NULL;
    } // delete_all;

    int length(){
      DL_ListElemContainer<T> *tmp=first;
      int i=0;
      while (tmp) {
	tmp=tmp->next; i++;
      }
      return i;
    } // length
  
    T* getfirst(DL_ListElemContainer<T> **cont){
       // get the first element in the list
       // (or NULL if the list is empty)
       // return the container pointer in cont
         if (*cont=first) return first->elem;
         else return NULL;
    } // getfirst
    
    T* getnext(DL_ListElemContainer<T> **cont){
       // after having called getfirst, this method can be used
       // to traverse the list (the ListElemContainer references the
       // current position in the list). It returns NULL if there are no
       // more elements left
         if (*cont) {
           if (*cont=(*cont)->next) return (*cont)->elem;
           else return NULL;
         }
         else return NULL;
    } // getnext
    
    T* getlast(DL_ListElemContainer<T> **cont){
       // get the last element in the list
       // (or NULL if the list is empty)
       // return the container pointer in cont
         if (*cont=last) return last->elem;
         else return NULL;
    } // getlast

    T* getprev(DL_ListElemContainer<T> **cont){
       // for backwards traversing of the list
         if (*cont) {
           if (*cont=(*cont)->prev) return (*cont)->elem;
           else return NULL;
         }
         else return NULL;
    } // getprev

    T* element(int i) {
      // return the i-th element in the list (first element is element #0):
      DL_ListElemContainer<T> *tmp=first;
      
      while ((i>0) && (tmp)) {
	i--; tmp=tmp->next;
      }
      if (tmp) return tmp->elem;
      else return NULL;
    }
  
    void assign(DL_containerlist<T>* lst) {
    // assign lst to this list (note that the elements themselves are not copied
      // first, throw away everything we have:
      DL_ListElemContainer<T> *pos=first, *next;
      T *el;
      if (pos) {
        next=pos->next;
        while (next) { delete pos; pos=next; next=pos->next; }
        delete pos;
      }
      // next: add all elements of lst to myself:
      el=lst->getfirst(&pos);
      while (el) {
        addelem(el);
	el=lst->getnext(&pos);
      }
    }
    
    DL_containerlist() { first=last=NULL; };  // constructor
    
    ~DL_containerlist() {                     // destructor
       // throw away any containers that are left:
       rem_all();
    } // destructor

}; // class DL_ontainerlist

#endif
