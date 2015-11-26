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

/**************************************************************************
Filename    : list.h
Last update : jan 1997 
Author      : B Barenbrug, based on the list-class by H.M. Roelfs 
**************************************************************************/

#ifndef DL_LISTH
#define DL_LISTH

#include <string.h>
#include <iostream.h>
#include <stdlib.h>

//*************************************************************************

class DL_ListElem {
  public:
    DL_ListElem*  next;          // pointer to next element of list
    DL_ListElem*  prev;          // pointer to previous element of list

    DL_ListElem(){next=prev=NULL;}; // Constructor
    virtual ~DL_ListElem(){};       // Destructor Must be vitual for delete_all
                                    // to work!
};


//*********************************************************************

class DL_List {
  private :
	DL_ListElem* head;   // pointer to head of list
	DL_ListElem* tail;   // pointer to tail of list
	int       number;    // number elements of list

  public :
	void  addelem(DL_ListElem*);       // Add an element to the list
        void  addbefore(DL_ListElem*,DL_ListElem*); // add an element before
                                           // another element in the list
        void  addafter(DL_ListElem*,DL_ListElem*); // add an element after
                                           // another element in the list
	
	void  remelem(DL_ListElem*);       // Remove list element (not deleted)
	void  rem_all();                   // remove all list elements
	void  delete_all();                // delete all list elements
	int length();                      // return the length of the list
	
	DL_ListElem* getfirst();           // Get first list element
	DL_ListElem* getnext(DL_ListElem*);// Get next list element

	DL_ListElem* getlast();            // Get last list element
	DL_ListElem* getprev(DL_ListElem*);// Get previous list element
	
	DL_ListElem* element(int);         // Get a list element
	
        DL_List();              // Constructor
        ~DL_List(){};           // Destructor
}; 

inline DL_List::DL_List() { 
  head=NULL; tail=NULL;
  number=0;
};

inline int DL_List::length() {
  return number;
};

inline DL_ListElem* DL_List::getfirst() {
  return head;
};

inline DL_ListElem* DL_List::getnext(DL_ListElem* le) {
  if (le) return le->next;
  else return head;
};

inline DL_ListElem* DL_List::getlast() {
  return tail;
};

inline DL_ListElem* DL_List::getprev(DL_ListElem* le) {
  if (le) return le->prev;
  else return tail;
};

inline void DL_List::addelem(DL_ListElem* le) {
  DL_ListElem* prev=tail;
  tail=le;

  if (prev) {
	tail->prev=prev;
	tail->next=NULL;
	prev->next=tail;
  }
  else {
    tail->prev=NULL; // added by E. Peeters, 11 oct 93
    tail->next=NULL; // added by E. Peeters, 11 oct 93
	head=tail;
  };

  number++;
}

inline void DL_List::addafter(DL_ListElem *le, DL_ListElem *el) {
// add le after element el
  if (el) {
    if (el==tail) addelem(le); // add at the tail
    else {
      le->next=el->next;
      el->next=le;
      le->next->prev=le;
      le->prev=el;
      number++;
    }
  }
  else { // adding after NULL means adding add the head
    le->prev=NULL;
    le->next=head;
    if (head) {
      head->prev=le;
      head=le;
    }
    else head=tail=le;
    number++;
  }
}

inline void DL_List::addbefore(DL_ListElem *le, DL_ListElem *el) {
// add le before element el
  if (el) {
    if (el==head) addafter(le,NULL); // add at the head
    else {
      le->prev=el->prev;
      el->prev=le;
      le->prev->next=le;
      le->next=el;
      number++;
    }
  }
  else addelem(le); // adding before NULL means adding add the end
}

inline void DL_List::remelem(DL_ListElem* le) {
  // changed by E. Peeters, 11 oct 93 
  if ( (head == le) && (tail == le)) {
     head=NULL;
     tail=NULL;
  }
  else {
    if (head==le) {
       head=head->next;
       if (head) head->prev=NULL;
    }
    else {
      if (tail==le) {
	 if (tail->prev) tail->prev->next=NULL;
	 tail=tail->prev;
      }
      else {
	 le->prev->next=le->next;
	 if (le->next) le->next->prev=le->prev;
      };
    };
  };
  number--;
  // maybe delete le?
}

inline DL_ListElem* DL_List::element(int i) {
  if (i < number) {
    int count=0;	 
	DL_ListElem* le=head;
	while (count < i) {
          le=le->next;
	  count++;
	};
    return le;
  }
  else return NULL;
}

inline void DL_List::rem_all() {
  while (head) remelem(head);
  number=0;
}

inline void DL_List::delete_all() {
  DL_ListElem *tmp2,*tmp=head;
  while (tmp) {
     tmp2=tmp->next;
     delete tmp;
     tmp=tmp2;
  }
  head=tail=NULL;
  number=0;  
}

#endif
