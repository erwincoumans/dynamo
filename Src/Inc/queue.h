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
// filename	: dqueue.h
// description	: support for sparse largematrix routines
// author	: Mischa Courtin '96
//

#ifndef DL_Queue_H
#define DL_Queue_H

#include <stdlib.h>

template <class T> class DL_Queue {
   public:
      T get();
      void put(T);
      int size(void);
      int contains(T);
      DL_Queue();
      ~DL_Queue();
   private:
      typedef struct Elem {T x; struct Elem* next;} Elem;
      Elem* head;
      Elem* tail;
      int   nrelem;
};


template <class T> DL_Queue<T>::DL_Queue(){
   head = tail = NULL;
   nrelem = 0;
}


template <class T> DL_Queue<T>::~DL_Queue(){
   while (head) {
      Elem* tmp = head;
      head      = head->next;
      delete tmp;
   }
}


template <class T> inline int DL_Queue<T>::size(void){
   return nrelem;
}


template <class T> inline void DL_Queue<T>::put(T x)
{
   Elem* tmp = new Elem;
   tmp->x    = x;
   tmp->next = NULL;

   if (nrelem) tail= (  tail->next = tmp ) ;
   else tail=head=tmp;

   nrelem++;
}


template <class T> inline T DL_Queue<T>::get()
{
   if (nrelem) {
      Elem* tmp = head;
      head      = head->next;
      nrelem--;

      if (nrelem==0) tail = NULL;

      T ret = tmp->x;
      delete tmp;
      return ret;
   }
   else return NULL;
}


template <class T> int DL_Queue<T>::contains(T t)
{
   for(Elem* p=head;p;p=p->next) if (p->x == t) return 1;
   return 0;
}

#endif
