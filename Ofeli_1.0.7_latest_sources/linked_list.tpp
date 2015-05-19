/****************************************************************************
**
** Copyright (C) 2010-2012 Fabien Bessy.
** Contact: fabien.bessy@gmail.com
**
** This file is part of project Ofeli.
**
** http://www.cecill.info/licences/Licence_CeCILL_V2-en.html
** You may use this file under the terms of the CeCILL license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Fabien Bessy and its Subsidiary(-ies) nor the
**     names of its contributors may be used to endorse or promote products
**     derived from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
**
****************************************************************************/

namespace ofeli
{

///////////////////// list initialization/destruction /////////////////////

template <typename T>
list<T>::list() : head(new Node( T(), NULL ,NULL)), tail(head), listSize(0)
{
}

template <typename T>
list<T>::list(const list& copied)
{
    push_front(copied);
}

template <typename T>
list<T>& list<T>::operator=(const list& rhs)
{
    // 2 implementations

    // assignment with copy-and-swap idiom
    /*
    list<T> temp(rhs);
    list<T>::swap(*this,temp);
    return *this;
    */

    // optimized assignment which minimizes the number of new and delete
    if( this != &rhs ) // no auto-affectation
    {
        iterator       it_lhs = this->begin();
        const_iterator it_rhs = rhs.begin();

        // while a node exists in each list
        while( !it_lhs.end() && !it_rhs.end() )
        {
            // copies data without creating a new node
            *it_lhs = *it_rhs;

            ++it_lhs;
            ++it_rhs;
        }

        if( it_rhs.end() )
        {
            if( !it_lhs.end() ) // it means this->size() > rhs.size()
            {
                // function 'erase_after' is faster than 'erase'
                // specially for large object of type T
                // because it needn't copying data field of the next node
                while( !it_lhs.node->next->end() )
                {
                    // ==> it deletes nodes
                    it_lhs = erase_after(it_lhs);
                }
                erase(it_lhs); // deletes the last node too
            }
        }
        else
        {
            if( it_lhs.end() )
            {
                // it means this->size() < rhs.size()
                // ==> it creates new nodes with the value *it_rhs

                // it_lhs is at the end of list *this,
                // it_lhs is the sentinel node position
                it_lhs = insert_before(it_lhs,*it_rhs);
                ++it_rhs;

                for(  ; !it_rhs.end(); ++it_rhs )
                {
                    // now, it_lhs is the position of the last "true" node
                    // before the sentinel node
                    it_lhs = insert_after(it_lhs,*it_rhs);
                }
            }
        }
    }

    return *this;
}

template <typename T>
list<T>::~list()
{
    // while the list is not empty (including the sentinel node)
    while( head != NULL )
    {
        pop_front();
    }
    head = NULL;
    tail = NULL;
}

template <typename T>
void list<T>::clear()
{
    // while the list is not empty (excluding the sentinel node)
    while( !empty() )
    {
        pop_front();
    }

    return;
}


///////////////////////// list moving and reading /////////////////////////

template <typename T>
typename list<T>::iterator list<T>::begin()
{
    return head;
}
template <typename T>
typename list<T>::iterator list<T>::end()
{
    return tail;
}

template <typename T>
typename list<T>::const_iterator list<T>::begin() const
{
    return head;
}

template <typename T>
typename list<T>::const_iterator list<T>::end() const
{
    return tail;
}

////////////////////////////// list changes ///////////////////////////////

template <typename T>
inline void list<T>::assign(const T array[], int array_length)
{
    clear();

    iterator it_this = begin();
    int index = 0;

    if( index < array_length )
    {
        // list *this is empty so
        // it_this is at the end of list *this,
        // it_this is the sentinel node position
        it_this = insert_before(it_this,array[index]);
        index++;
    }

    for(  ; index < array_length; index++ )
    {
        // now, it_this is the position of the last "true" node
        // before the sentinel node
        it_this = insert_after(it_this,array[index]);
    }

    return;
}

template <typename T>
inline void list<T>::assign(int n, const T& value)
{
    clear();

    for( int iteration = 0; iteration < n; iteration++ )
    {
        push_front(value);
    }

    return;
}

template <typename T>
inline void list<T>::pop_front()
{
    bool last = false;
    if(head == tail)last=true;

    Link new_head = head->next;

    delete head;

    head = new_head;
    if(last) tail = head;

    listSize--;

    return;
}

template <typename T>
inline void list<T>::push_front(const T& value)
{
    Node* newNode = new Node(value,head,NULL);

    head->previous = newNode;

    Link aux = head;

    head = newNode;

    if(aux->next == NULL) tail = head;

    listSize++;

    return;
}

template <typename T>
inline void list<T>::push_front(const list& copied)
{
    copied.tail->next = head;
    head =  copied.head;
    listSize += copied.size();

    return;
}

template <typename T>
template <typename BinaryPredicate>
inline void list<T>::put_away(const T& value, BinaryPredicate compare)
{
    iterator position = begin();

    while( !position.end() )
    {
        if( compare(value,*position) )
        {
            break;
        }

        ++position;
    }

    insert_before(position,value);

    return;
}

template <typename T>
inline void list<T>::put_away(const T& value)
{
    put_away( value, less<T>() );
    return;
}

template <typename T>
inline typename list<T>::iterator list<T>::insert_before(iterator position, const T& value)
{
    Link previousNode = position.node->previous;

    Node* newNode = new Node(value,position.node,previousNode);
    previousNode->next = newNode;

    position.node->previous = newNode;

    listSize++;

    return ++position;
}

template <typename T>
inline typename list<T>::iterator list<T>::insert_after(iterator position, const T& value)
{
    Node* newNode = new Node(value,position.node->next, position.node);

    Link nextNode = position.node->next;
    nextNode->previous = newNode;
    position.node->next = newNode;

    if(tail == position.node) tail = newNode;

    listSize++;

    return newNode->next;
}

template <typename T>
inline typename list<T>::iterator list<T>::erase(iterator position)
{    
    Link nextNode = position.node->next;
    Link previousNode = NULL;

    if(position.node->previous != NULL )
    {
        previousNode = position.node->previous;
        previousNode->next = position.node->next;
    }
    else
    {
        head = nextNode;
    }
    nextNode->previous = previousNode;

    Link deleteNode = position.node;
    if(deleteNode == tail) tail = previousNode;

    listSize--;

    return nextNode;
}

template <typename T>
inline typename list<T>::iterator list<T>::erase_after(iterator position)
{
    Link nextNextNode = position.node->next->next;
    nextNextNode->previous = position.node;

    Link nextNode = position.node->next;
    position.node->next = nextNextNode;

    if(nextNode == tail) tail = position.node;

    listSize--;

    return nextNextNode;
}

template <typename T>
inline typename list<T>::iterator list<T>::set_dump_after(iterator position)
{
    Link dump_node = new Node( T(), NULL,position.node );

    //Save the iterator next position
    Link current_pos = position.node;
    ++position;

    current_pos->next = dump_node;

    return position;
}

template <typename T>
void list<T>::collectList(std::vector<list<int>* >* splitedList,Node*** &sublistHead, int** &sublistHeadPosition,int list, int numThreads)
{
    head = splitedList->at(0)->head;
    sublistHead[list][0]=head;
    listSize=0;
		
    for(int i=0; i < splitedList->size(); i++)
    {
        sublistHead[list][i]= splitedList->at(i)->head;
        sublistHeadPosition[list][i] = listSize +1;

        listSize+=  splitedList->at(i)->size();

        if(i < splitedList->size()-1)
		{

			splitedList->at(i)->tail->next = splitedList->at(i+1)->head;            
        	splitedList->at(i+1)->head->previous = splitedList->at(i)->tail;
        }
        else
        {
            tail = splitedList->at(i)->tail;
        }
    }

    // Recalculate sublists head pointers positions
    int nElements = (listSize/numThreads);
	if(nElements > numThreads)
	{
		int r=0;
		if((listSize % numThreads) >= 1) r=1;
		 nElements += r;
	}	
	
	for(int i=1; i < numThreads; i++)
    {
       
		int posHead = sublistHeadPosition[list][i];
        if(posHead < (i * nElements) +1)
        {
            for(; posHead <= (i * nElements); posHead++)
            {   
				sublistHead[list][i] = sublistHead[list][i]->next;
                sublistHeadPosition[list][i]++;
            }
        }
        else
        {
            for(; posHead > ((i * nElements)+1); posHead--)
            {
                sublistHead[list][i] = sublistHead[list][i]->previous;
                sublistHeadPosition[list][i]--;
            }
        }
    }

}


template <typename T>
inline void list<T>::splitList(std::vector< list<T> *>& splited_List, int numThreads, Node*** sublistHead, int** sublistHeadPosition,int list)
{
    int listNumber=0, cont=0,i=1;

    splited_List[listNumber]->head = sublistHead[list][0];
    listNumber++;

    for(; i< numThreads; i++)
    {
        Link position = sublistHead[list][i];
        position = position->previous;

        // Set the tail and the size of the previous sublist
        splited_List[listNumber-1]->tail = position;
        splited_List[listNumber-1]->listSize = sublistHeadPosition[list][i] -1 - cont;
        cont = cont + splited_List[listNumber-1]->listSize;

        position = set_dump_after(position).node;

        // Set the head
        splited_List[listNumber]->head = position;

        // Invalide the previous pointer
        position->previous = NULL;

        listNumber++;
    }

    // Set the tail and the size of the last sublist
    splited_List[listNumber-1]->tail = end().node;
    splited_List[listNumber-1]->listSize = size() -cont;
}

template <typename T>
void list<T>::initialize_sublist_control(int listNumber, int numThreads, Node*** sublistHead, int** sublistHeadPosition)
{
    //Calculation of the number of elements assigned to each thread
    int nElementsPerThread = (listSize/numThreads)+(listSize%2);

    int cont=0,j=0;
    Link position = head;

    sublistHead[listNumber][j] =  position;
    sublistHeadPosition[listNumber][j] = cont+1;
    position = position->next;
    j++;
    cont++;

    while( position->next != NULL )
    {
        if((cont % nElementsPerThread) ==0)
        {
            sublistHead[listNumber][j] =  position;
            sublistHeadPosition[listNumber][j] = cont+1;
            j++;
        }
        cont++;
        position = position->next;
    }
}

//////////////////////////// list information /////////////////////////////

template <typename T>
inline bool list<T>::empty() const
{
    if( head->end() ) // if the head node is the sentinel
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <typename T>
int list<T>::size() const
{
    return listSize;
}

template <typename U>
bool operator==(const list<U>& lhs, const list<U>& rhs)
{
    typename list<U>::const_iterator it_lhs = lhs.begin();
    typename list<U>::const_iterator it_rhs = rhs.begin();

    while( !it_lhs.end() && !it_rhs.end() )
    {
        if( *it_lhs != *it_rhs )
        {
            return false;
        }

        ++it_lhs;
        ++it_rhs;
    }

    return it_lhs.end() && it_rhs.end();
}

template <typename U>
bool operator!=(const list<U>& lhs, const list<U>& rhs)
{
    return !( lhs == rhs );
}

template <typename U>
std::ostream& operator<<(std::ostream& os, const list<U>& displayed)
{
    os << std::endl << " -------------" << std::endl;
    os << " index | value" << std::endl;
    os << " -------------" << std::endl;

    int index = 0;
    for( typename list<U>::const_iterator position = displayed.begin(); !position.end(); ++position )
    {
        os.width(6);
        os << std::right << index++ << " | " << *position << std::endl;
    }

    os << " -------------" << std::endl;

    return os;
}

template <typename T>
void list<T>::display() const
{
    std::cout << *this;
}

template <typename T>
T* list<T>::get_array(int array_length) const
{
    T* array = new T[array_length];

    int index = 0;
    for( const_iterator position = begin(); !position.end(); ++position )
    {
        array[index++] = *position;
    }

    return array;
}

template <typename T>
T list<T>::getFirst() const
{
   return **head;
}

template <typename T>
T list<T>::getLast() const
{
   return **tail;
}

}
