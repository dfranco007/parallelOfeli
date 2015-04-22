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

#ifndef LINKED_LIST_HPP
#define LINKED_LIST_HPP

#include <vector>
#include <cstdlib> // NULL macro
#include <iostream> // for operator<< overloading

namespace ofeli
{

template <typename T = int>
class list
{

    //////////////////////   Nested structure and class   /////////////////////



        //////////////////////////   Declarations   ///////////////////////////

public :

    class const_iterator;

    class iterator;

    //! Element of the list.
    struct Node;

    //! Link is a pointer to a Node.
    typedef Node* Link;

    //! This structure implements a node of class \a list. It is composed by an element of type \a T and a pointer on the next node.
    struct Node
    {
        //! Constructor.
        Node(const T& data1, Link next1, Link previous1) : data(data1), next(next1), previous(previous1)
        {
        }
        //! Constructor.
        Node(Node &node1) : data(node1.data), next(node1.next), previous(node1.previous)
        {
        }

        //! Gets the node data of the Node
        T operator*() const
        {
            return data;
        }
        //! Checks if the node is at the end of the list, i.e. if the node is the sentinel node.
        bool end() const
        {
            return ( next == NULL) ? true : false;
        }

        //! Element storage.
        T data;
        //! Pointer to the next node.
        Link next;
        //! Pointer to the previous node.
        Link previous;
    };

        ///////////////////////////////////////////////////////////////////////


        //////////////////////////   Definitions   ///////////////////////////

public :

    //! const iterator to read a const list.
    class const_iterator
    {

        //! list has access to the private members of \a const_iterator.
        friend class list;

    public :

        //! Gets the node data. \a const_iterator protects the data against writing.
        const T& operator*() const
        {
            return node->data;
        }

        //! Checks if the \a const_iterator is at the end of the list, i.e. if the node is the sentinel node.
        bool end() const
        {
            return ( node->end() ) ? true : false;
        }

        //! Increment \a const_iterator.
        const_iterator& operator++()
        {
            node = node->next;
            return *this;
        }

    private :

        //! Constructor with a node.
        const_iterator(Link node1) : node(node1)
        {
        }

        //! \a const_iterator encapsulates a pointer to a node.
        Link node;
    };

    //! Iterator to modify a list.
    class iterator
    {

        //! \a list has access to the private members of iterator.
        friend class list;

    public :

        //! Gets the node data. \a iterator does not protect the data against writing.
        T& operator*() const
        {
            return node->data;
        }

        //! Checks if the iterator is at the end of the list, i.e. if the node is the sentinel node.
        bool end() const
        {
            return ( node->end() ) ? true : false;
        }

        //! Increment the iterator.
        iterator& operator++()
        {
            node = node->next;
            return *this;
        }

    private :

        //! Constructor with a node.
        iterator(Link node1) : node(node1)
        {
        }
        //! \a const_iterator encapsulates a pointer to a node.
        Link node;
    };



        ///////////////////////////////////////////////////////////////////////




    ///////////////////////////////////////////////////////////////////////////

public :

    ///////////////////// list initialization/destruction /////////////////////

    //! Constructor.
    list();

    //! Copy constructor.
    list(const list& copied);

    //! Assignment operator overloading.
    list& operator=(const list& rhs);

    //! Destructor. All the elements in the list container are dropped (including the sentinel node) : their destructors are called, and then they are removed from the list container, leaving it with a size of 0.
    ~list();

    //! All the elements in the list container are dropped (except the sentinel node) : their destructors are called, and then they are removed from the list container, leaving it with a size of 0.
    void clear();

    ///////////////////////// list moving and reading /////////////////////////

    //! Returns the head link.
    iterator begin();

    //! Returns the head link.
    iterator end();

    const_iterator end() const;

    //! Returns the head link.
    const_iterator begin() const; // const overloading


    ////////////////////////////// list changes ///////////////////////////////

    //! Assigns new content to the container, dropping all the elements contained in the container object before the call and replacing them by those specified by the parameters:
    //! - the new content is the repetition \a n times of copies of element value.
    void assign(int n, const T& value);

    //! Assigns new content to the container, dropping all the elements contained in the container object before the call and replacing them by those specified by the parameters:
    //! - the new content is a copy of an array.
    void assign(const T array[], int array_length);

    //! Removes the first element in the \a list container, effectively reducing its size by one.
    void pop_front();

    //! Inserts a new element at the beginning of the list.
    void push_front(const T& value);

    //! Inserts copies of list \a copied elements at the beginning of the list \a *this.
    void push_front(const list& copied);

    //! Puts away a new element to maintain sorted list into the specified order.
    template <typename BinaryPredicate>
    void put_away(const T& value, BinaryPredicate compare);

    //! Puts away a new element to maintain sorted list into ascending order.
    void put_away(const T& value);

    //! Inserts a new element before \a position. Returns the position of the inserted element.
    iterator insert_before(iterator position, const T& value);

    //! Inserts a new element after \a position. Returns the position of the inserted element.
    iterator insert_after(iterator position, const T& value);

    //! Removes from the list container the element at \a position and returns a valid iterator, i.e. the position of the next element.
    iterator erase(iterator position);

    //! Removes from the list container the element at \a ++position and returns the position.
    iterator erase_after(iterator position);

    //! Sets a dump element after position
    iterator set_dump_after(iterator position);

    //! Collect all subList of splitedList into the list which calls the function
    void collectList(std::vector<list<int>* >* splitedList,Node*** &sublistHead, int** &sublistHeadPosition,int list, int numThreads);

    //! Splits de list in a number of parts equals number of threads
    //! without modifying the original list
    inline void splitList(std::vector< list<T> *>& splited_List, int numThreads, Node*** sublistHead, int** sublistHeadPosition,int list);

    //! Initializates necessary control of sublists to construct them
    void initialize_sublist_control(int listNumber, int numThreads, Node*** sublistHead, int** sublistHeadPosition);


    //////////////////////////// list information /////////////////////////////

    //! Returns whether the list container is empty, i.e. whether its size is 0.
    bool empty() const;

    //! Returns the number of elements in the list without counting the sentinel node.
    //! The complexity of this function is in \a O(n) so if you want to get often the size of the list, you should update the size in a variable after each modification.
    int size() const;

    //! \a Equal \a to operator overloading.
    template<class U>
    friend bool operator==(const list<U>& lhs, const list<U>& rhs);

    //! \a Not \a equal \a to operator overloading.
    template<class U>
    friend bool operator!=(const list<U>& lhs, const list<U>& rhs);

    //! Overloading of cout <<. It displays a linked list in the same way as integral-type variable.
    template<class U>
    friend std::ostream& operator<<(std::ostream& os, const list<U>& displayed);

    //! A second way to display a linked list.
    void display() const;

    //! Gets an array with the value of list \a *this.
    T* get_array(int array_length) const;

    //! Gets the first element
    T getFirst() const;

    //! Gets the last element
    T getLast() const;

    ///////////////////////////////////////////////////

private :

    //! Head of the list.
    Link head;

    //!Tail of the list
    Link tail;

    //! List size
    int listSize;
};

// function object definitions for the template function parameter BinaryPredicate
template <typename T = int>
struct equal_to { bool operator()(const T& x, const T& y) const { return x == y; } };

template <typename T = int>
struct not_equal_to { bool operator()(const T& x, const T& y) const { return x != y; } };

template <typename T = int>
struct greater { bool operator()(const T& x, const T& y) const { return x > y; } };

template <typename T = int>
struct less { bool operator()(const T& x, const T& y) const { return x < y; } };

template <typename T = int>
struct greater_equal { bool operator()(const T& x, const T& y) const { return x >= y; } };

template <typename T = int>
struct less_equal { bool operator()(const T& x, const T& y) const { return x <= y; } };

// function object definition for the template function parameter UnaryPredicate
template <typename T = int>
struct predicate_example
{
    bool operator()(const T& value) const
    {
        if( value > T(10) && value < T(20) && value != T(13) )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

}

// list definitions
#include "linked_list.tpp"

#endif
