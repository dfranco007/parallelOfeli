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

// comment/uncomment below, respectively,
// if you want to test the algorithm with a row/column wise image data buffer
//#define COLUMN_WISE_IMAGE_DATA_BUFFER
// functions affected by the define :
// - int find_offset(int x, int y) const
// - void find_xy_position(int offset, int& x, int& y) const
// for a performance reason, this functions are defined in the header file

// comment/uncomment below, respectively,
// if you want to compile the 4/8-connected neighborhood version of the algorithm
//#define VERSION_8_CONNECTED_NEIGHBORHOOD
// functions affected by the define :
// - void switch_in(int offset)
// - void switch_out(int offset)
// - bool isRedundantLinPoint(int offset) const
// - bool isRedundantLoutPoint(int offset) const

#include "activecontour.hpp"
#include <cstring> // for the function "std::memcpy"
#include <cmath> // for the function "std::exp"
#include <cstdlib> // for the function "std::abs"
#include <algorithm> // for the function "std::max"

namespace ofeli
{

ActiveContour::ActiveContour(const unsigned char* img_data1, int img_width1, int img_height1,
                             bool hasEllipse1, double init_width1, double init_height1, double center_x1, double center_y1,
                             bool hasCycle2_1, int kernel_length1, double sigma1, int Na1, int Ns1) :
    img_data(img_data1), img_width(img_width1), img_height(img_height1),
    hasCycle2(hasCycle2_1), kernel_length(kernel_length1), sigma(sigma1), Na_max(Na1), Ns_max(Ns1),
    img_size(img_width1*img_height1), kernel_radius( (kernel_length1-1)/2 ), phi( new char[img_width1*img_height1] ),
    gaussian_kernel( make_gaussian_kernel(kernel_length1,sigma1) ), iteration_max( 5*std::max(img_width1,img_height1) ),
    iteration(0), Na(0), Ns(0), previous_lists_length(99999999), lists_length(0), oscillations_in_a_row(0),
    Lout( ), Lin( )
{
    if( img_data1 == NULL )
    {
        std::cerr << std::endl <<
        "The pointer img_data1 must be a non-null pointer, it must be allocated."
        << std::endl;
    }

    initialize_phi_with_a_shape(hasEllipse1, init_width1, init_height1, center_x1, center_y1);
    initialize_lists();

    //Taking the number of threads
    #pragma omp parallel
    {
       numThreads = omp_get_num_threads();
    }

    //Create the sublists
    for(int i=0; i < numThreads; i++)
    {
        Splited_Lout.push_back(new list<int>());
        Splited_Lin.push_back(new list<int>());
    }

    //Create the sublists control pointers to change splitList to constant complexity
    sublistHead = new list<int>::Node**[2];
    for(int i=0; i < 2; i++) sublistHead[i] = new list<int>::Node*[numThreads];
    sublistHeadPosition = new int*[2];
    for(int i=0; i < 2; i++) sublistHeadPosition[i] = new int[numThreads];
    Lout.initialize_sublist_control(0,numThreads,sublistHead, sublistHeadPosition);
    Lin.initialize_sublist_control(1,numThreads,sublistHead, sublistHeadPosition);

    //Initiliaze some extra variable
    hasOutwardEvolutionInThread = new bool[numThreads];
    hasInwardEvolutionInThread = new bool[numThreads];
}

ActiveContour::ActiveContour(const unsigned char* img_data1, int img_width1, int img_height1,
                             const char* phi_init1,
                             bool hasCycle2_1, int kernel_length1, double sigma1, int Na1, int Ns1) :
    img_data(img_data1), img_width(img_width1), img_height(img_height1),
    hasCycle2(hasCycle2_1), kernel_length(kernel_length1), sigma(sigma1), Na_max(Na1), Ns_max(Ns1),
    img_size(img_width1*img_height1), kernel_radius( (kernel_length1-1)/2 ), phi( new char[img_width1*img_height1] ),
    gaussian_kernel( make_gaussian_kernel(kernel_length1,sigma1) ), iteration_max( 5*std::max(img_width1,img_height1) ),
    iteration(0), Na(0), Ns(0), previous_lists_length(99999999), lists_length(0), oscillations_in_a_row(0),
    Lout(  ), Lin( )
{
    if( img_data1 == NULL )
    {
        std::cerr << std::endl <<
        "The pointer img_data1 must be a non-null pointer, it must be allocated."
        << std::endl;
    }

    if( phi_init1 == NULL )
    {
        std::cerr << std::endl <<
        "The pointer phi_init1 must be a non-null pointer, it must be allocated."
        << std::endl;
    }

    std::memcpy(phi,phi_init1,img_size); // initialize phi with a buffer copy
    initialize_lists();

    //Taking the number of threads
    #pragma omp parallel
    {
       numThreads = omp_get_num_threads();
    }

    //Create the sublists
    for(int i=0; i < numThreads; i++)
    {
        Splited_Lout.push_back(new list<int>());
        Splited_Lin.push_back(new list<int>());
    }

    //Create the sublists control pointers to change splitList to constant complexity
    sublistHead = new list<int>::Node**[2];
    for(int i=0; i < 2; i++) sublistHead[i] = new list<int>::Node*[numThreads];
    sublistHeadPosition = new int*[2];
    for(int i=0; i < 2; i++) sublistHeadPosition[i] = new int[numThreads];
    Lout.initialize_sublist_control(0,numThreads,sublistHead, sublistHeadPosition);
    Lin.initialize_sublist_control(1,numThreads,sublistHead, sublistHeadPosition);

    //Initiliaze some extra variable
    hasOutwardEvolutionInThread = new bool[numThreads];
    hasInwardEvolutionInThread = new bool[numThreads];
}

ActiveContour::ActiveContour(const ActiveContour& ac) :
    img_data(ac.img_data), img_width(ac.img_width), img_height(ac.img_height),
    hasCycle2(ac.hasCycle2), kernel_length(ac.kernel_length), sigma(ac.sigma), Na_max(ac.Na_max), Ns_max(ac.Ns_max),
    img_size(ac.img_size), kernel_radius(ac.kernel_radius), phi( new char[ac.img_size] ),
    gaussian_kernel( make_gaussian_kernel(ac.kernel_length,ac.sigma) ), iteration_max(ac.iteration_max),
    iteration(ac.iteration), Na(ac.Na), Ns(ac.Ns),
    previous_lists_length(ac.previous_lists_length), lists_length(ac.lists_length),
    oscillations_in_a_row(ac.oscillations_in_a_row),
    isStopped(ac.isStopped), hasLastCycle2(ac.hasLastCycle2), hasListsChanges(ac.hasListsChanges),
    hasOscillation(ac.hasOscillation),
    hasOutwardEvolution(ac.hasOutwardEvolution), hasInwardEvolution(ac.hasInwardEvolution),
    Lout(ac.Lout), Lin(ac.Lin), Splited_Lout(ac.Splited_Lout), Splited_Lin(ac.Splited_Lin), numThreads(ac.numThreads),
    sublistHead(ac.sublistHead), sublistHeadPosition(ac.sublistHeadPosition),
    hasInwardEvolutionInThread(ac.hasInwardEvolutionInThread),  hasOutwardEvolutionInThread(ac.hasOutwardEvolutionInThread)
{
    if( ac.img_data == NULL )
    {
        std::cerr << std::endl <<
        "The pointer img_data of the copied active contour must be a non-null pointer, it must be allocated."
        << std::endl;
    }

    std::memcpy(phi,ac.phi,img_size); // buffer copy

}

void ActiveContour::initialize_phi_with_a_shape(bool hasEllipse, double init_width, double init_height,
                                                double center_x, double center_y)
{
    int x, y; // position of the current pixel

    // performs an ellipse
    if( hasEllipse )
    {
        for( int offset = 0; offset < img_size; offset++ )
        {
            find_xy_position(offset,x,y); // x and y passed by reference

            // ellipse inequation
            if(   square( double(y)-(1.0+2.0*center_y)*double(img_height)/2.0 )
                  / square( init_height*double(img_height)/2.0 )
                + square( double(x)-(1.0+2.0*center_x)*double(img_width)/2.0 )
                  / square( init_width*double(img_width)/2.0 )
                > 1.0
              )
            {
                phi[offset] = 1; // outside boundary value
            }
            else
            {
                phi[offset] = -1; // inside boundary value
            }
        }
    }

    // performs a rectangle
    else
    {
        for( int offset = 0; offset < img_size; offset++ )
        {
            find_xy_position(offset,x,y); // x and y passed by reference

            if(    double(y) > ((1.0-init_height)*double(img_height)/2.0 + center_y*double(img_height))
                && double(y) < ((double(img_height)-(1.0-init_height)*double(img_height)/2.0) + center_y*double(img_height))
                && double(x) > ((1.0-init_width)*double(img_width)/2.0 + center_x*double(img_width))
                && double(x) < ((double(img_width)-(1.0-init_width)*double(img_width)/2.0) + center_x*double(img_width))
              )
            {
                phi[offset] = -1; // inside boundary value
            }
            else
            {
                phi[offset] = 1; // outside boundary value
            }
        }
    }

    return;
}

void ActiveContour::initialize_lists()
{
    // each point of a list must be connected at least by one point of the other list
    // eliminate redundant points in phi if needed, and initialize Lout and Lin
    // so you can pass to the constructor phi_init1, a binarized buffer
    // with 1 for outside region and -1 for inside region to simplify your task
    for( int offset = 0; offset < img_size; offset++ )
    {
        if( phi[offset] == 1 ) // outside boundary value
        {
            if( isRedundantLoutPoint(offset) )
            {
                phi[offset] = 3; // exterior value
            }
            else
            {
                Lout.push_front(offset);
            }
        }

        if( phi[offset] == -1 ) // inside boundary value
        {
            if( isRedundantLinPoint(offset) )
            {
                phi[offset] = -3; // interior value
            }
            else
            {
                Lin.push_front(offset);
            }
        }
    }

    return;
}

void ActiveContour::initialize_for_each_frame()
{
    if( Lout.empty() && Lin.empty() )
    {
        std::cerr << std::endl <<
        "The both lists Lout and Lin are empty so the algorithm could not converge. The active contour is initialized with an ellipse."
        << std::endl;

        initialize_phi_with_a_shape(true, 0.65, 0.65, 0.0, 0.0);
        initialize_lists();
    }

    isStopped = false;
    hasLastCycle2 = false;

    // 3 stopping conditions (re)initialized
    hasListsChanges = true;

    hasOscillation = false;
    oscillations_in_a_row = 0;

    iteration = 0;

    return;
}

const int* const ActiveContour::make_gaussian_kernel(int kernel_length1, double sigma1)
{
    // kernel_length impair and strictly positive
    if( kernel_length1 % 2 == 0 )
    {
        kernel_length1--;
    }

    if( kernel_length1 < 1 )
    {
        kernel_length1 = 1;
    }

    // protection against /0
    if( sigma1 < 0.000000001 )
    {
        sigma1 = 0.000000001;
    }

    int x, y; // position of the current pixel

    const int kernel_size1 = kernel_length1*kernel_length1;
    const int kernel_radius1 = (kernel_length1-1)/2;

    int* const gaussian_kernel1 = new int[kernel_size1];

    for( int offset = 0; offset < kernel_size1; offset++ )
    {
        // offset = x+y*kernel_length so
        y = offset/kernel_length1;
        x = offset-y*kernel_length1;

        gaussian_kernel1[offset] = int(  0.5 +
                                         100000.0
                                         * std::exp( -( double( square(y-kernel_radius1)
                                                               + square(x-kernel_radius1) )
                                                     ) / (2.0*square(sigma1)) )  );
    }

    return gaussian_kernel1;
}

ActiveContour::~ActiveContour()
{
    delete[] gaussian_kernel;
    delete[] phi;
}

void ActiveContour::do_one_iteration_in_cycle1()
{/*
    std::clock_t startTime, stopTime,ya,yi;
    double elapsedTime;
    startTime = std::clock();
    ya = std::clock();
    std::cout << "/////////////////////////////////" << std::endl;
*/
    // means of the Chan-Vese model for children classes ACwithoutEdges and ACwithoutEdgesYUV
    calculate_means(); // virtual function for region-based models

    hasOutwardEvolution = false;
    hasInwardEvolution = false;

    Lout.splitList(Splited_Lout,numThreads,sublistHead, sublistHeadPosition,0);
    Lin.splitList(Splited_Lin,numThreads,sublistHead, sublistHeadPosition,1);
/*
    stopTime = std::clock();
    elapsedTime = (stopTime - startTime) / double(CLOCKS_PER_SEC);
    std::cout <<"split-"<< elapsedTime << std::endl;
*/
    int tid;
    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        hasOutwardEvolutionInThread[tid]=false;     

        //startTime = std::clock();

        for( list<int>::iterator Lout_point = Splited_Lout[tid]->begin(); !Lout_point.end(); )
        {
            if( compute_external_speed_Fd(*Lout_point, tid) > 0 )
            {
                hasOutwardEvolutionInThread[tid]=true;

                // updates of the variables to calculate the means Cout and Cin
                updates_for_means_in1(tid); // virtual function for region-based models
                Lout_point = switch_in(Lout_point,tid); // outward local movement
                // switch_in function returns a new Lout_point
                // which is the next point of the former Lout_point
            }
            else
            {
                ++Lout_point;
            }
        }
/*
        stopTime = std::clock();
        elapsedTime = (stopTime - startTime) / double(CLOCKS_PER_SEC);
        std::cout << tid << "bucle1: "<<elapsedTime << std::endl;
        startTime = std::clock();
*/
        clean_Lin(tid); // eliminate Lin redundant points
/*
        stopTime = std::clock();
        elapsedTime = (stopTime - startTime) / double(CLOCKS_PER_SEC);
        std::cout << tid << "cleanLin: " << elapsedTime << std::endl;
        startTime = std::clock();
*/
        hasInwardEvolutionInThread[tid]=false;

        for( list<int>::iterator Lin_point = Splited_Lin[tid]->begin(); !Lin_point.end(); )
        {
            if( compute_external_speed_Fd(*Lin_point, tid) < 0 )
            {
                hasInwardEvolutionInThread[tid]=true;

                // updates of the variables to calculate the means Cout and Cin
                updates_for_means_out1(tid); // virtual function for region-based models
                Lin_point = switch_out(Lin_point,tid); // inward local movement
                // switch_out function returns a new Lin_point
                // which is the next point of the former Lin_point
            }
            else
            {
                ++Lin_point;
            }
        }
/*
        stopTime = std::clock();
        elapsedTime = (stopTime - startTime) / double(CLOCKS_PER_SEC);
        std::cout << tid << "bucle2: "<<elapsedTime << std::endl;
        startTime = std::clock();
*/
        clean_Lout(tid); // eliminate Lout redundant points
/*
        stopTime = std::clock();
        elapsedTime = (stopTime - startTime) / double(CLOCKS_PER_SEC);
        std::cout << tid << "cleanLout: " << elapsedTime << std::endl;
        startTime = std::clock();
*/
    }
    //startTime = std::clock();

    Lout.collectList(&Splited_Lout,sublistHead, sublistHeadPosition,0,numThreads );
    Lin.collectList(&Splited_Lin,sublistHead, sublistHeadPosition,1,numThreads);

    //Recalculate thread private variables
    for(int i=0; i < numThreads; i++)
    {
        if(hasOutwardEvolutionInThread[i])
        {
            hasOutwardEvolution = true;
            break;
        }
    }
    for(int i=0; i < numThreads; i++)
    {
        if(hasInwardEvolutionInThread[i])
        {
            hasInwardEvolution = true;
            break;
        }
    }
    //Collect all counters of each thread
    update_for_means_global();
/*
    stopTime = std::clock();
    elapsedTime = (stopTime - startTime) / double(CLOCKS_PER_SEC);
    std::cout << "COLLECT: " << elapsedTime << std::endl;
    startTime = std::clock();
    std::cout << "/////////////////////////////////" << std::endl;

    yi = std::clock();
    elapsedTime = (yi - ya) / double(CLOCKS_PER_SEC);
    std::cout << "BLUCE-ENTERO: " << elapsedTime << std::endl;
*/
    iteration++;
    return;
}

void ActiveContour::do_one_iteration_in_cycle2()
{
    int offset,tid,listSize;
    lists_length = 0;

    Lout.splitList(Splited_Lout,numThreads,sublistHead, sublistHeadPosition,0);
    Lin.splitList(Splited_Lin,numThreads,sublistHead, sublistHeadPosition,1);

    #pragma omp parallel private(tid,offset) reduction(+:listSize)
    {
        listSize=0;
        tid = omp_get_thread_num();

        // scan through Lout with a conditional increment
        for( list<int>::iterator Lout_point = Splited_Lout[tid]->begin(); !Lout_point.end(); )
        {
            offset = *Lout_point;

            if( compute_internal_speed_Fint(offset) > 0 )
            {
                // updates of the variables to calculate the means Cout and Cin
                updates_for_means_in2(offset,tid); // virtual function for region-based models
                Lout_point = switch_in(Lout_point,tid); // outward local movement
                // switch_in function returns a new Lout_point
                // which is the next point of the former Lout_point
            }
            else
            {
                listSize++;
                ++Lout_point;
            }
        }

        clean_Lin(tid); // eliminate Lin redundant points

        // scan through Lin with a conditional increment
        for( list<int>::iterator Lin_point = Splited_Lin[tid]->begin(); !Lin_point.end(); )
        {
            offset = *Lin_point;

            if( compute_internal_speed_Fint(offset) < 0 )
            {
                // updates of the variables to calculate the means Cout and Cin
                updates_for_means_out2(offset,tid); // virtual function for region-based models
                Lin_point = switch_out(Lin_point,tid); // inward local movement
                // switch_out function returns a new Lin_point
                // which is the next point of the former Lin_point
            }
            else
            {
                listSize++;
                ++Lin_point;
            }
        }
        clean_Lout(tid); // eliminate Lout redundant points
    }

    Lout.collectList(&Splited_Lout,sublistHead, sublistHeadPosition,0,numThreads );
    Lin.collectList(&Splited_Lin,sublistHead, sublistHeadPosition,1,numThreads);

    lists_length = listSize;

    //Collect all counters of each thread
    update_for_means_global();

    iteration++;
    return;
}

ActiveContour& ActiveContour::operator++()
{

    // Fast Two Cycle algorithm
    while( !isStopped )
    {

        ////////   cycle 1 : Na_max times, data dependant evolution   ////////
        while( Na < Na_max && !hasLastCycle2 )
        {
            do_one_iteration_in_cycle1();
            calculate_stopping_conditions1(); // it computes hasLastCycle2
            Na++;
            return *this; // just one iteration is performed
        }
        //////////////////////////////////////////////////////////////////////

        if( hasCycle2 )
        {
            ////   cycle 2 : Ns_max times, regularization of the active contour   ////
            while( Ns < Ns_max )
            {
                do_one_iteration_in_cycle2();
                Ns++;
                return *this; // just one iteration is performed
            }
            //////////////////////////////////////////////////////////////////////////

            if( hasLastCycle2 ) // a last cycle 2 has been performed before
            {
                isStopped = true;
            }
            else
            {
                calculate_stopping_conditions2(); // it computes isStopped
            }
        }
        else
        {
            if( hasLastCycle2 )
            {
                isStopped = true;
            }
        }

        Na = 0;
        Ns = 0;
    }

    return *this;
}

void ActiveContour::evolve()
{

    // Fast Two Cycle algorithm

    while( !isStopped )
    {

        ////////   cycle 1 : Na_max times, data dependant evolution   ////////
        while( Na < Na_max && !hasLastCycle2 )
        {
            do_one_iteration_in_cycle1();
            calculate_stopping_conditions1(); // it computes hasLastCycle2
            Na++;
        }
        //////////////////////////////////////////////////////////////////////

        if( hasCycle2 )
        {
            ////   cycle 2 : Ns_max times, regularization of the active contour   ////
            while( Ns < Ns_max )
            {
                do_one_iteration_in_cycle2();
                Ns++;
            }
            //////////////////////////////////////////////////////////////////////////

            if( hasLastCycle2 ) // a last cycle 2 has been performed before
            {
                isStopped = true;
            }
            else
            {
                calculate_stopping_conditions2(); // it computes isStopped
            }
        }
        else
        {
            if( hasLastCycle2 )
            {
                isStopped = true;
            }
        }

        Na = 0;
        Ns = 0;
    }

    return;
}

void ActiveContour::add_Rout_neighbor_to_Lout(int neighbor_offset,int tid)
{
    bool flag = false;
    // if a neighbor ∈ Rout
    if( phi[neighbor_offset] == 3 ) // exterior value
    {
        #pragma omp critical
        {
            if( phi[neighbor_offset] == 3 ) // exterior value
            {
                phi[neighbor_offset] = 1; // outside boundary value
                flag = true;
            }
        }
        if(flag)
        {
            Splited_Lout[tid]->push_front(neighbor_offset);
        }
    }

    return;
}

void ActiveContour::add_Rin_neighbor_to_Lin(int neighbor_offset,int tid)
{
    bool flag = false;
    // if a neighbor ∈ Rin
    if( phi[neighbor_offset] == -3 ) // interior value
    {    
       #pragma omr critical
       {
            if( phi[neighbor_offset] == -3 ) // interior value
            {
                phi[neighbor_offset] = -1; // inside boundary value
                flag = true;
            }
       }
        if(flag)
        {
            Splited_Lin[tid]->push_front(neighbor_offset);
        }
    }

    return;
}

list<int>::iterator ActiveContour::switch_in(list<int>::iterator Lout_point, int tid)
{
    int offset, x, y;
    offset = *Lout_point;
    find_xy_position(offset,x,y); // x and y passed by reference

    // Outward local movement

    #ifndef VERSION_8_CONNECTED_NEIGHBORHOOD

    //==========   4-connected neighborhood   =========
    if( y > 0 )
    {
        add_Rout_neighbor_to_Lout( find_offset(x,y-1), tid );
    }
    if( x > 0 )
    {
        add_Rout_neighbor_to_Lout( find_offset(x-1,y), tid );
    }
    if( x < img_width-1 )
    {
        add_Rout_neighbor_to_Lout( find_offset(x+1,y), tid );
    }
    if( y < img_height-1 )
    {
        add_Rout_neighbor_to_Lout( find_offset(x,y+1), tid );
    }
    //=================================================

    #else

    //==========   8-connected neighborhood   =========
    // if not in the image's border, no neighbors' tests
    if( x > 0 && x < img_width-1 && y > 0 && y < img_height-1 )
    {
        // scan through a 8-connected neighborhood
        for( int dy = -1; dy <= 1; dy++ )
        {
            for( int dx = -1; dx <= 1; dx++ )
            {
                if( !( dx == 0 && dy == 0 ) )
                {
                    add_Rout_neighbor_to_Lout( find_offset(x+dx,y+dy), tid );
                }
            }
        }
    }
    // if in the border, neighbors' tests
    else
    {
        // scan through a 8-connected neighborhood
        for( int dy = -1; dy <= 1; dy++ )
        {
            for( int dx = -1; dx <= 1; dx++ )
            {
                if( !( dx == 0 && dy == 0 ) )
                {
                    // existence tests
                    if( x+dx >= 0 && x+dx < img_width && y+dy >= 0 && y+dy < img_height )
                    {
                        add_Rout_neighbor_to_Lout( find_offset(x+dx,y+dy), tid );
                    }
                }
            }
        }
    }
    //=================================================

    #endif

    phi[offset] = -1; // 1 ==> -1

    int value = *Lout_point;

    Lout_point = Splited_Lout[tid]->erase(Lout_point);

    Splited_Lin[tid]->push_front(value);

    return Lout_point;
}

list<int>::iterator ActiveContour::switch_out(list<int>::iterator Lin_point,int tid)
{
    int offset, x, y;
    offset = *Lin_point;
    find_xy_position(offset,x,y); // x and y passed by reference

    // Inward local movement

    #ifndef VERSION_8_CONNECTED_NEIGHBORHOOD

    //==========   4-connected neighborhood   =========
    if( y > 0 )
    {
        add_Rin_neighbor_to_Lin( find_offset(x,y-1), tid );
    }
    if( x > 0 )
    {
        add_Rin_neighbor_to_Lin( find_offset(x-1,y), tid );
    }
    if( x < img_width-1 )
    {
        add_Rin_neighbor_to_Lin( find_offset(x+1,y), tid );
    }
    if( y < img_height-1 )
    {
        add_Rin_neighbor_to_Lin( find_offset(x,y+1), tid );
    }
    //=================================================

    #else

    //==========   8-connected neighborhood   =========
    // if not in the image's border, no neighbors' tests
    if( x > 0 && x < img_width-1 && y > 0 && y < img_height-1 )
    {
        // scan through a 8-connected neighborhood
        for( int dy = -1; dy <= 1; dy++ )
        {
            for( int dx = -1; dx <= 1; dx++ )
            {
                if( !( dx == 0 && dy == 0 ) )
                {
                    add_Rin_neighbor_to_Lin( find_offset(x+dx,y+dy), tid );
                }
            }
        }
    }
    // if in the border, neighbors' tests
    else
    {
        // scan through a 8-connected neighborhood
        for( int dy = -1; dy <= 1; dy++ )
        {
            for( int dx = -1; dx <= 1; dx++ )
            {
                if( !( dx == 0 && dy == 0 ) )
                {
                    // existence tests
                    if( x+dx >= 0 && x+dx < img_width && y+dy >= 0 && y+dy < img_height )
                    {
                        add_Rin_neighbor_to_Lin( find_offset(x+dx,y+dy), tid );
                    }
                }
            }
        }
    }
    //=================================================

    #endif

    phi[offset] = 1; // -1 ==> 1

    int value = *Lin_point;

    Lin_point = Splited_Lin[tid]->erase(Lin_point);

    Splited_Lout[tid]->push_front(value);

    return Lin_point;
}

int ActiveContour::compute_internal_speed_Fint(int offset)
{
    int x, y;
    find_xy_position(offset,x,y); // x and y passed by reference

    int Fint = 0;

    // if not in the image's border, no neighbors' tests
    if(    x > kernel_radius-1 && x < img_width-kernel_radius
        && y > kernel_radius-1 && y < img_height-kernel_radius
      )
    {
        for( int dy = -kernel_radius; dy <= kernel_radius; dy++ )
        {
            for( int dx = -kernel_radius; dx <= kernel_radius; dx++ )
            {
                Fint += gaussian_kernel[ (kernel_radius+dx)+(kernel_radius+dy)*kernel_length ]
                        * signum_function( -phi[find_offset(x+dx,y+dy)] );
            }
        }
    }
    // if in the border of the image, tests of neighbors
    else
    {
        for( int dy = -kernel_radius; dy <= kernel_radius; dy++ )
        {
            for(int dx = -kernel_radius; dx <= kernel_radius; dx++ )
            {
                if( x+dx >= 0 && x+dx < img_width && y+dy >= 0 && y+dy < img_height )
                {
                    Fint += gaussian_kernel[ (kernel_radius+dx)+(kernel_radius+dy)*kernel_length ]
                            * signum_function( -phi[find_offset(x+dx,y+dy)] );
                }
                else
                {
                    Fint += gaussian_kernel[ (kernel_radius+dx)+(kernel_radius+dy)*kernel_length ]
                            * signum_function( -phi[offset] );
                }
            }
        }
    }

    return Fint;
}

bool ActiveContour::isRedundantLinPoint(int offset) const
{
    int x, y;
    find_xy_position(offset,x,y); // x and y passed by reference

    #ifndef VERSION_8_CONNECTED_NEIGHBORHOOD

    //==========   4-connected neighborhood   =========
    // if ∃ a neighbor ∈ Lout | ∈ Rout
    if( y > 0 )
    {
        if( phi[ find_offset(x,y-1) ] >= 0 )
        {
            return false; // is not redundant point of Lin
        }
    }
    if( x > 0 )
    {
        if( phi[ find_offset(x-1,y) ] >= 0 )
        {
            return false; // is not redundant point of Lin
        }
    }
    if( x < img_width-1 )
    {
        if( phi[ find_offset(x+1,y) ] >= 0 )
        {
            return false; // is not redundant point of Lin
        }
    }
    if( y < img_height-1 )
    {
        if( phi[ find_offset(x,y+1) ] >= 0 )
        {
            return false; // is not redundant point of Lin
        }
    }
    //=================================================

    #else

    //==========   8-connected neighborhood   =========
    // if not in the image's border, no neighbors' tests
    if( x > 0 && x < img_width-1 && y > 0 && y < img_height-1 )
    {
        // scan through a 8-connected neighborhood
        for( int dy = -1; dy <= 1; dy++ )
        {
            for( int dx = -1; dx <= 1; dx++ )
            {
                if( !( dx == 0 && dy == 0 ) )
                {
                    // if ∃ a neighbor ∈ Lout | ∈ Rout
                    if( phi[ find_offset(x+dx,y+dy) ] >= 0 )
                    {
                        return false; // is not redundant point of Lin
                    }
                }
            }
        }
    }
    // if in the border of the image, tests of neighbors
    else
    {
        // scan through a 8-connected neighborhood
        for( int dy = -1; dy <= 1; dy++ )
        {
            for( int dx = -1; dx <= 1; dx++ )
            {
                if( !( dx == 0 && dy == 0 ) )
                {
                    // neighbors tests
                    if( x+dx >= 0 && x+dx < img_width && y+dy >= 0 && y+dy < img_height )
                    {
                        // if ∃ a neighbor ∈ Lout | ∈ Rout
                        if( phi[ find_offset(x+dx,y+dy) ] >= 0 )
                        {
                            return false; // is not redundant point of Lin
                        }
                    }
                }
            }
        }
    }
    //=================================================

    #endif

    // ==> ∀ neighbors ∈ Lin | ∈ Rin
    return true; // is redundant point of Lin
}

bool ActiveContour::isRedundantLoutPoint(int offset) const
{
    int x, y;
    find_xy_position(offset,x,y); // x and y passed by reference

    #ifndef VERSION_8_CONNECTED_NEIGHBORHOOD

    //==========   4-connected neighborhood   =========
    // if ∃ a neighbor ∈ Lin | ∈ Rin
    if( y > 0 )
    {
        if( phi[ find_offset(x,y-1) ] < 0 )
        {
            return false; // is not redundant point of Lout
        }
    }
    if( x > 0 )
    {
        if( phi[ find_offset(x-1,y) ] < 0 )
        {
            return false; // is not redundant point of Lout
        }
    }
    if( x < img_width-1 )
    {
        if( phi[ find_offset(x+1,y) ] < 0 )
        {
            return false; // is not redundant point of Lout
        }
    }
    if( y < img_height-1 )
    {
        if( phi[ find_offset(x,y+1) ] < 0 )
        {
            return false; // is not redundant point of Lout
        }
    }
    //=================================================

    #else

    //==========   8-connected neighborhood   =========
    // if not in the image's border, no neighbors' tests
    if( x > 0 && x < img_width-1 && y > 0 && y < img_height-1 )
    {
        // scan through a 8-connected neighborhood
        for( int dy = -1; dy <= 1; dy++ )
        {
            for( int dx = -1; dx <= 1; dx++ )
            {
                if( !( dx == 0 && dy == 0 ) )
                {
                    // if ∃ a neighbor ∈ Lin | ∈ Rin
                    if( phi[ find_offset(x+dx,y+dy) ] < 0 )
                    {
                        return false; // is not redundant point of Lout
                    }
                }
            }
        }
    }
    // if in the border of the image, tests of neighbors
    else
    {
        // scan through a 8-connected neighborhood
        for( int dy = -1; dy <= 1; dy++ )
        {
            for( int dx = -1; dx <= 1; dx++ )
            {
                if( !( dx == 0 && dy == 0 ) )
                {
                    // neighbors tests
                    if( x+dx >= 0 && x+dx < img_width && y+dy >= 0 && y+dy < img_height )
                    {
                        // if ∃ a neighbor ∈ Lin | ∈ Rin
                        if( phi[ find_offset(x+dx,y+dy) ] < 0 )
                        {
                            return false; // is not redundant point of Lout
                        }
                    }
                }
            }
        }
    }
    //=================================================

    #endif

    // ==> ∀ neighbors ∈ Lout | ∈ Rout
    return true; // is redundant point of Lout
}

void ActiveContour::clean_Lin(int tid)
{
    int offset;

    // scan through Lin with a conditional increment
    for( list<int>::iterator Lin_point =  Splited_Lin[tid]->begin(); !Lin_point.end(); )
    {
        offset = *Lin_point;

        // if ∀ neighbors ∈ Lin | ∈ Rin
        if( isRedundantLinPoint(offset) )
        {
            phi[offset] = -3; // -1 ==> -3
            Lin_point =  Splited_Lin[tid]->erase(Lin_point); // Lin_point ∈ Lin ==> Lin_point ∈ Rin
            // erase function returns a new Lin_point
            // which is the next point of the former Lin_point

        }
        else
        {
            ++Lin_point;
        }
    }
    return;
}

void ActiveContour::clean_Lout(int tid)
{
    int offset;

    // scan through Lout with a conditional increment
    for( list<int>::iterator Lout_point =  Splited_Lout[tid]->begin(); !Lout_point.end(); )
    {
        offset = *Lout_point;

        // if ∀ neighbors ∈ Lout | ∈ Rout
        if( isRedundantLoutPoint(offset) )
        {
            phi[offset] = 3; // 1 ==> 3
            Lout_point = Splited_Lout[tid]->erase(Lout_point); // Lout_point ∈ Lout ==> Lout_point ∈ Rout
            // erase function returns a new Lout_point
            // which is the next point of the former Lout_point
        }
        else
        {
            ++Lout_point;
        }
    }
    return;
}

int ActiveContour::compute_external_speed_Fd(int offset, int tid)
{
    // this class should never be instantiated

    return -1;
    // always an inward movement in each point of the active contour,
    // this speed is not very discriminant...

    // reimplement a better and data-dependent speed function in a child class
}

// at the end of each iteration in the cycle 1
void ActiveContour::calculate_stopping_conditions1()
{
    if( !hasInwardEvolution && !hasOutwardEvolution )
    {
        hasListsChanges = false;
    }


    if( !hasListsChanges || iteration >= iteration_max )
    {
        hasLastCycle2 = true;
    }

    return;
}

// at the end of the cycle 2
void ActiveContour::calculate_stopping_conditions2()
{
    // if the relative difference of active contour length between two cycle2 is less of 1%
    if(   double( std::abs(previous_lists_length-lists_length) )
        / double(previous_lists_length)
        < 0.01
      )
    {
        oscillations_in_a_row++;
    }
    else
    {
        oscillations_in_a_row = 0;
    }
    // keep last length to compare after
    previous_lists_length = lists_length;

    // if 3 times consecutively
    if( oscillations_in_a_row == 3 )
    {
        hasOscillation = true;
    }

    if( hasOscillation || iteration >= iteration_max )
    {
        isStopped = true;
    }

    return;
}

// to display the active contour position in the standard output
// std::cout << ac << std::endl;
std::ostream& operator<<(std::ostream& os, const ActiveContour& ac)
{
    os << std::endl << " -----------------------" << std::endl;
    os << "         Lout(" << ac.get_iteration() << ")" << std::endl;
    os << " -----------------------" << std::endl;
    os << "  index |     x | y     " << std::endl;
    os << " -----------------------" << std::endl;

    int x, y;

    int index = 0;
    for( ofeli::list<int>::const_iterator Lout_point = ac.get_Lout().begin(); !Lout_point.end(); ++Lout_point )
    {
        ac.find_xy_position(*Lout_point,x,y); // x and y passed by reference

        os.width(7);
        os << std::right << index << " |";

        os.width(6);
        os << std::right << x << " | " << y << std::endl;
        index++;
    }
    os << " -----------------------" << std::endl;


    os << std::endl << " -----------------------" << std::endl;
    os << "         Lin(" << ac.get_iteration() << ")" << std::endl;
    os << " -----------------------" << std::endl;
    os << "  index |     x | y     " << std::endl;
    os << " -----------------------" << std::endl;

    index = 0;
    for( ofeli::list<int>::const_iterator Lin_point = ac.get_Lin().begin(); !Lin_point.end(); ++Lin_point )
    {
        ac.find_xy_position(*Lin_point,x,y); // x and y passed by reference

        os.width(7);
        os << std::right << index << " |";

        os.width(6);
        os << std::right << x << " | " << y << std::endl;
        index++;
    }
    os << " -----------------------" << std::endl;

    return os;
}

void ActiveContour::display() const
{
    std::cout << *this;

    return;
}



void ActiveContour::calculate_means()
{
    return;
}

void ActiveContour::updates_for_means_in1(int tid)
{
    return;
}

void ActiveContour::updates_for_means_out1(int tid)
{
    return;
}

void ActiveContour::updates_for_means_in2(int offset,int tid)
{
    return;
}

void ActiveContour::updates_for_means_out2(int offset,int tid)
{
    return;
}
void ActiveContour::update_for_means_global()
{
    return;
}


void ActiveContour::initialize_for_each_frame(const unsigned char* img_data1)
{
    if( img_data1 == NULL )
    {
        std::cerr << std::endl << "The pointer img_data1 must be a non-null pointer, it must be allocated." << std::endl;
    }

    img_data = img_data1;

    initialize_for_each_frame();

    return;
}



const char* ActiveContour::get_phi() const
{
    return phi;
}

const ofeli::list<int>& ActiveContour::get_Lout() const
{
    return Lout;
}

const ofeli::list<int>& ActiveContour::get_Lin() const
{
    return Lin;
}



bool ActiveContour::get_hasListsChanges() const
{
    return hasListsChanges;
}

bool ActiveContour::get_hasOscillation() const
{
    return hasOscillation;
}

int ActiveContour::get_iteration() const
{
    return iteration;
}

int ActiveContour::get_iteration_max() const
{
    return iteration_max;
}

bool ActiveContour::get_isStopped() const
{
    return isStopped;
}

void ActiveContour::calculateCovering() const
{
    int sum=0,size=img_width*img_height;

    #pragma omp parallel for reduction(+:sum)
        for(int i= 0; i < size; i++)
        {
           if(phi[i] > 0 )
           {
                sum++;
           }
        }
    double percentage = (sum * 100) / size;
    std::cout << "Recubrimiento: " << percentage << "%" << std::endl;
}


std::vector<std::vector<int>* >* ActiveContour::isolateIsland()
{
    int cont=0;
    int size=img_width*img_height, i=0;
    list<int>* islandPoints = new list<int>();
    bool* visitedIslands = new bool[size];

    for(int a=0; a < size; a++) visitedIslands[a] = false;

    //Find a island
    for(; i < size; i++)
    {
        if(phi[i] == -1) break;
    }

    int previousPoint=-1,x,y,min_x=99999999,min_y=99999999,max_x=-1,max_y=-1;
    bool stop = false;

    while(1)
    {
        find_xy_position(i,x,y);

        //Take borders to make the image later
        if(x < min_x) min_x = x;
        if(x > max_x) max_x = x;
        if(y < min_y) min_y = y;
        if(y > max_y) max_y = y;

        int previousPosition = i;

        //No se comprueba si está fuera de la imagen
        for(int dx=-1; dx <= 1; dx++)
        {
            for(int dy=-1; dy <= 1; dy++)
            {
                int offset = find_offset(x+dx,y+dy);

                //Test for the next point of the border of the island
                if(phi[offset] == -1 && previousPoint != offset && offset != i)
                {
                    stop=true;

                    //Save the point
                    previousPoint= i;
                    islandPoints->push_front(i);
                    visitedIslands[i] = true;

                    //Continue with the next point
                    i = offset;
                }
                if(stop) break;
            }
            if(stop)break;
        }
        stop = false;

        //When we've make a round or the island
        if(visitedIslands[i] || previousPosition==i) break;
        cont++;

    }//end while

    //Create the isolated island
    std::vector< std::vector<int>* >* island = new std::vector< std::vector<int>* >(max_x-min_x +1);
    for(int i=0; i < island->size(); i++)
       island->at(i) = new std::vector<int>(max_y-min_y+1);

    //Fill up the matrix
    for(int a=0; a < island->size(); a++)
        for(int b=0; b < island->at(a)->size(); b++)
            island->at(a)->at(b) = 0;

    for(list<int>::iterator position = islandPoints->begin(); !position.end(); ++position)
    {
        int X,Y;
        find_xy_position(*position,X,Y);

        //Reajust the positions
        X = X - min_x;
        Y = Y - min_y;
        island->at(X)->at(Y) = 1;
    }

    /*for(int a=0; a < island->at(0)->size(); a++)
    {
        for(int b=0; b < island->size(); b++)
        {
            std::cout << ", " << island->at(b)->at(a);
        }std::cout << std::endl;
    }*/

    return island;
}


}
