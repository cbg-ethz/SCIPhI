/**
 * SCIPhI: Single-cell mutation identification via phylogenetic inference
 * <p>
 * Copyright (C) 2018 ETH Zurich, Jochen Singer
 * <p>
 * This file is part of SCIPhI.
 * <p>
 * SCIPhI is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * SCIPhI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with SCIPhI. If not, see <http://www.gnu.org/licenses/>.
 *
 * @author: Jochen Singer
 */

#ifndef RAND_H
#define RAND_H

#include <vector>
#include "rand.h"
#include <iostream>
#include <random>
#include <stdlib.h>
#include <time.h>

using namespace std;


/*****    functions for sampling random numbers inside C++  *****/
inline 
void 
initRand(){
	time_t t;
	time(&t);
	srand((unsigned int)t);              // initialize random number generator
	//srand(1);
}

inline
int
sampleRandomMove(std::array<double, 4> const & prob)
{ 

    double percent = rand() % 100;
    double sum = prob[0];
    for(std::size_t i=0; i<prob.size()-1; i++){
        if(percent <= sum*100){
          return i+1;
        }
        sum += prob[i+1];
    }
    return prob.size();
}

#endif
