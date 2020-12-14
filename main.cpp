/**
 * @file main.cpp
 * @author Alex Lu (luh60@mcmaster.ca)
 * @brief An Nbody simulation that applies a version of Barnes-Hut Tree code and the Velocity-Verlet algorithm 
 * @version 0.1
 * @date 2020-12-12
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <chrono>
#include <cmath>
#include <fstream>

#include "bodygen.hpp"


using namespace std;

/**
 * @brief The main function. Runs the Spacetree and bodygen constructors and the simulate function based on inputs. Also checks for correct inputs and returns error messages if command line inputs are incorrect.
 * 
 * @param argc Number of inputs - must be 3 or 4 (in addition to ./bodygen.exe)
 * @param argv The inputs
 * @return int 
 */
int main(int argc, char* argv[])
{
    chrono::time_point start_time{chrono::steady_clock::now()};
    if(argc > 5 or argc < 4)
    {
        std::cout << "Incorrect number of inputs\n";
        return 0;
    }
    else if(argc == 4)
    {
        for(size_t j{2}; j < 4; j++)
        {
            if((size_t) atoi(argv[j]) <= 0)
            {
                std::cout << "Invalid inputs detected - please input a positive integer, long double and postive integer, or a string, long double, and positive integer.\n";
                return 0;
            }
        }
        if((size_t) atoi(argv[1]) == 0)
        {
            ifstream file(argv[1]);
            if(!file.is_open())
            {
                std::cout << "Input file not found\n";
                return 0;
            }
            else
            {
                string str = (string) argv[1];
                long double ld = (long double) atoi(argv[2]);
                size_t st = (size_t) atoi(argv[3]);
                bodygen gen{str,ld,st};
                gen.simulate();
                chrono::time_point end_time{chrono::steady_clock::now()};
                chrono::duration<double> elapsed_time_seconds{end_time - start_time};
                chrono::duration<double, milli> elapsed_time_milli{end_time - start_time};
                cout << "Elapsed time: " << elapsed_time_seconds.count() << " seconds, ";
                return 0;
            }
        }
    }
    else
    {
        for(size_t j{1}; j < 5; j++)
        {
            if((size_t) atoi(argv[j]) <= 0 && j != 2)
            {
                std::cout << "Invalid inputs detected - please input a positive integer, long double and postive integer, or a string, long double, and positive integer.\n";
                return 0;
            }
        }
        size_t st1 = (size_t) atoi(argv[1]);
        string str = (string) argv[2];
        long double ld = (long double) atoi(argv[3]);
        size_t st = (size_t) atoi(argv[4]);
        bodygen gen{st1,str,ld,st};
        gen.simulate();
        chrono::time_point end_time{chrono::steady_clock::now()};
        chrono::duration<double> elapsed_time_seconds{end_time - start_time};
        chrono::duration<double, milli> elapsed_time_milli{end_time - start_time};
        cout << "Elapsed time: " << elapsed_time_seconds.count() << " seconds, ";
        return 0;
    }
}