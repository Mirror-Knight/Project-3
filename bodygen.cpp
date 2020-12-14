/**
 * @file bodygen.cpp
 * @author Alex Lu (luh60@mcmaster.ca)
 * @brief An Nbody simulation that applies a version of Barnes-Hut Tree code and the Velocity-Verlet algorithm 
 * @version 0.1
 * @date 2020-12-12
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <chrono>
#include <fstream>
#include <direct.h>
#include <random>

#include "bodygen.hpp"

using namespace std;

/**
 * @brief Global variable declared here. Unconventional, but in the context of a physical simulation unchanging physical constants should be constant.
 * 
 */
long double G = 6.67E-11;

/**
 * @brief Overloaded operator + that adds the elements of two arrays to produce a third array
 * 
 * @tparam T 
 * @tparam L 
 * @param v1 Input array 1
 * @param v2 Input array 2
 * @return array<T,L> 
 */
template <typename T, size_t L>
array<T,L> operator+(const array<T,L> &v1, const array<T,L> &v2)
{
    array<T,L> sum;
    for(size_t i{0}; i < L; i++)
    {
        sum[i] = v1[i] + v2[i];
    }
    return(sum);
}

/**
 * @brief Overloaded operator - that subtracts the elements of two arrays to produce a third array
 * 
 * @tparam T 
 * @tparam L 
 * @param v1 Input array 1
 * @param v2 Input array 2
 * @return array<T,L> 
 */
template <typename T, size_t L>
array<T,L> operator-(const array<T,L> &v1, const array<T,L> &v2)
{
    array<T,L> sum;
    for(size_t i{0}; i < L; i++)
    {
        sum[i] = v1[i] - v2[i];
    }
    return(sum);
}

/**
 * @brief Overloaded operator * that takes the dot product of two arrays
 * 
 * @tparam T 
 * @tparam L 
 * @param v1 Input array 1
 * @param v2 Input array 2
 * @return array<T,L> 
 */
template <typename T, size_t L>
T operator*(const array<T,L> &v1, const array<T,L> &v2)
{
    T sum{0};
    for(size_t i{0}; i < L; i++)
    {
        sum = sum + v1[i]*v2[i];
    }
    return(sum);
}

/**
 * @brief Overloaded operator * that multiples an array by a scalar
 * 
 * @tparam T 
 * @tparam L 
 * @param s Input scalar (must be of same type as the elements in the array)
 * @param v Input array
 * @return array<T,L> 
 */
template <typename T, size_t L>
array<T,L> operator*(const T &s, array<T,L> &v)
{
    array<T,L> returnval;
    for(size_t i{0}; i < L; i++)
    {
        returnval[i] = s*v[i];
    }
    return(returnval);
}

/**
 * @brief Overloaded operator >= that compares two arrays and determines if each element in the first array is greater than or equal to the corresponding element in the second array
 * 
 * @tparam T 
 * @tparam L 
 * @param v1 Input array 1
 * @param v2 Input array 1
 * @return true Returns true if all comparisons satisfy the condition
 * @return false Returns false if at least one comparison does not satisfy the condition
 */
template <typename T, size_t L>
bool operator>=(const array<T,L> &v1, const array<T,L> &v2)
{
    for(size_t i{0}; i < L; i++)
    {
        if(v1[i] >= v2[i])
        {
            continue;
        }
        else
        {
            return(false);
        }  
    }
    return(true);
}

/**
 * @brief Overloaded operator < that compares two arrays and determines if each element in the first array is less than the corresponding element in the second array
 * 
 * @tparam T 
 * @tparam L 
 * @param v1 Input array 1
 * @param v2 Input array 1
 * @return true Returns true if all comparisons satisfy the condition
 * @return false Returns false if at least one comparison does not satisfy the condition
 */
template <typename T, size_t L>
bool operator<(const array<T,L> &v1, const array<T,L> &v2)
{
    for(size_t i{0}; i < L; i++)
    {
        if(v1[i] < v2[i])
        {
            continue;
        }
        else
        {
            return(false);
        }  
    }
    return(true);
}

/**
 * @brief Overloaded operator == that compares two strings and determine if they are equal
 * 
 * @param a Input string 1
 * @param b Input string 2
 * @return true 
 * @return false 
 */
bool operator==(const string &a, const string &b)
{
    for(size_t i{0}; i < a.size(); i++)
    {
        if(a[i] != b[i])
        {
            return(false);
        }
    }
    return(true);
}
/**
 * @brief Overloaded operator << that writes out body data - position, velocity, mass and radius. 
 * 
 * @param out ostream object
 * @param b body
 * @return ostream& 
 */
ostream &operator<<(ostream &out, body &b)
{
    out << b.index << ',';
    for(size_t i{0}; i < 3; i++)
    {
        out << b.position[i] << ',';
    }
    for(size_t j{0}; j < 3; j++)
    {
        out << b.velocity[j] << ',';
    }
    out << b.mass << ',' << b.radius;
    return(out);
}

/**
 * @brief Overloaded operator << that writes out data for a vector of bodies - positions and radii. Used in a different context than above
 * 
 * @param out ostream object
 * @param v vector of bodies
 * @return ostream& 
 */
ostream &operator<<(ostream &out, vector<body> &v)
{
    out << "x coord" << ',' << "y coord" << ',' << "z coord" << ',' << "scalar\n";
    for(size_t i{0}; i < v.size(); i++)
    {
        out << v[i].position[0] << ',' << v[i].position[1] << ',' << v[i].position[2] << ',' << v[i].radius << '\n';
    }
    return(out);
}

/**
 * @brief Calculates the modulus of a vector (in this case array)
 * 
 * @tparam T type
 * @tparam L size_t
 * @param v input array
 * @return T resulting modulus
 */
template <typename T, size_t L>
T moodulus(const array<T,L> &v)
{
    T sum{0};
    for(size_t i{0}; i < v.size(); i++)
    {
        sum = sum + v[i]*v[i];
    }
    return(sqrt(sum));
}

/**
 * @brief returns the sign of some value
 * 
 * @tparam T 
 * @param val input value
 * @return int 1 or -1
 */
template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

/**
 * @brief returns cross product of two arrays of length 3
 * 
 * @tparam T 
 * @param v1 input array 1
 * @param v2 input arrat 2
 * @return array<T,3> 
 */
template <typename T>
array<T,3> crossprod(array<T,3> &v1, array<T,3> &v2)
{
    array<T,3> result = {v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]};
    return(result);
}



/**
 * @brief Construct a new Spacetree:: Spacetree object
 * 
 * @param inputreg Initializes private member regi to this
 */
Spacetree::Spacetree(region inputreg)
    : regi{inputreg} {}

/**
 * @brief Makes a tree given the input region regi
 * 
 * @return Node* returns the tree
 */
Node* Spacetree::treegen()
{
    Node* root = new Node;
    regi.checkcol = false;
    root = makeatree(regi);
    return(root);
}

/**
 * @brief Adds null nodes to a leaf node via a simple loop through all nodes the leaf node leads to
 * 
 * @param newnode Input leaf node
 * @return Node* Returns the leaf node
 */
Node* Spacetree::addnulls(Node* newnode)
{
    for(size_t i{0}; i < 8; i++)
    {
        newnode->Nodelist[i] = NULL;
    }
    return(newnode);
} 

/**
 * @brief Updates the velocities of all bodies in a vector if there are collisions. Two bodies will elastically collide if the distance between them is less than the sum of their radii
 * 
 * @param bodyvector Input vector of bodies. Only two bodies can interact at a given moment.
 * @param indices Stores the used body indices. Bodies with these indices will be ignored 
 * @return vector<body> 
 */
vector<body> Spacetree::updatecollision(vector<body> &bodyvector)
{
    vector<size_t> indices(0);
    for(size_t i{0}; i < bodyvector.size(); i++)
    {
        for(size_t j{0}; j < bodyvector.size(); j++)
        {
            int check{0};
            if(i == j)
            {
                continue;
            }
            for(size_t k{0}; k < indices.size(); k++)
            {
                if(i == indices[k] or j == indices[k])
                {
                    check = 1;
                    break;
                }
            }
            if(check == 1)
            {
                continue;
            }
            if(moodulus(bodyvector[i].position - bodyvector[j].position) < bodyvector[i].radius + bodyvector[j].radius)
            {
                const long double fact1 = (2*bodyvector[j].mass*((bodyvector[i].velocity - bodyvector[j].velocity)*(bodyvector[i].position - bodyvector[j].position))/((bodyvector[i].mass + bodyvector[j].mass)*moodulus(bodyvector[i].position - bodyvector[j].position)*moodulus(bodyvector[i].position - bodyvector[j].position)));
                const long double fact2 = (2*bodyvector[i].mass*((bodyvector[j].velocity - bodyvector[i].velocity)*(bodyvector[j].position - bodyvector[i].position))/((bodyvector[i].mass + bodyvector[j].mass)*moodulus(bodyvector[i].position - bodyvector[j].position)*moodulus(bodyvector[i].position - bodyvector[j].position)));
                array<long double, 3> posdiff1 = bodyvector[i].position - bodyvector[j].position;
                array<long double, 3> posdiff2 = bodyvector[j].position - bodyvector[i].position;
                bodyvector[i].velocity = bodyvector[i].velocity - fact1*posdiff1;
                bodyvector[j].velocity = bodyvector[j].velocity - fact2*posdiff2;
                indices.push_back(i);
                indices.push_back(j);
            }
        }
    }
    return(bodyvector);
}


/**
 * @brief Recursively makes a tree given some input region. Values such as mass, center of gravity, extent from the region are stored in the current node, before the region is subdivided into eight and the function is recursively called on each subdivision. Collisions are also updated here, but only if the extent of a region is 10 times the maximum radius of all bodies in that region
 * 
 * @param reg Input region that stores boundaries and a vector of bodies inside the region
 * @return Node* Returns the tree
 */
Node* Spacetree::makeatree(region reg)
{
    Node* pointer = new Node;
    pointer->nodepath = reg.regnodepath;
    if(reg.bodiesinregion.size() == 1)
    {
        pointer = addnulls(pointer);
        pointer->isleaf = true;
        pointer->solebody = reg.bodiesinregion[0];
        pointer->cog = pointer->solebody.position;
        pointer->cogmass = pointer->solebody.mass;
        pointer->extent = 0;
        return(pointer);
    }
    pointer->isleaf = false;
    
    long double totmass{0};
    long double extent{0};
    array<long double, 3> tempcog = {0,0,0};
    vector<body> regbods = reg.bodiesinregion;
    long double maxrad{0};
    for(size_t h{0}; h < regbods.size(); h++)
    {
        totmass = totmass + regbods[h].mass;
        if(regbods[h].radius > maxrad)
        {
            maxrad = regbods[h].radius;
        }      
    }
    for(size_t hh{0}; hh < regbods.size(); hh++)
    {
        tempcog = tempcog + ((regbods[hh].mass)/(totmass))*regbods[hh].position;
    }
    for(size_t hhh{0}; hhh < regbods.size(); hhh++)
    {
        extent = extent + moodulus(regbods[hhh].position - tempcog);
    }
    
    pointer->cogmass = totmass;
    pointer->cog = tempcog;
    pointer->extent = 2*extent/regbods.size();

    
    
    if(2*extent/regbods.size() < 10*maxrad && !reg.checkcol)
    {
        reg.bodiesinregion = updatecollision(reg.bodiesinregion);
        reg.checkcol = true;
    }
    
    long double xboxlength = (reg.xrange[1] - reg.xrange[0])/2;
    long double yboxlength = (reg.yrange[1] - reg.yrange[0])/2;
    long double zboxlength = (reg.zrange[1] - reg.zrange[0])/2;  

    array<region,8> regionlist; 
    array<string,8> pathnames = {"dll","dlr","dal","dar","ull","ulr","ual","uar"};
    
    regionlist[0].xrange = {reg.xrange[0], reg.xrange[0] + xboxlength};
    regionlist[0].yrange = {reg.yrange[0], reg.yrange[0] + yboxlength};
    regionlist[0].zrange = {reg.zrange[0], reg.zrange[0] + zboxlength};

    regionlist[7].xrange = {reg.xrange[0]+xboxlength, reg.xrange[1]};
    regionlist[7].yrange = {reg.yrange[0]+yboxlength, reg.yrange[1]};
    regionlist[7].zrange = {reg.zrange[0]+zboxlength, reg.zrange[1]};

    regionlist[1].xrange = {reg.xrange[0]+xboxlength, reg.xrange[1]};
    regionlist[1].yrange = {reg.yrange[0], reg.yrange[0] + yboxlength};
    regionlist[1].zrange = {reg.zrange[0], reg.zrange[0] + zboxlength};

    regionlist[2].xrange = {reg.xrange[0], reg.xrange[0] + xboxlength};
    regionlist[2].yrange = {reg.yrange[0]+yboxlength, reg.yrange[1]};
    regionlist[2].zrange = {reg.zrange[0], reg.zrange[0] + zboxlength};

    regionlist[3].xrange = {reg.xrange[0]+xboxlength, reg.xrange[1]};
    regionlist[3].yrange = {reg.yrange[0]+yboxlength, reg.yrange[1]};
    regionlist[3].zrange = {reg.zrange[0], reg.zrange[0] + zboxlength};

    regionlist[4].xrange = {reg.xrange[0], reg.xrange[0] + xboxlength};
    regionlist[4].yrange = {reg.yrange[0], reg.yrange[0] + yboxlength};
    regionlist[4].zrange = {reg.zrange[0]+zboxlength, reg.zrange[1]};

    regionlist[5].xrange = {reg.xrange[0]+xboxlength, reg.xrange[1]};
    regionlist[5].yrange = {reg.yrange[0], reg.yrange[0] + yboxlength};
    regionlist[5].zrange = {reg.zrange[0]+zboxlength, reg.zrange[1]};

    regionlist[6].xrange = {reg.xrange[0], reg.xrange[0] + xboxlength};
    regionlist[6].yrange = {reg.yrange[0]+yboxlength, reg.yrange[1]};
    regionlist[6].zrange = {reg.zrange[0]+zboxlength, reg.zrange[1]};    


    for(size_t i{0}; i < 8; i++)
    {
        vector<body> newbodiesinregion(0);
        regionlist[i].bodiesinregion = newbodiesinregion;
        regionlist[i].regnodepath = reg.regnodepath + pathnames[i];
        array<long double, 3> minranges = {regionlist[i].xrange[0], regionlist[i].yrange[0], regionlist[i].zrange[0]};
        array<long double, 3> maxranges = {regionlist[i].xrange[1], regionlist[i].yrange[1], regionlist[i].zrange[1]};
        for(size_t j{0}; j < reg.bodiesinregion.size(); j++)
        {
            if (reg.bodiesinregion[j].position >= minranges && reg.bodiesinregion[j].position < maxranges)
            {
                regionlist[i].bodiesinregion.push_back(reg.bodiesinregion[j]);
            }
            else
            {
                continue;
            }
        }

        if(regionlist[i].bodiesinregion.size() >= 1)
        {
            if(reg.checkcol)
            {
                regionlist[i].checkcol = true;
            }
            else
            {
                regionlist[i].checkcol = false;
            }
            pointer->Nodelist[i] = makeatree(regionlist[i]);
        }
        else if(regionlist[i].bodiesinregion.size() == 0)
        {
            pointer->Nodelist[i] = NULL;
        }
    }
    return(pointer);
}

/**
 * @brief Construct a new bodygen::bodygen object. This constructor is called if there is an input file.
 * 
 * @param inputstring Initializes filename to inputstring. This is the desired filename to read
 * @param tstepinput Initializes timestep to tstepinput. This is the desired timestep
 * @param iter Initializes iter to iterations. This is the total number of iterations
 */
bodygen::bodygen(const string inputstring, long double tstepinput, const size_t iter)
    : filename{inputstring}, timestep{tstepinput}, iterations{iter} 
{
    bodyvector.resize(0);
    writeinitfile = false;
}

/**
 * @brief Construct a new bodygen::bodygen object. This constructor is called if there is no input file.
 * 
 * @param inputcount Initializes count to inputcount. This is the desired number of bodies
 * @param st Initializes filename to st. This is the desired filename to write the generate initial data to
 * @param tstepinput Initializes timestep to tstepinput. This is the desired timestep
 * @param iter Initializes iter to iterations. This is the total number of iterations
 */
bodygen::bodygen(const size_t inputcount, const string st, long double tstepinput, const size_t iter)
    : count{inputcount}, filename{st}, timestep{tstepinput}, iterations{iter} 
{                                                                           
    bodyvector.resize(0);
    writeinitfile = true;
}

/**
 * @brief Recursively deletes a tree, and frees up the memory occupied by the tree
 * 
 * @param tree Input tree to delete
 */
void bodygen::deletetree(Node* tree)
{
    if(tree == NULL)
    {
        return;
    }
    for(size_t i{0}; i < 8; i++)
    {
        deletetree(tree->Nodelist[i]);
    }
    delete tree;
}

/**
 * @brief The main function that does the Nbody simulation. Either a file is read, or data is generated. Updating functions are run to update the bodies before the tree is deleted and remade.
 * 
 */
void bodygen::simulate()
{
    if(writeinitfile)
    {
        datatree = makebodies(); //Make bodies if theres no initial data file
    }
    else
    {
        string infolist[8];
        string radiuss;
        ifstream file(filename); //Opens initial data file
        array<long double, 6> minmax = {0,0,0,0,0,0};
        int previndex{-1};
        while(file.is_open() == true)
        {
            body newbody;
            for(size_t i{0}; i < 8; i++)
            {
                getline(file, infolist[i], ',');
            }
            getline(file, radiuss, '\n');
            if(stoi(infolist[0]) == previndex)
            {
                break;
            }
            newbody.index = stoi(infolist[0]);
            newbody.position = {stold(infolist[1]),stold(infolist[2]),stold(infolist[3])};
            newbody.velocity = {stold(infolist[4]),stold(infolist[5]),stold(infolist[6])};
            newbody.mass = stold(infolist[7]);
            newbody.radius = stold(radiuss);
            newbody.acceleration = {0,0,0};
            newbody.newacceleration = {0,0,0};
            bodyvector.push_back(newbody);
            previndex = newbody.index;

            size_t k2{0};
            for(size_t k{0}; k < 3; k++)
            {
                if(newbody.position[k] < minmax[k2])
                {
                    minmax[k2] = newbody.position[k];
                }
                else if(newbody.position[k] > minmax[k2+1])
                {
                    minmax[k2+1] = newbody.position[k];
                }
                k2 = k2 + 2;
            }
        }
        file.close();
        space.xrange = {minmax[0] - 1,minmax[1] + 1};
        space.yrange = {minmax[2] - 1,minmax[3] + 1};
        space.zrange = {minmax[4] - 1,minmax[5] + 1};
        space.bodiesinregion = bodyvector;
        Spacetree space_tree{space};
        datatree = space_tree.treegen();
    }
    string strdirname = filename.substr(0, filename.size()-4);
    const char* dirname = strdirname.c_str();
    mkdir(dirname);
    array<long double,6> minimaxi;
    size_t j{0};
    int ccount{0};
    for(size_t i{0}; i < iterations; i++)
    {
        datatree = updateallacceleration(datatree, datatree);
        datatree = update(datatree);
        if(j == 100)
        {
            ofstream datafile;
            datafile.precision(30);
            filename = ".\\" + strdirname + "\\" + strdirname + ".csv." + to_string(ccount);
            datafile.open(filename);
            datafile << fixed << bodyvector;
            datafile.close();
            j = 0;
            ccount = ccount + 1;
        }
        minimaxi = calcminmax();
        space.xrange = {minimaxi[0] - 1,minimaxi[1] + 1};
        space.yrange = {minimaxi[2] - 1,minimaxi[3] + 1};
        space.zrange = {minimaxi[4] - 1,minimaxi[5] + 1};
        space.bodiesinregion = bodyvector;
        Spacetree space_tree{space};
        deletetree(datatree);
        datatree = space_tree.treegen();
        j = j + 1;
    }
}

/**
 * @brief Calculates the boundaries of a region
 * 
 * @return array<long double,6> 
 */
array<long double,6> bodygen::calcminmax()
{
    array<long double,6> minmax{0,0,0,0,0,0};
    for(size_t i{0}; i < bodyvector.size(); i++)
    {
        size_t k2{0};
        for(size_t k{0}; k < 3; k++)
        {
            if(bodyvector[i].position[k] < minmax[k2])
            {
                minmax[k2] = bodyvector[i].position[k];
            }
            else if(bodyvector[i].position[k] > minmax[k2+1])
            {
                minmax[k2+1] = bodyvector[i].position[k];
            }
            k2 = k2 + 2;
        }
    }
    return(minmax);
}

/**
 * @brief Randomly generates a set of coordinates within some annulus with inner radius r1 and outer radius r2
 * 
 * @param r1 Inner radius
 * @param r2 Outer radius
 * @return array<long double,2> Returns the coordinates as array
 */
array<long double,2> bodygen::randcircgen(long double r1, long double r2)
{
    random_device rd;
    mt19937_64 mt64(rd());
    long double scale = ((long double) mt64())/((long double) mt64.max());
    long double scale2 = ((long double) mt64())/((long double) mt64.max());
    long double scale3 = ((long double) mt64())/((long double) mt64.max());
    long double chosenr = r1 + (r2 - r1)*scale;
    array<long double, 2> xd = {0,0};
    xd[0] = (2*scale3 - 1)*chosenr;
    xd[1] = sgn(2*scale2 - 1)*sqrt(chosenr*chosenr - xd[0]*xd[0]);
    return(xd);
}

/**
 * @brief Randomly generates a set of coordinates within some spherical shell with inner radius r1 and outer radius r2
 * 
 * @param r1 Inner radius
 * @param r2 Outer radius
 * @return array<long double,3> Returns coordinates of arrays
 */
array<long double,3> bodygen::randspheregen(long double r1, long double r2)
{
    random_device rd;
    mt19937_64 mt64(rd());
    long double scale = ((long double) mt64())/((long double) mt64.max());
    long double scale2 = ((long double) mt64())/((long double) mt64.max());
    long double scale3 = ((long double) mt64())/((long double) mt64.max());
    long double scale4 = ((long double) mt64())/((long double) mt64.max());
    long double chosenr = r1 + (r2 - r1)*scale;
    array<long double, 3> xd = {0,0,0};
    xd[0] = (2*scale2 - 1)*chosenr;
    xd[1] = (2*scale3 - 1)*sqrt(chosenr*chosenr - xd[0]*xd[0]);
    xd[2] = sgn(2*scale4 - 1)*sqrt(chosenr*chosenr - xd[0]*xd[0] - xd[1]*xd[1]);
    return(xd);
}

/**
 * @brief Generates body data. randpheregen and randcircgen are used to help generate random body data when needed. The uncommented code randomly generates bodies in some cube, while the commented code more elaborate generates bodies within some sphere or annulus with appropriate initial velocites such that they orbit, using the cross product. The body data is then placed into a region, then placed into a tree
 * 
 * @return Node* Returns the tree
 */
Node* bodygen::makebodies()
{
    bodyvector.resize(count);
    ofstream datafile;
    datafile.precision(30);
    datafile.open(filename);
    array<long double, 8> randlist;
    for(size_t i{0}; i < count; i++)
    {
        for(size_t j{0}; j < 8; j++)
        {
            random_device rd;
            mt19937_64 mt64(rd());
            randlist[j] = 2*((long double) mt64()/((long double) mt64.max())) - 1;
        }
        //array<long double,3> spherevars = randspheregen(10, 1E7);
        //array<long double,2> circvars = randcircgen(1.2E9, 4E9);

        //bodyvector[i].position = {spherevars[0], spherevars[1], spherevars[2]};
        //bodyvector[i].position = {circvars[0], circvars[1], 1E5*randlist[2]};

        //array<long double,3> randrotvec = {randlist[0],randlist[1],randlist[2]};
        //randrotvec = (1/moodulus(randrotvec))*randrotvec;
        //array<long double,3> perpvec = crossprod(bodyvector[i].position, randrotvec);
        //perpvec = (1/moodulus(perpvec))*perpvec;
        //const long double extrafactor = sqrt(G*1E30/(moodulus(bodyvector[i].position)))/moodulus(bodyvector[i].position);

        //bodyvector[i].velocity = crossprod(bodyvector[i].position,perpvec);
        //bodyvector[i].velocity = extrafactor*bodyvector[i].velocity;
        bodyvector[i].position = {10E15*randlist[0], 10E15*randlist[1], 10E15*randlist[2]};
        bodyvector[i].velocity = {1000*randlist[3], 1000*randlist[4], 1000*randlist[5]};
        bodyvector[i].mass = 3E30*(randlist[6]+1)/2;
        bodyvector[i].radius = 1E9*(randlist[7]+1)/2;
        bodyvector[i].index = i;


        if(i == count - 1)
        {
            /*
            bodyvector[i].position = {0,0,0};
            bodyvector[i].velocity = {0,0,0};
            bodyvector[i].mass = 1E40;
            bodyvector[i].radius = 0;
            bodyvector[i].index = i;
            */
            datafile << fixed << bodyvector[i];
        }
        else
        {
            datafile << fixed << bodyvector[i] << '\n';
        }
        
    }
    datafile.close();
    array<long double,6> minimaxi = calcminmax();
    space.xrange = {minimaxi[0] - 1,minimaxi[1] + 1};
    space.yrange = {minimaxi[2] - 1,minimaxi[3] + 1};
    space.zrange = {minimaxi[4] - 1,minimaxi[5] + 1};
    space.bodiesinregion = bodyvector;
    Spacetree space_tree{space};
    return(space_tree.treegen());
}

/**
 * @brief Determines whether or not some node "root" is a child of some node "tree" by checking the nodepath of each node. Function is recursively called through all child nodes of "tree"
 * 
 * @param root Input Node
 * @param tree Input Node
 * @return true 
 * @return false 
 */
bool bodygen::comparetree(Node* root, Node* tree)
{
    if(root->nodepath == tree->nodepath)
    {
        return(true);
    }
    else
    {
        for(size_t i{0}; i < 8; i++)
        {
            if(tree->Nodelist[i] != NULL)
            {
                return(comparetree(root,tree->Nodelist[i]));
            }
        }
        return(false); 
    }
    
}

/**
 * @brief Updates all the accelerations all leafs of the tree recursively. Takes in two inputs since a function updatesingleacceleration is called which needs a leaf node and the whole tree as inputs
 * 
 * @param tree Input node
 * @param wholetree Input node
 * @return Node* 
 */
Node* bodygen::updateallacceleration(Node* tree, Node* wholetree)
{
    if(tree == NULL)
    {
        return(tree);
    }
    if(tree->isleaf)
    {
        tree = updatesingleacceleration(tree,wholetree);
    }
    else
    {
        for(size_t i{0}; i < 8; i++)
        {
            tree->Nodelist[i] = updateallacceleration(tree->Nodelist[i],wholetree);
        }
    }
    return(tree);
}

/**
 * @brief Updates the acceleration of a single leaf node "root", recursively from the whole tree. If some node is far enough away and bodies span a small enough angular size, the leaf node sees a single mass at that node's center of gravity instead of a collection of bodies.
 * 
 * @param root Input leaf node
 * @param tree Input tree
 * @return Node* 
 */
Node* bodygen::updatesingleacceleration(Node* root, Node* tree)
{
    if(tree == NULL)
    {
        return(root);
    }
    if(!comparetree(root,tree) && (tree->extent/(moodulus(root->solebody.position - tree->cog)) < 0.3 or tree->isleaf))
    {
        long double A = G*tree->cogmass/(moodulus(tree->cog - root->solebody.position)*moodulus(tree->cog - root->solebody.position)*moodulus(tree->cog - root->solebody.position));
        array<long double,3> B = tree->cog - root->solebody.position;
        root->solebody.newacceleration = root->solebody.newacceleration + A*B;
    }
    else
    {
        for(size_t i{0}; i < 8; i++)
        {
            root = updatesingleacceleration(root,tree->Nodelist[i]);
        }
    }
    return(root);
}

/**
 * @brief Updates the leaf nodes of the tree recursively, applying the velocity-verlet algorithm.
 * 
 * @param tree Input tree to update
 * @return Node* 
 */
Node* bodygen::update(Node* tree)
{
    if(tree == NULL)
    {
        return(tree);
    }
    if(tree->isleaf)
    {
        array<long double,3> sumacc = tree->solebody.acceleration + tree->solebody.newacceleration; 

        tree->solebody.position = tree->solebody.position + timestep*tree->solebody.velocity + (0.5*timestep*timestep)*tree->solebody.acceleration;
        tree->solebody.velocity = tree->solebody.velocity + (0.5*timestep)*sumacc;
        tree->solebody.acceleration = tree->solebody.newacceleration;
        tree->solebody.newacceleration = {0,0,0};
        bodyvector[tree->solebody.index] = tree->solebody;
    }
    else
    {
        for(size_t i{0}; i < 8; i++)
        {
            tree->Nodelist[i] = update(tree->Nodelist[i]);
        }
    }
    return(tree);
}


