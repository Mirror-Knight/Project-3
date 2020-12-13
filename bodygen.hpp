#include <iostream>
#include <string>
#include <array>
#include <vector>

using namespace std;

/**
 * @brief Body class that stores the information of a body - position, velocity, acceleration, mass, and radius. Index keeps track of the body in updating "bodyvector" (in another class), and newacceleration is needed in Velocity-Verlet
 * 
 */
class body
   {
        public:
            array<long double, 3> position, velocity, acceleration, newacceleration;
            long double mass, radius;
            int index;
    };

/**
 * @brief Region class stores boundaries of a region, a vector of bodies in the region, a nodepath and a boolean that checks if collisions have been computed
 * @param regnodepath This is a string that encodes the path one takes to reach some child node from the root node of the tree.
 * @param checkcol This checks if collision has been computed. This is set to false initially, and set to true if collisions are checked, preventing needless extra computations at child nodes.
 */
class region
    {   
        public:
            array<long double, 2> xrange, yrange, zrange;
            vector<body> bodiesinregion;
            string regnodepath;
            bool checkcol;
    };

/**
 * @brief Node struct that stores various elements - main important detail is that this Node points to eight other node pointers.
 * @param isleaf Is this node a leaf
 * @param nodepath Same as regnodepath from region. The way my code works requires this to exist in both the region class and the node struct
 * @param cog This is the center of gravity coordinates
 * @param cogmass Center of mass of the node
 * @param extent Approximate "size" or "spread" of bodies in a node
 * @param solebody The body stored in a node if it is a leaf. Otherwise this stores garbage data that will not be referenced
 * 
 */
struct Node
{
    bool isleaf;
    string nodepath;
    array<long double,3> cog;
    long double cogmass;
    long double extent;
    body solebody;
    array<Node*,8> Nodelist; 
};

/**
 * @brief A class that basically makes the tree
 * 
 */
class Spacetree
{
    public:
        Spacetree(region);
        Node* treegen();
    private:
        region regi;
        Node* addnulls(Node*);
        Node* makeatree(region);
        vector<body> updatecollision(vector<body> &);
        vector<body> mergebodies(vector<body> &);
};

/**
 * @brief The main class that runs the simulation or builds the bodies. Most paramters are straightforward
 * @param bodyvector This vector stores the information of all bodies in the simulation. When each leaf is updated, the corresponding index in this bodyvector is also updated
 * 
 */
class bodygen
{
    private:
        Node* update(Node*);
        Node* makebodies();
        Node* updatesingleacceleration(Node*, Node*);
        Node* updateallacceleration(Node*, Node*);
        void deletetree(Node*);

        bool comparetree(Node*, Node*);
        array<long double,6> calcminmax();
        array<long double,2> randcircgen(long double, long double);
        array<long double,3> randspheregen(long double, long double);
        
        size_t count{100};
        string filename;
        long double timestep{1};
        size_t iterations{100};

        Node* datatree;
        region space;
        bool writeinitfile;
        vector<body> bodyvector;
    public:
        bodygen(string, long double, size_t);
        bodygen(size_t, string, long double, size_t);
        void simulate();
};






