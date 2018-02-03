#include <vector>
#include <algorithm>
#include "Vector_old.hpp"

// Computes the distances of all points in the range [begin,end) to a given query point.
//
// The results are stored in the Range starting with distanceBegin having the same length
// as the range of points. The i-th element of the range stores the distance of the i-th point
// in the range to the query point as a double value.
//
// The supplied iterators can be assumed to be bidirectional
template<class SetIterator, class DistanceIterator>
void computeDistances(
	SetIterator begin, SetIterator end,//pointers to the start and end of the set of points
	Vector const& query, //the query point
	DistanceIterator distanceBegin //iterator to the beginning of the range storing the distances
){
    while(begin != end) //iterate through set of points and calculate the norm of difference between that point and query point
        {*distanceBegin++ = (*begin++ - query).CalculateNorm();} //distanceBegin is the output iterator to store the results

}

// class representing a Pair of a distance and an iterator to the point
// with the distance to a query point.
template<class Iterator>
struct Pair{
	double distance;
	Iterator iterator;
};

// Computes the distances of all points in the range [begin,end) to a given query point
// and returns a sorted array of pairs: the distance to the query point and the iterator
// in the range that has this distance.
//
// Elements are returned in ascending order, i.e. the pair with the smallest distance
// is at the front, next is the pair with the second smallest distance, etc.
template<class SetIterator>
std::vector<Pair<SetIterator> > createSortedPointDistancePairs(
	SetIterator begin, SetIterator end,//pointers to the start and end of the set of points
	Vector const& query //the query point
){
    int l = std::distance(begin, end); //get length of SetIterator
    std::vector<double> distances(l); //make a vector called distances same length as SetIterator to store the results
    computeDistances(begin, end, query, distances.begin()); //calculate distances using above function and store results in distances vector

    std::vector<Pair<SetIterator> > v(l); //to access members in struct Pair and declare an std::vector v that is the length of SetIterator
    for (int i = 0; i < l; i++) //iterate through the length of the vector
        v[i] = Pair<SetIterator>{distances[i], begin+i;} //initialises a Pair struct of distance and iterator pointer, which in this case is SetIterator with begin and end
    //begin+i is similar to begin[i]
    //v.begin() returns iterator to start of vector v
    //v.end() returns iterator to end of vector v
        //comes with the containter class std::vector
    std::sort(v.begin(), v.end(), [](Pair<SetIterator> a, Pair<SetIterator> b){return a.distance < b.distance;});
    //third parameter is a lambda function to compare 2 pairs since there is not a native way to compare 2 Pair structs
    //third is lambda expression and return a < b gives ascending order
    //(can also use bool cmp(const Type1 &a, const Type2 &b) which returns True if first argument is less than 2nd argmenet
    return v;
}

//Returns the Iterators to the k nearest neighbours of a given query point.
//
// The iterators are sorted in ascending order by distance to the query point.
//
// The set of points is given in the range [begin,end].
template<class SetIterator>
std::vector<SetIterator> kNearestNeighbour(
	SetIterator begin, SetIterator end,//pointers to the start and end of the set of points
	Vector const& query, //the query point
	unsigned int k // the k nearest neighbours that are going to be returned
	//unsigned so only nonnegative integers
){
    auto all = createSortedPointDistancePairs(begin, end, query); //auto => any type and all stores all the sorted iterators in ascending order
    //declare std::vector that is of type SetIterator so k_nearest is a vector of iterators
    std::vector<SetIterator> k_nearest(k); //new vector that is same length as k to store results of k nearest neighbours
    for (unsigned i = 0; i < k; i++)
        k_nearest[i] = all[i].iterator; //takes all iterators from sorted pairs and puts them in a vector of iterators
    return k_nearest;
}

