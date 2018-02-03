#include <vector>
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

    for (SetIterator pos = begin; pos!=end; ++pos, ++distanceBegin)
    {
        *distanceBegin = len(*pos - query);
    }

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
){// access members in struct Pair
    std::vector<Pair<SetIterator> > v;
    for (SetIterator pos = begin; pos!=end; ++pos)
    {
        computeDistances(*pos, end, query, *pos);
        std::sort(v.begin(), v.end());
    }
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
){
    std::vector <SetIterator> k_nearest;
    std::vector <Pair<SetIterator> > v;
    for (SetIterator pos = begin; pos!=end; ++pos)
    {
     v = createSortedPointDistancePairs(pos, end, query);
    }

    for (int i = 0; i < k; i++)
    {
        k_nearest = v[i];
    }
    return k_nearest;
}
