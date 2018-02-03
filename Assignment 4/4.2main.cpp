#include "Vector_old.hpp"
#include "assignment4-2.hpp"
#include <iostream>

int main()
{
    std::vector<Vector> v{Vector(3), Vector(3), Vector(3), Vector(3)};
    for (int i = 0; i < 4; i++)
    {
        v[i] = Vector(3);
        v[i][0] = i;
        v[i][1] = i+1;
        v[i][2] = i*2;
        std::cout << "v[" << i << "] = {" << i << " " << i+1 << " " << i*2 << "}" << std::endl;
    }

    Vector vec(3);
    vec[0] = 10;
    vec[1] = 5;
    vec[2] = 7;
    std::cout << "v[] = {" << vec[0] << " " << vec[1] << " " << vec[2] << "}" << std::endl;
    auto n = kNearestNeighbour<decltype(v.begin())>(v.begin(), v.end(), vec, 2);
    for (auto i : n)
        std::cout << "v[] = {" << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << "}" << std::endl;
    return 0;
}

