#include <iostream>
#include <fstream>
#include "OutOfRangeException.hpp"

OutOfRangeException::OutOfRangeException(std::string probString):Exception("OUT OF RANGE", probString)
{

}

