#include <iostream>
#include <fstream>
#include "FileNotOpenException.hpp"

FileNotOpenException::FileNotOpenException(std::string probString):Exception("FILE NOT OPEN", probString)
{

}

