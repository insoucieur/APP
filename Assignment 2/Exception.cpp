#include <iostream>
#include <fstream>
#include "Exception.hpp"

/*The constructors for each of the two classes should take only the probString argument
and should set the tagString member to an appropriate string.
Write a catch block which is able to catch a generic exception
but can also differentiate between these two types of error.*/

//Constructor
Exception::Exception(std::string tagString, std::string probString)
{
    mTag = tagString;
    mProblem = probString;
}

void Exception::PrintDebug () const
{
    std::cerr << "** Error ("<<mTag<<") ** \n";
    std::cerr << "Problem: " << mProblem << "\n\n";
}
