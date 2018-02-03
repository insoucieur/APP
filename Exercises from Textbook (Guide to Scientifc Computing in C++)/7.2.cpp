#include <iostream>

class AbstractPerson
{
public:
    void virtual Print() = 0; //preferred method of making an abstract class
    //with pure virtual method (so cannot be instantiated)
};
class Mother: public AbstractPerson
{
public:
    void Print(){std::cout<<"Mother\n";}
};
class Daughter: public Mother
{
public:
    void Print(){std::cout<<"Daughter\n";}
};

int main(int argc, char* argv[])
{
    AbstractPerson* p_mother = new Mother;
    AbstractPerson* p_daughter = new Daughter;
    p_mother -> Print();
    p_daughter -> Print();
    delete p_mother;
    delete p_daughter;
    /*//4.
    AbstractPerson* p_abstract = new AbstractPerson;
    p_abstract->Print();
    delete p_abstract;*/
}
//2. remove public from inheritance declaration of derived class (lines 9 and 15)
//then get error: cannot cast derived class to its private base class

//3. remove virtual keyword from line 6:
    //then just get the first fn in the private class so just get
    //Never instantiate
    //Never instantiate
//remove virutal keyword from line 12 or adding virtual on line 18
    //still get the same output as originally

//4. just get Never instantiate printed as well

//6. remove virtual keyword from line 6:
    //error: 'Print' is not virtual and cannot be declared pure
//remove virtual keywrod from line 12:
    //still get original output

//7. doesn't work
