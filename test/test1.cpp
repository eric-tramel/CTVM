#include <iostream>
#include "ctvm.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

int main()
{
    using namespace boost::numeric::ublas;
    vector<double> v (10);

    // Test Library Linking
    foo();

    for (unsigned i = 0; i < v.size (); ++ i)
    v (i) = i;

    std::cout << v << std::endl;

    //Using functions of vector
    v.resize(15,1);
    std::cout << "nsize increased to 15 "<< v << std::endl;

    v.resize(10,1);
    std::cout << v << std::endl;

    v.insert_element(0,10);
    std::cout << "nInserted 10 " << v << std::endl;

    v.erase_element(0);
    std::cout << "nRemoving " << v << std::endl;

    v.clear();
    std::cout << "nClearing " << v << std::endl;
}