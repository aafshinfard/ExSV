#ifndef EXTRATOOLS_H
#define EXTRATOOLS_H


#include <string>
using namespace std;




string convertNumToStr(long long number) {
    ostringstream ostr;
    ostr<<number;
    return ostr.str();
}




#endif // EXTRATOOLS_H
