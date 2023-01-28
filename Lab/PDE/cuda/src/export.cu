#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "export.hpp"

void vectorParse(std::ostream& stream, std::vector<double>& vec, unsigned int J)
{
    std::stringstream strStream;

    /* Opening Bracket */
    strStream << "[";

    /* Put each value into the stream */
    for (auto i = 0; i < J - 1; i++)
    {
        strStream << std::fixed << std::setprecison(DOUBLEPRECISION) << vec[i] << ", ";
    }
    /* Closing Bracket */
    strStream << std::fixed << std::setprecison(DOUBLEPRECISION) << vec[J - 1] << "]";
    std::string parseString = strStream.str();

    stream << parseString;
}
