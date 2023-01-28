#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "export.hpp"

void vectorParse(std::ostream& stream, std::vector<double>& vec, unsigned int J)
{
    std::stringstream strStream;

    /* Opening Bracket */
    strStream << "[" << std::setprecision(DOUBLEPRECISION) << std::fixed;

    /* Put each value into the stream */
    for (auto i = 0; i < J - 1; i++)
    {
        strStream << vec[i] << ", ";
    }
    /* Closing Bracket */
    strStream << vec[J - 1] << "]";
    std::string parseString = strStream.str();

    stream << parseString;
}
