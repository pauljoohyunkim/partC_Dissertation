#include "export.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

JsonExporter::JsonExporter(std::string afilename)
{
    filename = afilename;
    outfile.open(filename);
    outfile << "[" << std::setprecision(DOUBLEPRECISION) << std::fixed;
}

JsonExporter::~JsonExporter()
{
    outfile << "]";
    outfile.close();
}

void JsonExporter::operator << (std::string arg)
{
    outfile << arg;
}

void JsonExporter::operator << (std::vector<double>& vec)
{
    auto len = vec.size();

    if (!firstEntry)
    {
        outfile << ",\n";
    }

    outfile << "[";
    bool firstValueOfVector { true };
    for (auto val : vec)
    {
        if (!firstValueOfVector)
        {
            outfile << ",";
        }
        outfile << val;
        firstValueOfVector = false;
    }
    outfile << "]";
    
    firstEntry = false;
}
