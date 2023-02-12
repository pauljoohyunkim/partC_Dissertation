#include "export.hpp"
#include <iostream>
#include <fstream>

Exporter::Exporter(std::string afilename)
{
    filename = afilename;
    outfile.open(filename);
}

Exporter::~Exporter()
{
    outfile.close();
}

void Exporter::operator << (std::string arg)
{
    outfile << arg;
}
