#ifndef __EXPORT_HPP__
#define __EXPORT_HPP__

#include <iostream>
#include <fstream>
#include <vector>

#define DOUBLEPRECISION 15

class JsonExporter
{
    public:
        JsonExporter(std::string afilename);
        ~JsonExporter();
        
        void operator << (std::string arg);
        void operator << (std::vector<double>& vec);

    private:
        std::string filename { };
        std::ofstream outfile;
        bool firstEntry { true };

};

#endif  // export.hpp
