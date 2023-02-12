#ifndef __EXPORT_HPP__
#define __EXPORT_HPP__

#include <iostream>
#include <fstream>

class Exporter
{
    public:
        Exporter(std::string afilename);
        ~Exporter();
        
        void operator << (std::string arg);

    private:
        std::string filename { };
        std::ofstream outfile;

};

#endif  // export.hpp
