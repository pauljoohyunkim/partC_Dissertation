#ifndef __EXPORT_CUH__
#define __EXPORT_CUH__

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

#endif  // export.cuh
