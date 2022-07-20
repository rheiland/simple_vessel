#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdio>

int main()
{
    std::vector<double> position = {1.,2.,3};

    FILE* fp; 
    std::string filename = "test.mat";
    fp = fopen( filename.c_str() , "wb" );
    // fp = fopen( "test.mat" , "wb" );
    if( fp == NULL )
    {
        std::cout << "Error: could not open file " << filename << "!" << std::endl;
        return 0;
    }

    // fwrite( (char*) &( position ) , sizeof(double) , 3 , fp );
    std::fwrite( position.data(), sizeof(position[0]) , 3 , fp );
    return 0;
}