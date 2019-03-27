#ifndef WRITE_BINVOX_H
#define WRITE_BINVOX_H

#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <cstdint>

// Write a binvox file format
//
// Inputs:
//	 S      side(0)*side(1)*side(2) vector of binary values
//   side   3 vector with dimensions of the voxel grid
inline void writeBINVOX(
    std::string filename,
	Eigen::VectorXf &S,
	Eigen::Vector3i &side,
    const int isovalue);

inline void writeBINVOX(
    std::string filename,
	Eigen::VectorXf &S,
	Eigen::Vector3i &side,
    const int isovalue)
{
    int width = side(0);
    int height = side(1);
    int depth = side(2);
    int largest = std::max(std::max(width, height),depth);

    std::cout << "writing file" << std::endl;
    std::ofstream binvox_file;
    binvox_file.open(filename, std::ios::out | std::ios::binary);
    

    binvox_file << "#binvox 1" << std::endl;
    binvox_file << "dim " << side(0) << " " 
        << side(1) << " "
        << side(2) << std::endl;
    binvox_file << "translate " << 0.0f << " "
        << 0.0f << " "
        << 0.0f << std::endl;
    binvox_file << "scale " << 1.0f << std::endl;
    binvox_file << "data" << std::endl;

    S = (S.array() >= isovalue).select(1, S);

    std::cout << "wrote header" << std::endl;

    // std::vector<uint8_t> vec;
    // std::cout << vec.size() << std::endl;
    // for(int i = 0; i < largest*largest*largest; i++)
    // {
    //     vec.push_back(0);
    //     vec.push_back(1);
    // }
    // std::cout << "init array" << std::endl;


    // keep a sort of state machine for writing run length encoding
    uint8_t state = S.row(0).value();
    uint8_t ctr = 0;


    for(int x = 0; x < width; x++)
    {
        for(int y = 0; y < height; y++)
        {
            for(int z = 0; z < depth; z++)
            {
                auto get_index = [&](int x, int y, int z) 
                { 
                    int index = (x * width*height) + (z * width) + y;  
                    return index; 
                };
                int index = get_index(x, y, z);
                std::cout << index << " of " << S.rows() << std::endl;
                float visibility = S.row(index).value();
                uint8_t v = (uint8_t)visibility;
                
                if(v == state)
                {
                    ctr++;
                    // if ctr hits max, dump
                    if(ctr == 255)
                    {
                        binvox_file.write(reinterpret_cast<const char *>(&state), sizeof(state));
                        binvox_file.write(reinterpret_cast<const char *>(&ctr), sizeof(ctr));
                        // std::cout << state << " " << ctr << std::endl;
                        ctr = 0;
                    }
                }
                else
                {
                    // if switch state, dump
                    binvox_file.write(reinterpret_cast<const char *>(&state), sizeof(state));
                    binvox_file.write(reinterpret_cast<const char *>(&ctr), sizeof(ctr));
                    // std::cout << state << " " << ctr << std::endl;

                    state = v;
                    ctr = 1;
                }
            }
        }
    }

            // for(int index = 0; index < S.rows(); index++)
            // {
                // float visibility = S.row(index).value();
                // uint8_t v = (uint8_t)visibility;
                
                // if(v == state)
                // {
                //     ctr++;
                //     // if ctr hits max, dump
                //     if(ctr == 255)
                //     {
                //         binvox_file.write(reinterpret_cast<const char *>(&state), sizeof(state));
                //         binvox_file.write(reinterpret_cast<const char *>(&ctr), sizeof(ctr));
                //         std::cout << state << " " << ctr << std::endl;
                //         ctr = 0;
                //     }
                // }
                // else
                // {
                //     // if switch state, dump
                //     binvox_file.write(reinterpret_cast<const char *>(&state), sizeof(state));
                //     binvox_file.write(reinterpret_cast<const char *>(&ctr), sizeof(ctr));
                //     std::cout << state << " " << ctr << std::endl;

                //     state = v;
                //     ctr = 1;
                // }
    // }
    
    // uint8_t arr[4] = {vec[0], vec[1], vec[2], vec[3]};
    // std::cout << (char*)&arr << std::endl;
    // flush out remainders
    if (ctr > 0)
    {
        binvox_file.write(reinterpret_cast<const char *>(&state), sizeof(state));
        binvox_file.write(reinterpret_cast<const char *>(&ctr), sizeof(ctr));
    }
    // std::cout << vec.size() << std::endl;
    // binvox_file.write((char*)&vec[0], vec.size());
    binvox_file.close();
}

#endif