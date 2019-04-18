#ifndef WRITE_PGM_H
#define WRITE_PGM_H

#include <Eigen/Core>
#include <iostream>

inline bool write_pgm(
  const std::string & filename,
  const Eigen::VectorXi &data,
  const Eigen::Vector3i &side,
  const int isovalue);

inline bool write_pgm(
  const std::string & filename,
  const Eigen::VectorXi &data,
  const Eigen::Vector3i &side,
  const int isovalue)
{
  Eigen::VectorXi S_iso = (data.array() >= isovalue).select(isovalue, data.array()-data.array());

  std::ofstream s(filename);
  if(!s.is_open())
  {
    std::cerr<<"Failed to open "<<filename<<std::endl;
    return false;
  }
  
  s<<"P5\n"<<
    side(0) <<" "<<side(1)<<" "<<side(2) <<"\n"<<
    S_iso.maxCoeff()<<"\n";

  for(int index = 0; index < S_iso.rows(); index ++)
  {
    int val = S_iso.row(index).value();
    s<< val <<" ";
  }

  s<<"\n";

  return true;
}


#endif