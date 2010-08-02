#include "periodic_table.h"

namespace fluidPhysics {
  struct element {
    char name[3] ;
    double mass ;
  } ;

  element pt[] = {
    {"e",    0.00055},  // electron

    {"H" ,   1.00797},

    {"He",   4.003},

    {"Li",   6.939},
    {"Be",   9.012},

    {"B",   10.811},
    {"C",   12.011},
    {"N",   14.007},
    {"O",   15.999},
    {"F",   18.998},
    {"Ne",  20.183},

    {"Na",  22.99},
    {"Mg",  24.31},

    {"Al",  26.981},
    {"Si",  28.09},
    {"P",   30.974},
    {"S",   32.064},
    {"Cl",  35.45},
    {"Ar",  39.948},

    {"K",   39.102},
    {"Ca",  40.08},

    {"Sc",  44.96},
    {"Ti",  47.90},
    {"V",   50.94},
    {"Cr",  51.996},
    {"Mn",  54.94},
    {"Fe",  55.847},
    {"Co",  58.93},
    {"Ni",  58.71},
    {"Cu",  63.54},
    {"Zn",  65.37},
    {"Ga",  69.72},
    {"Ge",  72.59},
    {"As",  74.92},
    {"Se",  78.96},
    {"Br",  79.91},
    {"Kr",  83.80},

    {"Rb",  85.47},
    {"Sr",  87.62},

    {"Y",   88.905},
    {"Zr",  91.22},
    {"Nb",  92.906},
    {"Mo",  95.94},
    {"Tc",  99.},
    {"Ru", 101.07},
    {"Rh", 102.905},
    {"Pd", 106.4},
    {"Ag", 107.87},
    {"Cd", 112.40},
    {"In", 114.82},
    {"Sn", 118.69},
    {"Sb", 121.75},
    {"Te", 127.60},
    {"I",  126.90},
    {"Xe", 131.3},

    {"Cs", 132.91},
    {"Ba", 137.34},

    {"La", 138.91},
    {"Ce", 140.12},
    {"Pr", 140.907},
    {"Nd", 144.24},
    {"Pm", 145.},
    {"Sm", 150.35},
    {"Eu", 151.96},
    {"Gd", 157.25},
    {"Tb", 158.924},
    {"Dy", 162.50},
    {"Ho", 164.93},
    {"Er", 167.26},
    {"Tm", 168.93},
    {"Yb", 173.04},
    {"Lu", 174.97},

    {"Hf", 178.5},
    {"Ta", 180.95},
    {"W",  183.85},
    {"Re", 186.2},
    {"Os", 190.2},
    {"Ir", 192.2},
    {"Pt", 195.09},
    {"Au", 196.97},
    {"Hg", 200.59},
    {"Tl", 204.37},
    {"Pb", 207.19},
    {"Bi", 208.98},
    {"Po", 210.},
    {"At", 210.},
    {"Rn", 222.0},

    {"Fr", 223.},
    {"Ra", 226.05},

    {"Ac", 227.},
    {"Th", 232.04},
    {"Pa", 231.},
    {"U",  238.03},
    {"Np", 237.},
    {"Pu", 242.},
    {"Am", 243.},
    {"Cm", 245.},
    {"Bk", 249.},
    {"Cf", 249.},
    {"Es", 253.},
    {"Fm", 254.},
    {"Md", 256.},
    {"No", 254.},
    {"Lw", 257.},

    {"_" , 0.},    // User defined species (Eg. _Air)

    {"M",  0.},    // Used for M-body equations
  } ;

  periodic_table::periodic_table()
  {
    for(unsigned int i=0;i<sizeof(pt)/sizeof(element);++i) {
      elements[pt[i].name] = i ;
    }
  }

  bool periodic_table::valid_element(string el)
  {
    return elements.find(el) != elements.end() ;
  }

  int periodic_table::find_table_entry(string el)
  {
    if(valid_element(el))
      return elements[el] ;
    else
      return -1 ;
  }

  double periodic_table::molecular_mass(int entry)
  {
    return pt[entry].mass ;
  }    
}
