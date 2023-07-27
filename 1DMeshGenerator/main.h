#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <cmath>
using namespace std;
#define ThreeD 1

  #ifdef ThreeD
  #define DIM 3
     int GLO_Nodes[3] = {101,2,2};
     double GLO_MeshMin[3] = {0.0,0.0,0.0};
     double GLO_MeshMax[3] = {80.0,0.05,0.05};
  #elif TwoD 
  #define DIM 2
     int GLO_Nodes[2] = {2,2};
  #endif

#include "classes.h"

