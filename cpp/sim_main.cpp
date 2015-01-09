#include "sim.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main(int argc, const char* argv[])
{
  sim(atoi(argv[1]), atoi(argv[2]), argv[3]);
  return 0;
}
